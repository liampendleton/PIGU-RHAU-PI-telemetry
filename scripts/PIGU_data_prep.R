library(sp)
library(sf)
library(here)
library(dplyr)
library(gridExtra)
library(grid)
library(ggplot2)

###################################
######### DATA FORMATTING #########
###################################
# PIGU data formatting
file_list <- list.files(path = here("data", "PIGU_data", "PIGU_tags", "Processed_data", "CSV")) #get the list of files in the directory

dfs <- list() #create empty list to store each data frame

ids <- c(44067, 45657, 45658, 44072, 44372, 44505) #list of PIGU tag IDs
nest_ID <- c("MAR_3", "MAR_11", "MS_14", "MAR_2",  "MAR_10", "MAR_1")
nest_x <- as.numeric(c(48.127846, 48.12736, 48.124139, 48.127946, 48.127453, 48.128024))
nest_y <- as.numeric(c(-122.920373, -122.920626, -122.927807, -122.920335, -122.920551, -122.920301)) 
nest_locs <- data.frame(ID = ids,
                        nest_ID = nest_ID,
                        nest_x = nest_x, 
                        nest_y = nest_y)

# Loop through each file and read it, skipping the first 5 lines of data frame title
for (i in seq_along(file_list)) {
  df <- read.csv(file.path(here("data", "PIGU_data", "PIGU_tags", "Processed_data", "CSV"), file_list[i]), header = FALSE, skip = 5) #read in each CSV
  df <- mutate(df, ID = ids[i]) #add IDs to each data frame
  df <- filter(df, df[[9]] != 0.00000)  #remove rows where position is not reported
  dfs[[i]] <- df #store data frame in list
}

# Combine all data frames into a single data frame
PIGU_data <- bind_rows(dfs)
colnames(PIGU_data) <- c("DAY", "MONTH", "YEAR", "HOUR", "MINUTE", "SECOND", "SoD", "#SAT", "LAT", "LONG", "ALT", "CLOCK_OFFSET", "ACCURACY", "BATTERY", "NA1", "NA2", "ID")

# Create datetime string
PIGU_data_datetime <- lapply(PIGU_data[,1:6], as.character) #convert date/time columns to character type

datetime_string <- with(PIGU_data_datetime, sprintf("%s-%s-%s %s:%s:%s", YEAR, MONTH, DAY, HOUR, MINUTE, SECOND))

# Finalize proto PIGU data
PIGU_data <- data.frame(
  ID = as.character(PIGU_data[,17]),  
  x = PIGU_data[,9],  
  y = PIGU_data[,10],  
  time = as.POSIXct(datetime_string, format = "%Y-%m-%d %H:%M:%S"
  ))

#####################################################
### VARIOUS WAYS TO VALIDATE DATA; CHECK FOR GAPS ###
#####################################################
# Report indices of gaps in data
data_gaps <- function(data, threshold) {
  gaps <- c()  #empty vec to store gap indices
  for (i in 1:(nrow(data) - 1)) {  #iterate over each pair of consecutive points
    time_diff <- difftime(data$time[i + 1], data$time[i], units = "mins")  #diff between consecutive points
    if (time_diff > threshold) {  #if time diff exceeds the threshold,
      gaps <- c(gaps, i)  #add index of the gap to the vec
    }
  }
  return(gaps)
}

threshold <- 19  #threshold of minutes (19 is probably sweet spot to account for variation in 15 min interval)

PIGU_gap_indices <- data_gaps(PIGU_data, threshold) #use function
print(PIGU_gap_indices) #show where gaps exist so we can manually look later and validate how gaps were addressed

# Report mean and median amount of time between 
time_intervals <- function(data) {
  sorted <- data[order(data$ID, data$time), ] #sort first
  time_diffs <- diff(sorted$time) #get diffs
  time_diffs_minutes <- as.numeric(time_diffs, units = "mins") #convert to minutes
  
  avg_interval <- mean(time_diffs_minutes, na.rm = TRUE)
  median_interval <- median(time_diffs_minutes, na.rm = TRUE)
  
  return(list(average_interval = avg_interval, median_interval = median_interval))
}

PIGU_intervals <- time_intervals(PIGU_data)
print(PIGU_intervals) #avg should be higher than 15; median should be around 15

##########################
### DATA INTERPOLATION ###
##########################
# Make function
interpolate_points <- function(gap_start, gap_end, num_points) {
  delta_x <- (gap_end$x - gap_start$x) / (num_points + 1)  #change in x-coord per interpolated point
  delta_y <- (gap_end$y - gap_start$y) / (num_points + 1)  #change in y-coord per interpolated point
  
  interp_points <- list() #to store interpolated points
  
  # Generate interpolated points
  for (i in 1:num_points) {  
    interpolated_x <- gap_start$x + i * delta_x  #calc x-coord of interpolated point
    interpolated_y <- gap_start$y + i * delta_y  #calc y-coord of interpolated point
    interpolated_time <- gap_start$time + as.difftime(i * 15, units = "mins")  #calc timestamp of interpolated point
    
    # Structure like PIGU dataset!
    interp_points[[i + 1]] <- data.frame(  #use i + 1 as index to start from 1
      ID = gap_start$ID,
      x = interpolated_x,
      y = interpolated_y,
      time = interpolated_time
    )
  }
  return(interp_points)
}


# Format dataframe like original dataset
interp_data <- PIGU_data[0,]

# Set up data and feed into function
for (bird_id in unique(PIGU_data$ID)) {  
  bird_data <- PIGU_data[PIGU_data$ID == bird_id,]  #start with unique ID
  bird_data <- bird_data[order(bird_data$ID, bird_data$time),]  #order by ID and timestamp
  
  for (i in 1:(nrow(bird_data) - 1)) {  #iterate over each pair of consecutive points
    time_diff <- difftime(bird_data$time[i + 1], bird_data$time[i], units = "mins")  #time diff between points
    
    if (time_diff > 19) {  #19 seems like the sweet spot, since not every point was collected at EXACTLY 15 min intervals
      num_point_insert <- as.integer(time_diff / 15) - 1  #calc number of points to interpolate
      
      # Use the function!
      interp_points <- interpolate_points(
        gap_start = bird_data[i,],
        gap_end = bird_data[i + 1,],
        num_points = num_point_insert)
      
      for (j in 1:length(interp_points)) {
        interp_data <- rbind(interp_data, interp_points[[j]])  #add interpolated points to dataframe
      }
    }
  }
}

merged_data <- rbind(PIGU_data, interp_data) #append data

###############################
### PROJECT LOC DATA TO UTM ###
###############################

utm_proj <- "+proj=utm +zone=10 +north +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs" #define the PROJ string

# Convert x and y col to spdf
PIGU_data_sp <- SpatialPointsDataFrame(coords = merged_data[, c(3, 2)],
                                       data = merged_data,
                                       proj4string = CRS("+proj=longlat +datum=WGS84"))

nest_data_sp <- SpatialPointsDataFrame(coords = nest_locs[, c(4, 3)],
                                       data = nest_locs,
                                       proj4string = CRS("+proj=longlat +datum=WGS84"))

PIGU_data_utm <- spTransform(PIGU_data_sp, CRS(utm_proj)) #transform the coordinates to UTM zone 10N; around Protection Island
nest_data_utm <- spTransform(nest_data_sp, CRS(utm_proj)) #same thing

# Make nest_data UTM dataframe 
nest_data <- data.frame(
  ID <- nest_locs[,1],
  nest_x = nest_data_utm@coords[,1],
  nest_y = nest_data_utm@coords[,2]
)

colnames(nest_data) <- c("ID", "nest_x", "nest_y")

# Finalize PIGU dataframe design
PIGU_data <- data.frame(
  ID = merged_data[,1],  
  x = PIGU_data_utm@coords[,1],  
  y = PIGU_data_utm@coords[,2],  
  time = merged_data[,4]
)

PIGU_data$time <- PIGU_data$time + years(2000) #fix 

PIGU_nest_data <- merge(PIGU_data, nest_data, by = "ID", all.x = TRUE)

# Let's reorganize data
PIGU_data <- PIGU_nest_data %>%
  arrange(ID, time) %>%
  group_by(ID)

# Save the data!
write.table(PIGU_data, file = here("data", "PIGU_data", "PIGU_data.csv"), col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
