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

###############################
### PROJECT LOC DATA TO UTM ###
###############################

utm_proj <- "+proj=utm +zone=10 +north +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs" #define the PROJ string

# Convert x and y col to spdf
PIGU_data_sp <- SpatialPointsDataFrame(coords = PIGU_data[, c(3, 2)],
                                       data = PIGU_data,
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
  ID = PIGU_data[,1],  
  x = PIGU_data_utm@coords[,1],  
  y = PIGU_data_utm@coords[,2],  
  time = PIGU_data[,4]
)

PIGU_nest_data <- merge(PIGU_data, nest_data, by = "ID", all.x = TRUE)

# Let's reorganize data
PIGU_data <- PIGU_nest_data %>%
  group_by(ID) %>%
  arrange(ID, time)

# Save the data!
write.table(PIGU_data, file = here("data", "PIGU_data", "PIGU_data_UTM.csv"), col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)

