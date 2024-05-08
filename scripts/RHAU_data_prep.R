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
# RHAU data formatting
file_list <- list.files(path = here("data", "RHAU_data", "RHAU_tags", "Processed_data", "CSV")) #get the list of files in the directory

dfs <- list() #create empty list to store each data frame

ids <- c(45653, 45659, 45701, 44149, 45663, 45672, 45695) #list of RHAU tag IDs
nest_x <- as.numeric(rep(48.127757, length(ids)))
nest_y <- as.numeric(rep(-122.922553, length(ids))) 
nest_locs <- data.frame(ID = ids,
                        nest_x = nest_x, 
                        nest_y = nest_y)

# Loop through each file and read it, skipping the first 5 lines of data frame title
for (i in seq_along(file_list)) {
  df <- read.csv(file.path(here("data", "RHAU_data", "RHAU_tags", "Processed_data", "CSV"), file_list[i]), header = FALSE, skip = 5) #read in each CSV
  df <- mutate(df, ID = ids[i]) #add IDs to each data frame
  df <- filter(df, df[[9]] != 0.00000)  #remove rows where position is not reported
  dfs[[i]] <- df #store data frame in list
}

# Combine all data frames into a single data frame
RHAU_data <- bind_rows(dfs)
colnames(RHAU_data) <- c("DAY", "MONTH", "YEAR", "HOUR", "MINUTE", "SECOND", "SoD", "#SAT", "LAT", "LONG", "ALT", "CLOCK_OFFSET", "ACCURACY", "BATTERY", "NA1", "NA2", "ID")

# Create datetime string
RHAU_data_datetime <- lapply(RHAU_data[,1:6], as.character) #convert date/time columns to character type

datetime_string <- with(RHAU_data_datetime, sprintf("%s-%s-%s %s:%s:%s", YEAR, MONTH, DAY, HOUR, MINUTE, SECOND))

# Finalize proto RHAU data
RHAU_data <- data.frame(
  ID = as.character(RHAU_data[,17]),  
  x = RHAU_data[,9],  
  y = RHAU_data[,10],  
  time = as.POSIXct(datetime_string, format = "%Y-%m-%d %H:%M:%S"
  ))

# Let's reorganize data
RHAU_data <- RHAU_data %>%
  group_by(ID) %>%
  arrange(time, .by_group = TRUE)

# Save the latlon data!
write.table(RHAU_data, file = here("data", "RHAU_data", "RHAU_data_latlon.csv"), col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
###############################
### PROJECT LOC DATA TO UTM ###
###############################

utm_proj <- "+proj=utm +zone=10 +north +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" #define the PROJ string

# Convert x and y col to spdf
RHAU_data_sp <- SpatialPointsDataFrame(coords = RHAU_data[, c(3, 2)],
                                       data = RHAU_data,
                                       proj4string = CRS("+proj=longlat +datum=WGS84"))

nest_data_sp <- SpatialPointsDataFrame(coords = nest_locs[, c(3, 2)],
                                       data = nest_locs,
                                       proj4string = CRS("+proj=longlat +datum=WGS84"))

RHAU_data_utm <- spTransform(RHAU_data_sp, CRS(utm_proj)) #transform the coordinates to UTM zone 10N; around Protection Island
nest_data_utm <- spTransform(nest_data_sp, CRS(utm_proj)) #same thing

# Make nest_data UTM dataframe 
nest_data <- data.frame(
  ID <- nest_locs[,1],
  nest_x = nest_data_utm@coords[,1],
  nest_y = nest_data_utm@coords[,2]
)

colnames(nest_data) <- c("ID", "nest_x", "nest_y")

# Finalize RHAU dataframe design
RHAU_data <- data.frame(
  ID = RHAU_data[,1],  
  x = RHAU_data_utm@coords[,1],  
  y = RHAU_data_utm@coords[,2],  
  time = RHAU_data[,4]
)

RHAU_nest_data <- merge(RHAU_data, nest_data, by = "ID", all.x = TRUE)

# Let's reorganize data
RHAU_data <- RHAU_nest_data %>%
  group_by(ID) %>%
  arrange(time, .by_group = TRUE)

# Save the data!
write.table(RHAU_data, file = here("data", "RHAU_data", "RHAU_data_UTM.csv"), col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)
