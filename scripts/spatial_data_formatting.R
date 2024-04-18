install.packages("sp")
install.packages("rgdal")
library(sp)
library(rgdal)
library(here)
library(dplyr)


############################
### PIGU data formatting ###
############################
file_list <- list.files(path = here("data", "PIGU_data", "PIGU_tags", "Processed_data", "CSV")) #get the list of files in the directory

dfs <- list() #create empty list to store each data frame

ids <- c(44067, 45657, 45658, 44072, 44372, 44505) #list of PIGU tag IDs

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

# Project lat/long as UTM
utm_proj <- "+proj=utm +zone=10 +north +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs" #define the PROJ string

# Convert columns 9 and 10 to spatial points data frame
PIGU_data_sp <- SpatialPointsDataFrame(coords = PIGU_data[, c(10, 9)],
                                       data = PIGU_data,
                                       proj4string = CRS("+proj=longlat +datum=WGS84"))

PIGU_data_utm <- spTransform(PIGU_data_sp, CRS(utm_proj)) #transform the coordinates to UTM zone 10N; around Protection Island

# Create datetime string
PIGU_data_datetime <- lapply(PIGU_data[,1:6], as.character) #convert date/time columns to character type

datetime_string <- with(PIGU_data_datetime, sprintf("%s-%s-%s %s:%s:%s", YEAR, MONTH, DAY, HOUR, MINUTE, SECOND))

# Finalize PIGU data
PIGU_data <- data.frame(
  ID = PIGU_data[,17],  
  x = PIGU_data_utm@coords[,1],  
  y = PIGU_data_utm@coords[,2],  
  time = datetime_string
)





# Group the data by ID
grouped_data <- group_by(PIGU_data, ID)

# Initialize an empty list to store plots
plot_list <- list()

# Loop through each unique ID
for (id in unique(PIGU_data$ID)) {
  # Subset data for the current ID
  subset_data <- filter(PIGU_data, ID == id)
  
  # Create plot for the current ID
  p <- ggplot(data = subset_data, aes(x = x, y = y)) +
    geom_point() +
    labs(title = paste("Track for ID:", id),
         subtitle = paste("ID:", id)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 10, face = "bold"),
          plot.margin = margin(t = 30))
  
  # Add plot to the list
  plot_list[[id]] <- p
}

# Arrange plots in a grid
grid.arrange(grobs = plot_list, ncol = 3)  # Change ncol as needed







cat(paste0('"', paste(colnames(PIGU_data), collapse = '","'), '"'), file = here("data", "PIGU_data", "PIGU_data.csv"), sep = '\n')
write.table(PIGU_data, file = here("data", "PIGU_data", "PIGU_data.csv"), col.names = FALSE, row.names = FALSE, sep = ",", quote = FALSE, append = TRUE)
############################
### RHAU data formatting ###
############################
file_list <- list.files(path = here("data", "RHAU_data", "RHAU_tags", "Processed_data", "CSV")) #get the list of files in the directory

dfs <- list() #create empty list to store each data frame

ids <- c(44149, 44167, 44369, 44469, 44720, 45653, 45659, 45663, 45672, 45695, 45701) #list of RHAU tag IDs

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

# Convert columns 9 and 10 to spatial points data frame
RHAU_data_sp <- SpatialPointsDataFrame(coords = RHAU_data[, c(10, 9)],
                                       data = RHAU_data,
                                       proj4string = CRS("+proj=longlat +datum=WGS84"))

RHAU_data_utm <- spTransform(RHAU_data_sp, CRS(utm_proj)) #transform the coordinates to UTM zone 10N; around Protection Island

# Create datetime string
RHAU_data_datetime <- lapply(RHAU_data[,1:6], as.character) #convert date/time columns to character type

datetime_string <- with(RHAU_data_datetime, sprintf("%s-%s-%s %s:%s:%s", YEAR, MONTH, DAY, HOUR, MINUTE, SECOND))

# Finalize RHAU data
RHAU_data <- data.frame(
  ID = RHAU_data[,17],  
  x = RHAU_data_utm@coords[,1],  
  y = RHAU_data_utm@coords[,2],  
  time = datetime_string
)

cat(paste0('"', paste(colnames(RHAU_data), collapse = '","'), '"'), file = here("data", "RHAU_data", "RHAU_data.csv"), sep = '\n')
write.table(RHAU_data, file = here("data", "RHAU_data", "RHAU_data.csv"), col.names = FALSE, row.names = FALSE, sep = ",", quote = FALSE, append = TRUE)

