library(sp)
library(sf)
library(here)
library(dplyr)
library(gridExtra)
library(grid)
library(ggplot2)

##################################################
### Get map data
state <- st_read(here("data", "states", "cb_2018_us_state_500k.shp")) #for some reason I had to use this read command. there is another, read_sf() that you may need to use in different cases 
wa <- state[state$NAME == "Washington", ] #filter to obtain WA state

# Let's focus in on region of interest
lon_min <- -122.96
lon_max <- -122.9
lat_min <- 48.1
lat_max <- 48.15
##################################################
### PIGU data formatting
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

# Create datetime string
PIGU_data_datetime <- lapply(PIGU_data[,1:6], as.character) #convert date/time columns to character type

datetime_string <- with(PIGU_data_datetime, sprintf("%s-%s-%s %s:%s:%s", YEAR, MONTH, DAY, HOUR, MINUTE, SECOND))

# Finalize PIGU data
PIGU_data <- data.frame(
  ID = as.character(PIGU_data[,17]),  
  x = PIGU_data[,9],  
  y = PIGU_data[,10],  
  time = as.POSIXct(datetime_string, format = "%Y-%m-%d %H:%M:%S"
))
##################################################
bird.list <- unique(PIGU_data$ID)

# Iterate over each unique bird ID
for (i in 1:length(bird.list)) {
  bird.data <- PIGU_data[which(PIGU_data$ID == bird.list[i]),] #go through each tracked individual's data
  bird.data <- bird.data[order(bird.data$time),] #order by time
  diff <- bird.data$time[2:nrow(bird.data)] - bird.data$time[1:(nrow(bird.data) - 1)] #get differences between all consecutive points
  sort(diff)
 
   # Threshold values (in minutes)
  thresholds <- c(35, 50, 65)
  for (threshold in thresholds) {  #iterate over each threshold value
    indices <- which(diff > threshold) #find gaps that exceed thresholds
    for (i in indices) {
      if (i < nrow(bird.data)) {  #ensure we don't go out of bounds
        rows <- c(i, i + 1)  #get front and end of gap
        gaps <- bird.data[rows, ]  #make a dataframe showing the gap
        
        # Start making figures!
        id_title <- as.character(gaps[1, "ID"]) #change ID to character so we can use it as title
        threshold_text <- paste("     >", threshold, "min") #add threshold text
        title <- paste("ID:", id_title, threshold_text) #combine ID and threshold text
            time1_sub <- paste("Time 1:", gaps[1, "time"], " (Row:", rows[1], ")") #subtitle showing the datetime and row 
            time2_sub <- paste("Time 2:", gaps[2, "time"], " (Row:", rows[2], ")") #subtitle showing the datetime and row 
        
            # Make a map zoomed in on PI
        map1 <-  ggplot() +
          geom_sf(data = wa, color = "lemonchiffon3", fill = "ivory2", lwd = 0.75) +   #make base land map
          coord_sf(expand = F,
                   xlim = c(lon_min, lon_max),
                   ylim = c(lat_min, lat_max)) +   #zoom in on area of interest
          theme(panel.background = element_rect(fill = "white"), #background
                panel.border = element_rect(color = "black", linewidth = 1,
                                            linetype = "solid", fill = NA)) + #border around map image
          geom_point(data = gaps[,2:3], aes(x = y, y = x), col = "red", shape = 4, size = 2)
        
        # Make a map zoomed out
        map2 <- ggplot() +
          geom_sf(data = wa, color = "lemonchiffon3", fill = "ivory2", lwd = 0.75) +   #make base land map
          coord_sf(expand = F,
                   xlim = c(lon_min - 0.1, lon_max + 0.1),
                   ylim = c(lat_min - 0.1, lat_max + 0.1)) +   #zoom in on area of interest
          theme(panel.background = element_rect(fill = "white"), #background
                panel.border = element_rect(color = "black", linewidth = 1,
                                            linetype = "solid", fill = NA)) + #border around map image
          geom_point(data = gaps[,2:3], aes(x = y, y = x), col = "red", shape = 4, size = 2)
        
        # Arrange subtitles vertically
        time_sub <- arrangeGrob(
          textGrob(time1_sub, gp = gpar(fontsize = 12, fontface = "italic")),
          textGrob(time2_sub, gp = gpar(fontsize = 12, fontface = "italic")),
          ncol = 1
        )
        
        # Arrange title, subtitle, and plots
        grid.arrange(
          arrangeGrob(textGrob(title, gp = gpar(fontsize = 14)), NULL,  # Add title at the top
                      time_sub,
                      ncol = 2),
          arrangeGrob(map1, map2, ncol = 2), #maps next to each other
          nrow = 2
        )
      }
    }
  }
}
