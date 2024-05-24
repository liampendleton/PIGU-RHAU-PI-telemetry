install.packages("ggmap")
install.packages("gganimate")
library(ggmap)
library(ggplot2)
library(gganimate)
library(dplyr)




install.packages("devtools")
library(devtools)
devtools::install_github("16EAGLE/moveVis")
install.packages("moveVis-0.9.9.tar.gz", repos = NULL)
library(moveVis)
library(move)

PIGU_data <- read.csv(here("data", "PIGU_data", "PIGU_data_latlon.csv"))
RHAU_data <- read.csv(here("data", "RHAU_data", "RHAU_data_latlon.csv"))
PIGU_data$time <- as.POSIXct(PIGU_data$time,tz="UTC") #convert times to POSIX
RHAU_data$time <- as.POSIXct(RHAU_data$time, tz="UTC")

# Let's reorganize data
PIGU_data <- PIGU_data %>%
  group_by(ID) %>%
  arrange(time, .by_group = TRUE)

RHAU_data <- RHAU_data %>%
  group_by(ID) %>%
  arrange(time, .by_group = TRUE)


PIGU_44067 <- PIGU_data[PIGU_data$ID == 44067,]
PIGU_44067$time <- as.POSIXct(PIGU_44067$time, format = "%Y-%m-%d %H:%M:%S")
PIGU_45657 <- PIGU_data[PIGU_data$ID == 45657,]
PIGU_45657$time <- as.POSIXct(PIGU_45657$time, format = "%Y-%m-%d %H:%M:%S")
RHAU_45663 <- RHAU_data[RHAU_data$ID == 45663,]
RHAU_45672 <- RHAU_data[RHAU_data$ID == 45672,]

# Google Key
register_google(key = "AIzaSyBHrNrmjoXOcEyBmXKiO_ensZUfICCznqU", write = TRUE)

map_bounds <- c(left = -123.074633, #PIGU
                bottom = 48.012008,
                right = -122.781436,
                top = 48.239799)

map_bounds <- c(left = -123.1, #RHAU
                bottom = 48.0,
                right = -122.5,
                top = 48.65)

sat_map <- get_map(location = map_bounds, #switch out for PIGU or RHAU
                   source = "google",
                   maptype = "satellite",
                   zoom = 10)


p <- ggmap(sat_map) +
     geom_point(data = RHAU_45672, aes(x = y, y = x), color = "cyan", size = 3) +
     geom_path(data = RHAU_45672, aes(x = y, y = x), color = "steelblue1", lwd = 1, alpha = 0.5) +
     labs(title = 'Movement of RHAU 45672', x = 'Longitude', y = 'Latitude') +
     theme_minimal()


anim <- p + transition_reveal(along=time) +
  ease_aes('linear')

animate(anim, duration = 25, fps = 5, width = 800, height = 600)
              
anim_save("RHAU_45672.gif", animation=anim, path=here("results", "RHAU_results"), fps = 5, duration = 25)








# Create a cumulative index for each observation in the correct order
PIGU_data$index <- seq_along(PIGU_data$ID)

g <- ggmap(sat_map) +
  geom_point(data = PIGU_data, aes(x = y, y = x, group = ID, color = factor(ID)), size = 3) +
  geom_path(data = PIGU_data, aes(x = y, y = x, group = ID, color = factor(ID)), lwd = 1, alpha = 0.5) +
  scale_color_viridis_d() +  # Ensures distinct colors for each ID
  labs(title = 'Tracking Movements of Individuals', x = 'Longitude', y = 'Latitude') +
  theme_minimal()

# Use transition_reveal to progressively reveal points and paths
anim <- g + 
  transition_reveal(along = as.numeric(PIGU_data$time)) +
  ease_aes('linear')

# Attempt to animate again
animate(anim, duration = 30, fps = 10)






