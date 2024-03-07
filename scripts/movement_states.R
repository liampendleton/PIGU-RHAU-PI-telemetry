## assign behavioral state to PIGU GPS setwd("/Users/stellasolasz/Desktop/PIGU_GPS_2023_COPY/complete_tags_GPS")

#install.packages("fitdistrplus")
require(fitdistrplus)
require(dplyr)
#install.packages("datatable")
library(data.table)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("marmap")
library(marmap)
#install.packages("raster")
library(raster)
#install.packages("viridis")
library(viridis)
#install.packages("scales")
library(scales)
#install.packages("moveHMM")
library(moveHMM)
#install.packages("lubridate")
library(lubridate)
library(raster)
#install.packages("metR")
library(metR)

library(ggplot2)
library(sf)
library(dplyr)
library(ggplot2)
library(ggmap)
library(sp)
library(readr)
install.packages("pak")
pak::pak('thomasp85/gganimate')
library(gganimate)
install.packages ("gapminder")
library(gapminder)
install.packages("magrittr")
library(magrittr)
install.packages("gifski")
library(gifski)
install.packages("ggmap")
library(ggmap)
### read in all of the data
## batch read in all files
# clean up GPS data ===========================================================================================================
## batch read in all files

setwd("/Users/stellasolasz/Desktop/PIGU_GPS_2023_COPY/complete_tags_GPS")

files <- list.files(path = "/Users/stellasolasz/Desktop/PIGU_GPS_2023_COPY/complete_tags_GPS",pattern = ".csv")
result <- rbindlist(sapply(files, fread,simplify = FALSE,skip = 0, header = TRUE), idcol = 'filename')
print(result) 
##colnames(result) <- c("file","loc","site", "bandcombo", "mass_deployment", "mass_retireval", "sex", "day","month","year","hour","min","sec","second_of_day","number_of_satellites","latitude","longitude","altitude","clock.offset", "accuracy", "batt", "processing_para", "processing_para")
colnames(result) <- c("ID", "loc", "site", "bandcombo", "mass_d", "mass_r", "sex", "day", "month", "year","hour", "min", "sec", "second_d", "number_of_satellites", "latitude", "longitude", "altitude", "clock_offset", "accouracy", "battery", "processing_para", "processing_para2")
view(result)
## make date column
result$date <- paste(result$day,result$month,result$year,sep="/")
result$time <- paste(result$hour,result$min,result$sec,sep=":")
result$date.time <- paste(result$date,result$time,sep=" ")
#result$date.time <- as.Date(result$date.time, format =  "%d/%m/%y %H:%M:%S")
result$date.time <- as.POSIXct(result$date.time,format="%d/%m/%y %H:%M:%S",tz=Sys.timezone())
print(result)
colnames(result)
names(result)
print(result)
## filter for accuracy 
filtered <- result %>% 
  filter(number_of_satellites >= 7,latitude > 37 & latitude < 39, longitude < -122.8, longitude > -123.2)

## create tag label
filtered$ID <- sub(".csv$", "", filtered$ID) ##get ride of how many spaces over i want to keep 
print(filtered)

summary(filtered$latitude)
summary(filtered$longitude)

## convert lat lon to utm and add to dataframe
latlon <- filtered
coordinates(latlon) <- c("longitude", "latitude") ##duplicating data frame
proj4string(latlon) <- CRS("+proj=longlat +datum=WGS84")  ## for example change CRS to WGS84
#view(latlon)

res <- spTransform(latlon, CRS("+proj=utm +zone=10 ellps=WGS84")) ##changing to a new latlong
utms <- as.data.frame(res) 
#view(res)

##plot points here--> MAYBE TRY LATER 
res_sf <- st_as_sf(res)

# Create a scatter plot of ALL tagged birds GPS points
ggplot() +
  geom_sf(data = res_sf, color = "blue", size = 0.5) +
  labs(xlab = "Longitude", ylab = "Latitude") +
  ggtitle("Scatter Plot of GPS Points")

#Enter in Diet Data============
diet <- read_csv("completediet.csv")

diet$date <- as.Date(paste(diet$day, diet$month, diet$Year, sep = "/"))
diet$time <- sprintf("%02d:%02d:00", diet$hours, diet$minutes)  # Ensure leading zeros and add seconds
diet$date_time <- as.POSIXct(paste(diet$date, diet$time), format = "%Y-%m-%d %H:%M:%S", tz = "Your_Timezone_Here")
# View the data
View(diet)

#FILTER OUT to look more closely=========================================
subset_data <- filtered[filtered$ID == "Tag_45621_LHH48", ]

# Create a scatter plot of GPS points with labels
ggplot(subset_data, aes(x = longitude, y = latitude, label = point_number)) +
  geom_point() +
  labs(xlab = "Longitude", ylab = "Latitude") +
  ggtitle("Subset Plot of GPS Points for Tag LHH48 with Point Numbers")

### Read in Meta Data ========================================================== I gotta fix my meta data

meta_deployments <- read.csv("/Users/stellasolasz/Desktop/PIGU_GPS_2023_COPY/META_DATA\ /Pigu_deployment_Summer\ 2023..csv")

all.dat <- filtered %>% full_join(meta_deployments, by = 'file') ##new data frame to be working from


## remove bad or incomplete tracks ===============================================
## This is filtering out all the tags/files from all.dat that are incomplete. Thats why data good only has one track. 
##I think All my tracks are complete. Idk how i would know which file names would be incomplete?

data.good <- filtered %>% filter(ID != "Tag_45415_GG14") %>% 
                                  mutate(date = as_datetime(date.time)) #reclassifying as date time column 
unique(data.good$ID)
dat.sub <- data.good

write.csv(dat.sub, "full_gps_all.csv", row.names = FALSE)

### loop through and save each track object for each id
## make track data 
##amt is animal movement tools
install.packages("amt")
library(amt)

## do this to start things off === I'm confused what this section is doing? 
i = "Tag_44169_GG01"
sub <- dat.sub %>% filter(ID == i)

sub <- sub[complete.cases(sub[, c("longitude", "latitude", "date.time", "ID")]), ]
print(sub)

default_value <- 0 
sub$date.time[is.na(sub$date.time)] <- default_value

# Print the summary of date.time
print(summary(sub$date.time))

trk <- make_track(sub, longitude, latitude, date.time, id = file,all_cols = TRUE)

### assign bursts for gaps in 30 seconds, need a constant step duration! Is this the correct Rates?
trk2 <- track_resample(trk, rate = seconds(30), tolerance = seconds(10))
trk2$ID <- "start"

resamp.dat <- trk2

## now the loop
trips <- unique(dat.sub$ID)

for(i in trips) {
  sub <- dat.sub %>% filter(ID == i)
  
  trk <- make_track(sub, longitude, latitude, date.time, id = ID,all_cols = TRUE)
  
  ### assign bursts for gaps in 5 minutes, need a constant step duration
  trk2 <- track_resample(trk, rate = seconds(30), tolerance = seconds(10))
  trk2$ID <- i
  
  resamp.dat <- rbind(resamp.dat,trk2)
}
#Overall, this loop processes each unique ID in the dat.sub dataset, creates resampled tracks for each ID, and aggregates them into a single dataset resamp.dat.

## get rid of the first step
resamp.data <- resamp.dat %>% filter(ID != "start")


## identify individual tracks ==========================================================================
library(ggplot2)
library(sf)
library(maps)
library(mapdata)
library(dplyr)
library(stringr)
library(lubridate)
library(argosfilter)
install.packages("argosfilter")
install.packages("stringr")
install.packages("lubridate")
# install my package 
install.packages("devtools") # for installing packages from github
devtools::install_github("abfleishman/trakR",upgrade = "ask") # install my package
install.packages("trakR") 
library(trakR)

colnames(resamp.data)[1:4] <- c("Longitude", "Latitude", "DateTime", "CaptureID")
print(resamp.data)

resamp.data$Dist2Col<-trakR::Dist2Colony(tracks = resamp.data, 
                                    dataLat = "Latitude",
                                    dataLon = "Longitude",
                                    ColonyLat = 37.698355,
                                    ColonyLong = -123.003786)

tracks_w_trips<-trakR::MakeTrip(tracks = resamp.data,
                                ID = "CaptureID",
                                DistCutOff = 0.50,
                                Dist2Colony = "Dist2Col",
                                NumLocCut = 15)

resamp.data <- resamp.data %>%
 mutate(Dist2Col = trakR::Dist2Colony(tracks = resamp.data, 
                                      dataLat = "Latitude",
                                      dataLon = "Longitude",
                                      ColonyLat = 37.698355,
                                      ColonyLong = -123.003786))

head(tracks_w_trips)

png(file="tracks.png", units="in",width=14,height=16,res=500)
ggplot(tracks_w_trips,
       aes(Longitude,
           Latitude,
           col=factor(TripNum)))+
  geom_path(linewidth=.7)+
  geom_point(data=tracks_w_trips[tracks_w_trips$ColonyMovement%in%c("Out","In"),])+
  theme_classic(base_size = 16)+
  labs(color="TripNum")+
  facet_wrap(~CaptureID,scales = "free")
dev.off()

#i want to look at each individual bird and zoom into the craziness so alter x and y limits 
subset_gg01 <- tracks_w_trips[tracks_w_trips$CaptureID == "Tag_44169_GG01", ]#subset out data to look at individual birds

subset_gg01 <- subset_gg01[subset_gg01$TripNum %in% c(0),]

# Create a scatter plot of GPS points with labels
ggplot(subset_gg01,
       aes(Longitude,
           Latitude,
           col=factor(TripNum)))+
  geom_path(linewidth=.7)+
  geom_point(data=subset_gg01[subset_gg01$ColonyMovement %in% c("Out", "In"),])+
  theme_classic(base_size = 16)+
  labs(color="TripNum")

##animate
farallon_shape <- st_read("/Users/stellasolasz/Desktop/farallon_dtm_alt_0m_poly")

# Check coordinate system of shapefile
st_crs(farallon_shape)

# Check coordinate system of GPS points
st_crs(subset_gg01)
#check class 
class(subset_gg01)
#animate 
subset_gg01$DateTime <- as.POSIXct(subset_gg01$DateTime)

mymap.paths <- ggplot() + 
  geom_point(data = subset_gg01, aes(x = Longitude, y = Latitude, colour = DateTime)) +
  geom_segment(data = create_segments(subset_gg01), 
               aes(x = Longitude, y = Latitude, xend = lead(Longitude), yend = lead(Latitude), colour = DateTime), 
               size = 1) +
  labs(x = "Longitude", y = "Latitude") +
  scale_colour_continuous(name = "Time") +
  theme(legend.position = "bottom")

# Display the map with paths
mymap.paths

# Define the order of frames
frames <- unique(subset_gg01$DateTime)

# Update plot to animate
path.animate.plot <- mymap.paths +
  transition_states(DateTime,
                    transition_length = 2, # Adjust transition length as needed
                    state_length = 2) +
  labs(title = 'Date: {closest_state}')  # Add a label on top to show the closest state date

# Adjust frames per second (fps) and number of frames (nframes) as needed
animate(path.animate.plot,
        fps = 1, # frames per second
        nframes = 500, # total number of frames
        end_pause = 10)

animate(path.animate.plot, renderer = gifski_renderer())

#looking at the movement of the birds at tip titles "0"
ggplot(subset_gg01, aes(x = DateTime, y = Latitude, group = 1)) +
  geom_path() +  # Connect the points with lines
  geom_point() +  # Add points for each location
  labs(x = "Datetime", y = "Latitude", title = "Movement during Trip 0")

#i want to look at each individual bird and zoom into the craziness so alter x and y limits 
subset_LHH61 <- tracks_w_trips[tracks_w_trips$CaptureID == "Tag_44258_LHH61",] #subset out data to look at individual birds

subset_LHH61 <- subset_LHH61[subset_LHH61$TripNum %in% c(0), ]

# Plot
ggplot(subset_LHH61,
       aes(Longitude,
           Latitude,
           col=factor(TripNum)))+
  geom_path(linewidth=.7)+
  geom_point(data=subset_LHH61[subset_LHH61$ColonyMovement %in% c("Out", "In"),])+
  theme_classic(base_size = 16)+
  labs(color="TripNum")


#========194
subset_LHH194 <- tracks_w_trips[tracks_w_trips$CaptureID == "Tag_45411_LHH194",]#subset out data to look at individual birds

subset_LHH194 <- subset_LHH194[subset_LHH194$TripNum %in% c(0), ]

# Plot
ggplot(subset_LHH194,
       aes(Longitude,
           Latitude,
           col=factor(TripNum)))+
  geom_path(linewidth=.7)+
  geom_point(data=subset_LHH194[subset_LHH194$ColonyMovement %in% c("Out", "In"),])+
  theme_classic(base_size = 16)+
  labs(color="TripNum")

#===== talk to mike / scott about this issue of tracks not being labeled properly and the trip zero being isolated around the island?

## add in and out BACK TO OUR NORMAL PROGRAM ===================== 
library(dplyr)

tracks_w_trips$INOut <- InOutPoints(
  tracks = tracks_w_trips,
  CaptureID = "CaptureID",
  DateTime = "DateTime",
  TripID = "TripNum",
  dist2colony = "Dist2Col",
  lag = 2,
  nPointsToSmooth = 10,
  minDist2Col = 10,
  Lon = "Longitude",
  Lat = "Latitude",
  Plot = TRUE,
  pdfName = "inout_plots.pdf"
)

##Plot Data testing in and out tracks 
ggplot(tracks_w_trips[tracks_w_trips$TripNum%in%c(1,2,3,4,5,6,7,8,9,10),],
       aes(Longitude,
           Latitude,
           col=factor(TripNum)))+
  geom_path(aes(lty=INOut),linewidth=.75)+
  scale_linetype_manual(values = c(1,2,3))+
  facet_wrap(~TripNum, scales = "free")+
  theme_classic(base_size = 16)+
  labs(color="TripNum")

# Step 2: Aggregate the maximum distances by CaptureID
max_distances <- resamp.data %>%
  group_by(CaptureID) %>%
  summarize(max_distance = max(Dist2Col))

# View the maximum distances
print(max_distances)
png(file="maxdistance.png", units="in",width=14,height=16,res=500)

ggplot(max_distances, aes(x = CaptureID, y = max_distance)) +
  geom_bar(stat = "identity", fill = "darkgrey", color = "black") +
  labs(x = "CaptureID",
       y = "Maximum Distance from Colony (Km)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())  # Remove minor grid lines
ggsave("plot.png", plot = last_plot(), width = 7, height = 5)

#===============================================================================
  
### create new dataframe excluding trips with less than 10????locations===
install.packages("dplyr")
library(dplyr)
##? MIKE/chatgpt 
new <- tracks_w_trips %>% 
  group_by(CaptureID, TripNum) %>% 
  mutate(Var=max(Dist2Col)) %>% 
  filter(Var>1) ##how do i deside how many locations to exclude?

## Calculate the number of locations in each trip
new <- new %>% 
  group_by(CaptureID, TripNum) %>% 
  mutate(NumLocations = n())

png(file="tracksNEW.png", units="in",width=14,height=16,res=500)
ggplot(new,
       aes(Longitude,
           Latitude,
           col=factor(TripNum)))+
  geom_path(linewidth=.7)+
  geom_point(data=tracks_w_trips[tracks_w_trips$ColonyMovement%in%c("Out","In"),])+
  theme_classic(base_size = 16)+
  labs(color="TripNum")+
  facet_wrap(~CaptureID,scales = "free")
dev.off()

  ggplot(new, aes(x = NumLocations)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Trip Lengths",
       x = "Number of Locations in Trip",
       y = "Frequency") +
  theme_minimal()  
  
  ggplot(new, aes(x = as.factor(CaptureID), y = Dist2Col, fill = as.factor(CaptureID))) +
    geom_boxplot() +
    labs(title = "Boxplot of Dist2Col by Capture ID",
         x = "Capture ID",
         y = "Dist2Col") +
    theme_minimal()
summary(new$NumLocations)

### how long did each trip take ============================
timing <- new %>%
  group_by(CaptureID, TripNum, INOut) %>%
  summarize(CaptureID = unique(CaptureID),
    min = min(DateTime, na.rm = TRUE),
    max = max(DateTime, na.rm = TRUE),
    dur = max(DateTime, na.rm = TRUE) - min(DateTime, na.rm = TRUE)
  )

timing$minutes <- as.numeric(timing$dur / 60)

mean_duration <- timing %>%
  group_by(CaptureID) %>%
  summarize(mean_duration_minutes = mean(minutes))

summary(mean_duration$mean_duration_minutes)

bar_chart <- ggplot(mean_duration, aes(x = CaptureID, y = mean_duration_minutes)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  labs(title = "Mean Duration for Each CaptureID",
       x = "CaptureID",
       y = "Mean Duration (minutes)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(bar_chart)

#save timing data as CSV.
write.csv(timing, "timing.csv", row.names = FALSE)

ggplot(timing)+
  geom_boxplot(aes(INOut,minutes,fill = CaptureID))

##HIST of timing
install.packages("lme4")
require(lme4)
mod <- glm(minutes ~ INOut * CaptureID , data = timing)
hist(timing$minutes)
summary(mod)

summary_file <- "model_summary.txt"
cat(capture.output(summary(mod)), file = summary_file)

### HMM TIME!!!! LETS GGOOOOOO. fit a momentumHMM model ======================================================

##================MIKE CODE===================================================
data.move <- prepData(new, type="UTM", coordNames=c("Longitude","Latitude"))
glimpse(new)

hist(data.move$step, xlab = "step length", main = "")

# Indices of steps of length zero
whichzero <- which(data.move$step == 0)
# Proportion of steps of length zero in the data set
length(whichzero)/nrow(data.move)

hist(data.move$angle, xlab = "angle length", main = "")

summary(data.move$step)
summary(data.move$angle)

ggplot(data.move, aes(x = DateTime, y = step)) +
  geom_line() +
  labs(title = "Step Length Over Time",
       x = "Date and Time",
       y = "Step Length (meters)") +
  theme_minimal()

### starting values===I THINK IM STRUGGLING WITH FIGURING OUT MY PERAMETERS. THEY ARE SO CLOSE TO THE ISLAND AND ALSO SO FAR
#Forage, travel, rest
mu0 <- c(.01,.01,0.1) # step mean (two parameters: one for each state)
sigma0 <- c(.015,.015,0.2) # step SD
zeromass0 <- c(0.01,0.01,0.01) # step zero-mass
stepPar0 <- c(mu0,sigma0,zeromass0)

angleMean0 <- c(pi,0,0) # angle mean
kappa0 <- c(1,9,10) # angle concentration
anglePar0 <- c(angleMean0,kappa0)

## to reproduce results
set.seed(123)

length(stepPar0)  # Should be 6 for 3 states
length(anglePar0) # Should be 2 for 3 states

hmm_2states <- fitHMM(data=data.move, nbStates=3, stepPar0=stepPar0, anglePar0=anglePar0, formula=~1)
summary(hmm_2states)

plot(hmm_2states)

data.move$states <- viterbi(hmm_2states)

## have a look at things

foraging <- data.move %>% filter(states == 1)

## let's look at the tracks
png(file="trackshmm.png", units="in",width=14,height=16,res=500)
ggplot(data.move,
       aes(x,
           y,
           col=factor(TripNum)))+
  geom_path(linewidth=.7)+
  geom_point(data=foraging,aes(x,y),color="black",size=0.2)+
  theme_classic(base_size = 16)+
  labs(color="TripNum")+
  facet_wrap(~CaptureID, scales = "free")
dev.off()

ggplot(data.move,
       aes(x,
           y,
           col=factor(hour),
           group=TripNum))+
  geom_path(size=.7)+
  #geom_point(data=foraging,aes(long,lat),color="black",size=0.2)+
  theme_classic(base_size = 16)+
  labs(color="TripNum")+
  facet_wrap(~CaptureID)


### data for plotting ------------------------------------------------------------------------

# location of islands

islands <- data.frame(long = c(-123.003786,-123.031378,-123.103514),
                      lat = c(37.698355,37.727854,37.768175),
                      name =c("south","middle","north"))


## land layer
setwd("/Users/mike/Dropbox/Analyses/MMC Crab Fleet")

states    <- c('California', 'Nevada', 'Oregon', 'Washington','Alaska')

us <- getData("GADM",country="USA",level=1)

us.states <- us[us$NAME_1 %in% states,]

setwd("/Users/mike/Dropbox/Analyses/CAAU issa")

# bathymetry data

#Load shapefile
shp <- raster("/Users/mike/Dropbox/Analyses/CAAU issa/Data/GEBCO_12_Oct_2023_259d0c65e220/gebco_2023_n38.3835_s37.1077_w-123.9955_e-122.3222.asc")
plot(shp)
class(shp)

conts <- as.data.frame(shp,xy=T)

## colors for bathymetry plotting
greys<-c("aliceblue","grey45","grey10")

ggplot()+
  # draw topography and bathymetry
  # draw topography and bathymetry
  geom_raster(data=conts,aes(x,y,fill=gebco_2023_n38.3835_s37.1077_w.123.9955_e.122.3222))+
  #geom_contour(data=conts, aes(x=lon,y=lat,z=value),color="grey30",size=.09, breaks=seq(0,-5000,-100))+ # sea contours
  # add state and province borders
  # plot boat locationa
  geom_path(data=data.move,aes(long,lat,group=ID,color=factor(year.y)),size=0.3,alpha=0.3)+
  geom_point(data=foraging,aes(long,lat,color=factor(year.y)),size=0.2)+
  geom_polygon(data=us.states,aes(x=long,y=lat,group=group),size=0.3)+
  
  geom_contour(data=conts, aes(x,y,z=gebco_2023_n38.3835_s37.1077_w.123.9955_e.122.3222),color="grey15", size=0.1,alpha=0.5, limits = c(0,-200),breaks=seq(-20,-200,-20))+
  geom_point(data=islands,aes(long,lat))+
  #scale_fill_viridis(option = "A")+
  scale_fill_gradientn(colors = rev(greys),limits=c(-200,0),na.value="aliceblue")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  #coord_sf(expand = FALSE) 
  #coord_sf(xlim=c(-123.1,-122.9),ylim=c(37.65,37.75),expand = FALSE) ## study region
  coord_sf(xlim=c(-123.6,-122.5),ylim=c(37.4,38.2),expand = FALSE) ## study region


## what hours did foraging happen

data.move$Time <- as.POSIXct(strptime(data.move$time, "%H:%M:%S"),tz = "UTC")
data.move$Time.pct <- with_tz(data.move$Time, "America/Los_Angeles")

  ggplot()+
  geom_point(data=data.move,aes(Time.pct,ID,color=factor(states)))+
  scale_x_datetime(breaks = breaks_width("1 hour"),labels=date_format("%H:%M",tz="America/Los_Angeles"))

  ggplot(data.move, aes(x=factor(hour), fill = factor(states))) +
    geom_bar(stat="count",position = "fill")+
    facet_wrap(~year.y)
  
  ### how much time spent foraging
  
 ggplot()+
   geom_point(data=data.move,aes(long,lat,color=class),size=0.3)+
   facet_wrap(~factor(hour))
 
 
 ###############################################################################################
 ### summarize trips ====================================================================
 ################################################################################################
 
 t <- data.move %>% 
   group_by(ID,TripNum) %>%
   summarize(dist=sum(step,na.rm=T),class = unique(class))
  
 
 ggplot(t,aes(class,dist))+
          geom_boxplot()

        