library(here)
library(raster)
library(ggplot2)
install.packages("scico")
library(scico)
install.packages("Brobdingnag")
library(Brobdingnag)
if(!requireNamespace("RandomFieldsUtils",quietly=TRUE)){
  install.packages("RandomFieldsUtils_1.2.5.tar.gz", repos = NULL, type = "source") # most recent archived version; required by RandomFields
}
if(!requireNamespace("RandomFields",quietly=TRUE)){
  install.packages("RandomFields_3.3.14.tar.gz", repos = NULL, type = "source") # most recent archived version; required by Rhabit
}
remotes::install_github("papayoun/Rhabit",dependencies = TRUE)
library(Rhabit)
if(!requireNamespace("momentuHMM",quietly=TRUE) || packageVersion("momentuHMM")<2){
  remotes::install_github("bmcclintock/momentuHMM@develop",dependencies = TRUE) # requires momentuHMM version >= 2.0.0
}
library(momentuHMM)
library(ncdf4)
library(gridExtra)

source(here("Scripts", "supportingScripts/utility.R"))

### DATA FORMATTING
## PIGU DATA
# load PIGU data
PIGU_data <- read.csv(here("data", "PIGU_data", "PIGU_data.csv"))

# create time of day covariate
PIGU_data$date_time <- as.POSIXct(PIGU_data$time,tz="UTC") #convert times to POSIX 
PIGU_data$tod <- as.numeric(format(PIGU_data$date_time, "%H")) + as.numeric(format(PIGU_data$date_time, "%M"))/60

## COVARIATE DATA
# load in bathymetry data
bathydata <- nc_open(here("data", "gebco_2023_n48.8518_s47.1651_w-124.837_e-122.1564.nc")) #obtained from NOAA's ETOPO 2022 Grid Extract on 4/8/2024

# Get the variables from the NetCDF file
lon <- ncvar_get(bathydata, "lon")
lat <- ncvar_get(bathydata, "lat")
bathy <- ncvar_get(bathydata, "elevation")

# Create a raster object
raster_bathy <- raster(t(bathy), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat), crs = "+proj=longlat +datum=WGS84")
raster_bathy <- flip(raster_bathy, "y") #flip to proper orientation

# Project to UTM coordinates
library(rgdal)
utm_proj <- "+proj=utm +zone=10 +north +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs" #define UTM projection string
raster_bathy_utm <- projectRaster(raster_bathy, crs = utm_proj) #project to UTM
crop <- extent(c(xmin = 475, xmax = 525, ymin = 5300, ymax = 5350)) #set up layer crop
bathy_crop <- crop(raster_bathy_utm, crop) #crop raster to specified boundaries

# Plotting tracks on top of bathymetry data; just for visualization
# Split PIGU_data by ID
grouped_data <- split(PIGU_data, PIGU_data$ID)

# Plot it!
for (id in names(grouped_data)) {
  plot(bathy_crop, asp = 1)
  points(grouped_data[[id]]$x, grouped_data[[id]]$y, col = "red", pch = 20, cex = 0.5)
  title(paste("Individual ID:", id))
}

#################################
# create a test track
track_44067 <- tracks[tracks$ID == 44067,]



# Plot points against time to identify gaps
plot(x = track_44067$date_time, 
     y = rep(1,nrow(track_44067)))
  








# Close the NetCDF file
nc_close(bathydata)

covlist <- list(bathy = bathy_crop)

# Prepare data for momentuHMM and calculate gradients
tracks <- prepData(data = PIGU_data,
                   type = "UTM",
                   coordNames = c("x", "y"),
                   covNames = "tod",
                   spatialCovs = covlist,
                   gradient = TRUE,
                   altCoordNames="mu") #if raster stack, must have z values (time, date, etc.)



nbStates <- 2
dist <- list(mu="rw_mvnorm2")
formula <- ~cosinor(tod,period=24)
stateNames <- c("encamped","exploratory")
stateCols <- c("#56B4E9","#E69F00")





# create a test track
# track_44067 <- tracks[tracks$ID == 44067,]

# # Multiple imputation to address temporal irregularity
# # fit crawl model
# crwOut <- crawlWrap(obsData = track_44067, 
#                     timeStep = seq(head(track_44067, 1), (tail(track_44067, 1)), "15 mins"), #Does this need to be looped through unique individuals in "tracks"?
#                     theta=c(6.855, -0.007), fixPar=c(NA,NA)) #IDK WHAT TO DO HERE
# IGNORE FOR NOW

dist <- list(mu="rw_mvnorm2")










# specify model
DM <- list(mu=list(mean.x=~0+mu.x_tm1+crw(mu.x_tm1)+langevin(bathy.x),
                   mean.y=~0+mu.y_tm1+crw(mu.y_tm1)+langevin(bathy.y),
                   sd.x=~1,
                   sd.y=~1,
                   corr.xy=~1))








