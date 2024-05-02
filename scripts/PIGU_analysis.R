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
library(sf)

source(here("Scripts", "supportingScripts/utility.R"))

############
## PIGU DATA
# Load PIGU data
PIGU_data <- read.csv(here("data", "PIGU_data", "PIGU_data_UTM.csv"))
PIGU_data$time <- as.POSIXct(PIGU_data$time,tz="UTC") #convert times to POSIX 

# Multiple imputation to address temporal irregularity
# fit crawl model
crwOut <- crawlWrap(obsData = PIGU_data, timeStep = "15 mins",
                    theta=c(0,0), fixPar=c(NA,NA),
                    retryFits=10) #retry fits until convergence

PIGU_data_crw <- crwOut$crwPredict[,c(3,6,9:18)] #isolate interpolated data

#Create time of day covariate
PIGU_data_crw$tod <- as.numeric(format(PIGU_data_crw$time, "%H")) + as.numeric(format(PIGU_data_crw$time, "%M"))/60 #decimal hours

#Finish data
PIGU_data <- PIGU_data_crw[,c(3,6,9:18)]


## CANT FIGURE OUT HOW TO RASTERIZE
# # Create "distance to nest" covariate
# utm_proj <- "+proj=utm +zone=10 +north +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs" #define PROJ string
# PIGU_origin <- st_as_sf(PIGU_data, coords = c("x", "y"), crs = 32610) #EPSG:32610 for WGS84 10M
# PIGU_dest <- st_as_sf(PIGU_data, coords = c("nest_x", "nest_y"), crs = 32610)
# PIGU_data$dist2nest <- st_distance(PIGU_origin, PIGU_dest, by_element = TRUE) #calculate dist2nest

#################
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
install.packages("rgdal")
utm_proj <- "+proj=utm +zone=10 +north +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs" #define UTM projection string
raster_bathy_utm <- projectRaster(raster_bathy, crs = utm_proj) #project to UTM
crop <- extent(c(xmin = 475, xmax = 525, ymin = 5300, ymax = 5350)) #set up layer crop
bathy_crop <- crop(raster_bathy_utm, crop) #crop raster to specified boundaries

# Close the NetCDF file
nc_close(bathydata)
#################################
## MODEL SETUP

covlist <- list(bathy = bathy_crop)

# Prepare data for momentuHMM and calculate gradients
tracks <- prepData(data = PIGU_data,
                   type = "UTM",
                   coordNames = c("mu.x", "mu.y"),
                   covNames = "tod",
                   spatialCovs = covlist,
                   gradient = TRUE,
                   altCoordNames="mu") #if raster stack, must have z values (time, date, etc.)


nbStates <- 3 #three hidden states
dist <- list(mu="rw_mvnorm2") #random walk movement model with bivariate normal distribution
formula <- ~cosinor(tod,period=24) #cosinor function to model rhythmic patterns. Use 24hr to align with daily rhythm
stateNames <- c("encamped","exploratory", "foraging") #the states we want to identify
stateCols <- c("#339FFF","#FFF333", "#FF3300") #colors for different states


# specify model
DM <- list(mu=list(mean.x=~0+mu.x_tm1+crw(mu.x_tm1)+langevin(bathy.x),
                   mean.y=~0+mu.y_tm1+crw(mu.y_tm1)+langevin(bathy.y),
                   sd.x=~1,
                   sd.y=~1,
                   corr.xy=~1))

fixPar <- list(mu=c(NA,1,2,
                    NA,3,4,
                    NA,1,2,
                    NA,3,4,5,6,8,NA,NA))

PIGU_Fit <- fitCTHMM(tracks,Time.name="time",nbStates=nbStates,dist=dist,DM=DM,formula=formula,
                     Par0=list(mu=c(1,0,0,1,0,0,1,0,0,1,0,0,-4,-2,-4,-2,0,0)),fixPar=fixPar,
                     optMethod="TMB",control=list(silent=TRUE,trace=1),stateNames=stateNames,mvnCoords="mu") ##TWO ERRORS: Dimension mismatch and time differences?
                               






