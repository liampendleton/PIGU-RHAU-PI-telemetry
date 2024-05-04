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

#########################################################
## PIGU DATA
# Load PIGU data
PIGU_data <- read.csv(here("data", "PIGU_data", "PIGU_data.csv"))
PIGU_data$time <- as.POSIXct(PIGU_data$time,tz="UTC") #convert times to POSIX 

# Multiple imputation to address temporal irregularity
# fit crawl model
crwOut <- crawlWrap(obsData = PIGU_data, timeStep = "15 mins",
                    theta=list(c(0.658,2.715),
                               c(0.644,2.890),
                               c(1.408,16.170),
                               c(-0.171,13.895),
                               c(0.249,2.991),
                               c(0.267,2.647)), #theta values obtained from iterating through different values and monitoring output; see below
                    fixPar=c(NA,NA))

crwOutFits <- crwOut$crwFits #observe parameter estimates and AIC/LL
print(crwOutFits) #iterate through theta values above (crawlWrap) using sigma and beta

PIGU_data_crw <- crwOut$crwPredict[,c(3,6,9:17)] #isolate interpolated data

#Create time of day covariate
PIGU_data_crw$tod <- as.numeric(format(PIGU_data_crw$time, "%H")) + as.numeric(format(PIGU_data_crw$time, "%M"))/60 #decimal hours

#Finish data
PIGU_data <- PIGU_data_crw


## CANT FIGURE OUT HOW TO RASTERIZE
# # Create "distance to nest" covariate
# utm_proj <- "+proj=utm +zone=10 +north +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs" #define PROJ string
# PIGU_origin <- st_as_sf(PIGU_data, coords = c("x", "y"), crs = 32610) #EPSG:32610 for WGS84 10M
# PIGU_dest <- st_as_sf(PIGU_data, coords = c("nest_x", "nest_y"), crs = 32610)
# PIGU_data$dist2nest <- st_distance(PIGU_origin, PIGU_dest, by_element = TRUE) #calculate dist2nest

#########################################################
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
#########################################################
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

# Background setup
nbStates <- 3 #three hidden states
dist <- list(mu="rw_mvnorm2") #random walk movement model with bivariate normal distribution
formula <- ~cosinor(tod,period=24) #cosinor function to model rhythmic patterns. Use 24hr to align with daily rhythm
stateNames <- c("encamped","exploratory", "foraging") #the states we want to identify
stateCols <- c("#339FFF","#FFF333", "#FF3300") #colors for different states

# Specify model
DM <- list(mu=list(mean.x=~0+mu.x_tm1+crw(mu.x_tm1)+langevin(bathy.x),
                   mean.y=~0+mu.y_tm1+crw(mu.y_tm1)+langevin(bathy.y),
                   sd.x=~1,
                   sd.y=~1,
                   corr.xy=~1))

# Fix parameters; transition probabilities?
fixPar <- list(mu=c(NA,1,2, #x state 1
                    NA,3,4, #x state 2
                    NA,5,6, #x state 3
                    NA,1,2, #y state 1
                    NA,3,4, #y state 2
                    NA,5,6, #y state 3
                    5,6,7, #SD for x 1:3
                    5,6,7, #SD for y 1:3
                    NA,NA,NA)) #corr xy 1:3?

# Fit the model
PIGU_Fit <- fitCTHMM(tracks,Time.name="time",nbStates=nbStates,dist=dist,DM=DM,formula=formula,
                     Par0=list(mu=c(1,0,0,
                                    1,0,0,
                                    1,0,0,
                                    1,0,0,
                                    1,0,0,
                                    1,0,0,
                                    -4,-2,-2,
                                    -4,-2,-2,
                                    0,0,0)),
                     fixPar=fixPar,
                     optMethod="TMB",control=list(silent=TRUE,trace=1),stateNames=stateNames,mvnCoords="mu")

###############################################
## Adapting McClintock code down here
# plot pseudo-residuals
plotPR(PIGU_Fit)

# plot stationary distribution as a function of time of day
stpr <- plotStationary(PIGU_Fit,plotCI=TRUE,col=stateCols,return=TRUE)

# Viterbi-decoded states
st <- viterbi(PIGU_Fit)

# state 1 UD (encamped)
XB1 <- logUD1 <- bathy_crop * PIGU_Fit$CIbeta$mu$est[3]
values(logUD1) <- getValues(XB1) - log(sum(Brobdingnag::brob(getValues(XB1)))) #normalize

# state 2 UD (exploratory)
XB2 <- logUD2 <- bathy_crop * PIGU_Fit$CIbeta$mu$est[9]
values(logUD2) <- getValues(XB2) - log(sum(Brobdingnag::brob(getValues(XB2)))) #normalize

# state 3 UD (foraging)
XB3 <- logUD3 <- bathy_crop * PIGU_Fit$CIbeta$mu$est[9]
values(logUD3) <- getValues(XB2) - log(sum(Brobdingnag::brob(getValues(XB2)))) #normalize

par(mfrow=c(1,2))
raster::plot(logUD1,col=scico(palette="roma",256,direction=-1),main=stateNames[1],xlab="easting (km)",ylab="northing (km)")
points(zebraFit$data[,c("mu.x_tm1","mu.y_tm1")][which(st==1),],col=stateCols[1],pch=20)
raster::plot(logUD2,col=scico(palette="roma",256,direction=-1),main=stateNames[2],xlab="easting (km)")
points(zebraFit$data[,c("mu.x_tm1","mu.y_tm1")][which(st==2),],col=stateCols[2],pch=20)

# mixture of both UDs based on time spent in each state
tis <- timeInStates(zebraFit) # activity budgets
logUD <- log(tis$encamped*exp(logUD1) + tis$exploratory*exp(logUD2))
plotSpatialCov(zebraFit,logUD,colors=scico(palette="roma",256,direction=-1),col=stateCols)

# plot the habitat selection coefficients
beta <- rbind(as.data.frame(lapply(zebraFit$CIbeta$mu,function(x) x[,3:6])),as.data.frame(lapply(zebraFit$CIbeta$mu,function(x) x[,9:12])))
beta$cov <- c("grassland","bushed grassland","bushland","woodland")
beta$state <- rep(stateNames,each=4)

ggplot(beta, aes(x=factor(cov,level=beta$cov[1:4]), y=est, group=state, colour=state)) + scale_color_manual(values = stateCols) +
  geom_point(position=position_dodge(width=0.25)) + geom_errorbar(aes(ymin=lower, ymax=upper),position=position_dodge(width=0.25), width=.1) +
  xlab("distance to habitat type") + ylab(expression(beta)) +
  theme(legend.text = element_text(size = rel(1.15)),legend.title=element_text(size=rel(1.35)),axis.title=element_text(size = rel(1.25)),axis.text=element_text(size = rel(1.15)))

# steps and turns by state
par(mfrow=c(2,2))
hist(zebraFit$data$step[which(st==1)],breaks=seq(0,4.5,length=40),main="encamped step lengths")
hist(zebraFit$data$angle[which(st==1)],breaks=seq(-pi,pi,length=40),main="encamped turn angles")
hist(zebraFit$data$step[which(st==2)],breaks=seq(0,4.5,length=40),main="exploratory step lengths")
hist(zebraFit$data$angle[which(st==2)],breaks=seq(-pi,pi,length=40),main="exploratory turn angles")

par(mfrow=c(4,2))
for(j in c("grass","bushgrass","bush","wood")){
  hist(zebraFit$data[[j]][which(st==1)],breaks=seq(min(zebraFit$data[[j]]),max(zebraFit$data[[j]]),length=40),main=paste0("encamped: ",j))
  hist(zebraFit$data[[j]][which(st==2)],breaks=seq(min(zebraFit$data[[j]]),max(zebraFit$data[[j]]),length=40),main=paste0("exploratory: ",j))
}

# MALA simulation
mala <- malaSim(zebraFit,list(parIndex=list(c(3:6),c(9:12)),covNames=list(c("grass","bushgrass","bush","wood"),c("grass","bushgrass","bush","wood"))),niter=50,ssl=FALSE)
apply(mala,2,mean)







