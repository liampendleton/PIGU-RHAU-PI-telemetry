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
## RHAU DATA
# load RHAU data
RHAU_data <- read.csv(here("data", "RHAU_data", "RHAU_data.csv"))

# create time of day covariate
RHAU_data$date_time <- as.POSIXct(RHAU_data$time,tz="UTC")
RHAU_data$tod <- as.numeric(format(RHAU_data$date_time, "%H")) + as.numeric(format(RHAU_data$date_time, "%M"))/60

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

utm_proj <- "+proj=utm +zone=10 +north +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs" #define UTM projection string
raster_bathy_utm <- projectRaster(raster_bathy, crs = utm_proj) #project to UTM
crop <- extent(c(xmin = 400, xmax = 600, ymin = 5225, ymax = 5400)) #set up layer crop
bathy_crop <- crop(raster_bathy_utm, crop) #crop raster to specified boundaries

# Split RHAU_data by ID
grouped_data <- split(RHAU_data, RHAU_data$ID)

# Plot it!
for (id in names(grouped_data)) {
  plot(bathy_crop, xlim = c(400, 600), asp = 1)
  points(grouped_data[[id]]$x, grouped_data[[id]]$y, col = "red", pch = 20, cex = 0.5)
  title(paste("Individual ID:", id))
}

# Close the NetCDF file
nc_close(bathydata)