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

source(here("Scripts", "supportingScripts/utility.R"))

# load PIGU data
PIGU_data <- read.csv(here("data", "PIGU_data", "PIGU_data.csv"))

# create time of day covariate
PIGU_data$date_time <- as.POSIXct(PIGU_data$time,tz="UTC")
PIGU_data$tod <- as.numeric(format(PIGU_data$date_time, "%H")) + as.numeric(format(PIGU_data$date_time, "%M"))/60
####THIS IS DOWN TO LINE 24 IN "zebraAnalysis.R"####

