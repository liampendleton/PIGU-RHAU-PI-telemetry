library(raster)
library(ggplot2)
library(scico)
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

source("supportingScripts/utility.R")

# load zebra data from Michelot et al. (2020)
raw <- read.csv("data/zebra.csv")

# remove missing locations
raw <- raw[!is.na(raw$x),]

# create time of day covariate
raw$date_time <- as.POSIXct(raw$time,tz="UTC")
raw$tod <- as.numeric(format(raw$date_time, "%H")) + as.numeric(format(raw$date_time, "%M"))/60

# load habitat data
hb <- raster("data/vegetation2.grd")

# calculate distance to habitat type
grass <- hb
grass[hb>1] <- NA
grass <- raster::distance(grass)
names(grass) <- "grass"

bushgrass <- hb
bushgrass[hb!=2] <- NA
bushgrass <- raster::distance(bushgrass)
names(bushgrass) <- "bushgrass"

bush <- hb
bush[hb!=3] <- NA
bush <- raster::distance(bush)
names(bush) <- "bush"

wood <- hb
wood[hb<4] <- NA
wood <- raster::distance(wood)
names(wood) <- "wood"

covlist0 <- list(grass=grass,bushgrass=bushgrass,bush=bush,wood=wood)

# prepare data for momentuHMM and calculate gradients
tracks <- prepData(raw,covNames="tod",spatialCovs=covlist0,gradient=TRUE,altCoordNames="mu")

nbStates <- 2
dist <- list(mu="rw_mvnorm2")
formula <- ~cosinor(tod,period=24)
stateNames <- c("encamped","exploratory")
stateCols <- c("#56B4E9","#E69F00")

# specify model
DM <- list(mu=list(mean.x=~0+mu.x_tm1+crw(mu.x_tm1)+langevin(grass.x)+langevin(bushgrass.x)+langevin(bush.x)+langevin(wood.x),
                   mean.y=~0+mu.y_tm1+crw(mu.y_tm1)+langevin(grass.y)+langevin(bushgrass.y)+langevin(bush.y)+langevin(wood.y),
                   sd.x=~1,
                   sd.y=~1,
                   corr.xy=~1))

# use fixPar argument to fix and constrain parameters
fixPar <- list(mu=c(NA,1,2,3,4,5,NA,6,7,8,9,10,NA,1,2,3,4,5,NA,6,7,8,9,10,11,12,11,12,NA,NA))

zebraFit <- fitCTHMM(tracks,Time.name="date_time",nbStates=nbStates,dist=dist,DM=DM,formula=formula,
                     Par0=list(mu=c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,-4,-2,-4,-2,0,0)),fixPar=fixPar,
                     optMethod="TMB",control=list(silent=TRUE,trace=1),stateNames=stateNames,mvnCoords="mu")

# plot pseudo-residuals
plotPR(zebraFit)

# plot stationary distribution as a function of time of day
stpr <- plotStationary(zebraFit,plotCI=TRUE,col=stateCols,return=TRUE)

# Viterbi-decoded states
st <- viterbi(zebraFit)

# state 1 UD
XB1 <- logUD1 <- grass * zebraFit$CIbeta$mu$est[3] + bushgrass * zebraFit$CIbeta$mu$est[4] + bush * zebraFit$CIbeta$mu$est[5] + wood * zebraFit$CIbeta$mu$est[6]
values(logUD1) <- getValues(XB1) - log(sum(Brobdingnag::brob(getValues(XB1)))) #normalize

# state 2 UD
XB2 <- logUD2 <- grass * zebraFit$CIbeta$mu$est[9] + bushgrass * zebraFit$CIbeta$mu$est[10] + bush * zebraFit$CIbeta$mu$est[11] + wood * zebraFit$CIbeta$mu$est[12]
values(logUD2) <- getValues(XB2) - log(sum(Brobdingnag::brob(getValues(XB2)))) #normalize

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


