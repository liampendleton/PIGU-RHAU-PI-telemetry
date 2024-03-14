library(tidyverse)
library(rgdal)
library(terra)
library(sf)
library(sp)
library(marmap)
library(nngeo)
library(raster)
library(viridis)
library(scales)
library(ggplot2)
library(patchwork)
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

source("supportingScripts/sslCovariates.R")
#covlist <- getCovariates() # takes a while and a bit of memory, so covlist used in analysis can be loaded with "sslCovs.RData"
load("data/sslCovs.RData") # load the covariates (queried from ETOPO1 database hosted by NOAA on 6 June 2022)
load("data/tracks.RData") # load the tracks
load("data/Par0.RData") # initial parameter values for optimization

source("supportingScripts/utility.R") # utility functions

# standardize covariate rasters based on observations
covlist0 <- list()
cellInd <- cellFromXY(covlist[[1]],tracks[,c("mu.x","mu.y")])
for(i in names(covlist)){
  covlist0[[i]] <- (covlist[[i]]-mean(covlist[[i]][cellInd]))/sd(covlist[[i]][cellInd])
}

# clip raster to smaller area for plotting and extracting
covlistCrop <- lapply(covlist0,function(x) raster::crop(x,extent(-2247.282,-1212.282,198.0775,918.0775)))

# prepare data for momentuHMM and calculate gradients
tracks <- prepData(tracks,coordNames=c("mu.x","mu.y"),altCoordNames = "mu",spatialCovs=covlist0,gradient=TRUE)
tracks$intercept <- 1

# define data stream distribution
dist <- list(mu="rw_mvnorm2",dry="bern")

source("supportingScripts/DM.R") # model definitions 

fitLangevin <- UD <- list()

# original single-state model of Michelot et al. (2019)
fitLangevinOrig <- fitCTHMM(tracks,Time.name="date_time",Time.unit="hours",dist=dist,nbStates=1,
                            DM=list(mu = list(mean.x=~mu.x_tm1+langevin(depth.x)+langevin(slope.x)+langevin(d2site.x),
                                              mean.y=~mu.y_tm1+langevin(depth.y)+langevin(slope.y)+langevin(d2site.y),
                                              sd.x=~1,
                                              sd.y=~1,
                                              corr.xy=~1),
                                    dry = list(prob=~1)),
                            Par0=list(mu=c(1, -0.15059375, -0.02236392, -2.05770669,  1, -0.15059375, -0.02236392, -2.05770669, 1.00566450, 1.00566450, 0),
                                      dry=-0.4844352),
                            fixPar=list(mu=c(NA,1,2,3,NA,1,2,3,4,4,NA)),
                            workBounds=list(dry=matrix(c(-Inf,0),1,2)),
                            mvnCoords = "mu",
                            optMethod="TMB",control=list(silent=TRUE,trace=1),
                            stateNames="outbound",modelName="S=1 sigma(.) beta(depth+slope+d2site)")
UDO <- plotUD(nbUD=1,parIndex=list(2:4),covNames=list(c("depth","slope","d2site")),model=fitLangevinOrig,UDnames="outbound",UDstates=list(1),return=TRUE)

# fit models with additional states and covariates defined in `DM.R`
for(i in 1:length(DM)){
  message("\n",modelName[[i]])
  fitLangevin[[i]]<- fitCTHMM(tracks,Time.name="date_time",Time.unit="hours",dist=dist,nbStates=nbStates[[i]],DM=DM[[i]],
                               Par0=Par[[i]]$Par,
                               beta0=Par[[i]]$beta,
                               delta0=Par[[i]]$delta,
                               formula=formula[[i]],
                               fixPar=fixPar[[i]],
                               prior=prior[[i]],
                               workBounds=workBounds[[i]],
                               mvnCoords = "mu",
                               optMethod="TMB",stateNames=stateNames[[i]],modelName=modelName[[i]],kappa=kappa[[i]],
                               control=list(silent=TRUE,trace=1))
  UD[[i]] <- plotUD(inputUD[[i]]$nbUD,inputUD[[i]]$parIndex,inputUD[[i]]$covNames,model=fitLangevin[[i]],inputUD[[i]]$UDnames,inputUD[[i]]$UDstates,inputUD[[i]]$sign,return=TRUE)
}

# AIC excluding models with better fit but intermediary "foraging" state
penAIC(append(list(fitLangevinOrig),fitLangevin[c(1:5,7,9,12,14)]))

pdf("4stateUD.pdf",width=14,height=10)
statePlot(fitLangevin[[9]],UD[[9]])
dev.off()
pdf("5stateUD.pdf",width=14,height=14)
statePlot(fitLangevin[[14]],UD[[14]])
dev.off()

# pseudo-residuals
pdf("pseudoRes.pdf",width=21,height=7)
plotRes(fitLangevin[[2]])
plotRes(fitLangevin[[9]])
plotRes(fitLangevin[[14]])
dev.off()

# predicted and residuals
pdf("predRes.pdf",width=14,height=7)
plotPredRes(fitLangevin[[2]],inputUD[[2]])
plotPredRes(fitLangevin[[9]],inputUD[[9]])
plotPredRes(fitLangevin[[14]],inputUD[[14]])
dev.off()

# TPM plot (at median \Delta_t)
pdf("tpm.pdf",width=16,height=16)
tpmPlot(fitLangevin[[9]])
tpmPlot(fitLangevin[[14]])
dev.off()

# stationary distribution plot 
pdf("stationary.pdf",width=13,height=7)
statPlot(fitLangevin[[9]])
statPlot(fitLangevin[[14]])
dev.off()

# sojourn time
pdf("sojourn.pdf",width=14,height=7)
sojournTimePlot(fitLangevin[[9]])
sojournTimePlot(fitLangevin[[14]])
dev.off()

# compare simulated tracks from models 2, 9, and 14
sim <- expandUD <- list()
initPos <- mapply(function(x) c(tracks[which(tracks$ID==x)[1],c("mu.x")],tracks[which(tracks$ID==x)[1],c("mu.y")]),unique(tracks$ID),SIMPLIFY = FALSE)
for(i in c(2,9,14)){
  message("\n",modelName[[i]])
  set.seed(719647,kind="Mersenne-Twister",normal.kind = "Inversion")
  sim[[i]] <- simCTHMM(model=fitLangevin[[i]],spatialCovs=covlist0,initialPosition=initPos,states=TRUE,retrySims=5)
  expandUD[[i]] <- plotUD(inputUD[[i]]$nbUD,inputUD[[i]]$parIndex,inputUD[[i]]$covNames,model=fitLangevin[[i]],inputUD[[i]]$UDnames,inputUD[[i]]$UDstates,inputUD[[i]]$sign,cropRast=FALSE,return=TRUE)
}
slim <- apply(rbind(tracks[,c("mu.x","mu.y")],dplyr::bind_rows(sim)[,c("mu.x","mu.y")]),2,range)+c(-100,100,-50,50)

land <- (covlist0$depth*sd(covlist$depth[cellInd])+mean(covlist$depth[cellInd])-attr(covlist$depth,"absmin"))>=0

cextheme <- theme(legend.text = element_text(size = rel(1.15)),legend.title=element_text(size=rel(1.35)),axis.title=element_text(size = rel(1.25)),axis.text=element_text(size = rel(1.15)))
p1<-plotSpatialCov(tracks,crop(land,slim),colors=gray.colors(10, start = 1, end = 0.4, gamma = 2.2, alpha = NULL),alpha=0.5,col=stateCols[[4]],return=TRUE)+ggplot2::labs(fill = NULL,x=NULL,y="northing (km)")+guides(fill = "none")+theme(panel.border = element_blank())+cextheme+coord_fixed(ratio = 1.25,xlim=c(NA,-1250),ylim=c(200,650))
p2<-plotSpatialCov(sim[[2]],crop(expandUD[[2]]$outbound,slim),colors=viridis_pal()(100),alpha=0.5,col=stateCols[[4]],return=TRUE)+ggplot2::labs(fill = expression("log"(pi)),x=NULL,y=NULL)+guides(col = "none",shape="none")+theme(panel.border = element_blank())  + ggtitle("1 State")+ theme(plot.title = element_text(size=14,hjust=0.5))+cextheme+coord_fixed(ratio = 1.25,xlim=c(NA,-1250),ylim=c(200,650))
p3<-plotSpatialCov(sim[[9]],crop(expandUD[[9]]$foraging,slim),colors=viridis_pal()(100),alpha=0.5,col=stateCols[[4]],return=TRUE)+ggplot2::labs(fill = expression("log"(pi[2])),x="easting (km)",y="northing (km)")+guides(col = "none",shape="none")+theme(panel.border = element_blank())  + ggtitle("4 States")+ theme(plot.title = element_text(size=14,hjust=0.5))+cextheme+coord_fixed(ratio = 1.25,xlim=c(NA,-1250),ylim=c(200,650))
p4<-plotSpatialCov(sim[[14]],crop(expandUD[[14]]$foraging2,slim),colors=viridis_pal()(100),alpha=0.5,col=stateCols[[4]],return=TRUE)+ggplot2::labs(fill = expression("log"(pi[3])),x="easting (km)",y=NULL)+guides(col = "none",shape="none")+theme(panel.border = element_blank())  + ggtitle("5 States")+ theme(plot.title = element_text(size=14,hjust=0.5))+cextheme+coord_fixed(ratio = 1.25,xlim=c(NA,-1250),ylim=c(200,650))
pdf("simTracks.pdf",width=14,height=14)
plot(p1 + p2 + p3 + p4 + plot_layout(ncol = 2, nrow=2) )
dev.off()


# MALA simulations
mala <- list()
mala[[2]] <- malaSim(fitLangevin[[2]],inputUD[[2]],niter=50)
apply(mala[[2]] ,2,mean)
mala[[9]] <- malaSim(fitLangevin[[9]],inputUD[[9]],niter=50)
apply(mala[[9]] ,2,mean)
mala[[14]] <- malaSim(fitLangevin[[14]],inputUD[[14]],niter=50)
apply(mala[[14]] ,2,mean)
