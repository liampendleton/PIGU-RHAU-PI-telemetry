options("rgdal_show_exportToProj4_warnings"="none") # suppress annoying warnings
library(tidyverse)
library(raster)
library(rasterVis)
library(viridis)
library(ggplot2)
library(doFuture)
library(doRNG)
if(!requireNamespace("RandomFieldsUtils",quietly=TRUE)){
  install.packages("RandomFieldsUtils_1.2.5.tar.gz", repos = NULL, type = "source") # most recent archived version; required by RandomFields
}
if(!requireNamespace("RandomFields",quietly=TRUE)){
  install.packages("RandomFields_3.3.14.tar.gz", repos = NULL, type = "source") # most recent archived version; required by Rhabit
}
remotes::install_github("papayoun/Rhabit@31ddf44",dependencies = TRUE) # last commit before RandomFields was removed from dependencies
library(Rhabit)
if(!requireNamespace("momentuHMM",quietly=TRUE) || packageVersion("momentuHMM")<2){
  remotes::install_github("bmcclintock/momentuHMM@develop",dependencies = TRUE) # requires momentuHMM version >= 2.0.0
}
library(momentuHMM)

nsims <- 100
nbAnimals <- 20
lambda <- 100 # observation rate (1/lambda is expected time between successive observations)
obsPerAnimal <- 50000
ncores <- 7

lim <- c(-1, 1, -1, 1)*100
cropExtent <- extent(lim)
resol <- 1
ncov <- 3

beta1 <- c(6*resol,-4*resol,-5*resol,-0.1*resol) # state1 resource selection coefficients for the spatial covariates
sd_1 <- sqrt(5)
beta2 <- c(-4*resol,6*resol,5*resol,-0.1*resol) # state2 resource selection coefficients for the spatial covariates
sd_2 <- sqrt(7.5)

beta0 <- matrix(c(-2,-2),1,2)

fixPar2 <- list(mu=c(NA,1,2,3,4,NA,5,6,7,8,NA,1,2,3,4,NA,5,6,7,8,9,9,10,9,9,10,NA,NA))
fixPar1 <- list(mu=c(NA,1,2,3,4,NA,1,2,3,4,5,5,NA))

# subsample "high resolution" data
trimSeq <- obsPerAnimal/c(1,2,5,10,20,50,100)

set.seed(1,kind="Mersenne-Twister",normal.kind = "Inversion")

results <- list()
covlist <- list()
for(isim in 1:nsims){
  cat("Simulation",isim,"\n")
  #######################
  ## Define covariates ##
  #######################
  # Generate ncov spatial covariates
  covlist[[isim]] <- list()
  for(i in 1:ncov) {
    covlist[[isim]][[i]] <- Rhabit::simSpatialCov(lim = lim, nu = 0.6, rho = 50, sigma2 = 0.1, 
                                          resol = resol, raster_like = TRUE)
  }
  
  # Include squared distance to origin as covariate
  xgrid <- seq(lim[1], lim[2], by=resol)
  ygrid <- seq(lim[3], lim[4], by=resol)
  xygrid <- expand.grid(xgrid,ygrid)
  dist2 <- ((xygrid[,1])^2+(xygrid[,2])^2)/100
  covlist[[isim]][[4]] <- list(x=xgrid, y=ygrid, z=matrix(dist2, length(xgrid), length(ygrid)))
  
  # Compute utilization distribution for states 1 and 2
  UD1 <- Rhabit::getUD(covariates=covlist[[isim]], beta=beta1,log=TRUE)
  UD2 <- Rhabit::getUD(covariates=covlist[[isim]], beta=beta2,log=TRUE)
  
  # Plot covariates
  ggtheme <- theme(axis.title = element_text(size=12), axis.text = element_text(size=12),
                   legend.title = element_text(size=12), legend.text = element_text(size=12))
  c1plot <- Rhabit::plotRaster(rhabitToRaster(covlist[[isim]][[1]]), scale.name = expression(c[1])) + ggtheme
  c2plot <- Rhabit::plotRaster(rhabitToRaster(covlist[[isim]][[2]]), scale.name = expression(c[2])) + ggtheme
  UD1plot <- Rhabit::plotRaster(rhabitToRaster(UD1), scale.name = expression(pi)) + ggtheme
  UD2plot <- Rhabit::plotRaster(rhabitToRaster(UD2), scale.name = expression(pi)) + ggtheme
  
  names(covlist[[isim]]) <- c("cov1","cov2","cov3","d2c")
  spatialCovs <- lapply(lapply(covlist[[isim]],rhabitToRaster),function(x) {proj4string(x) <- CRS("+init=epsg:3416");return(x)})
  names(spatialCovs) <- c("cov1","cov2","cov3","d2c")
  for(i in names(spatialCovs)){
    names(spatialCovs[[i]]) <- i
  }
  
  dist <- list(mu="rw_mvnorm2")   # bivariate normal random walk
  
  # specify 2-state langevin pseudo-design matrix
  DM2 <- list(mu=list(mean.x=~0+mu.x_tm1+langevin(cov1.x)+langevin(cov2.x)+langevin(cov3.x)+langevin(d2c.x),
                      mean.y=~0+mu.y_tm1+langevin(cov1.y)+langevin(cov2.y)+langevin(cov3.y)+langevin(d2c.y),
                      sd.x=~state2(intercept),
                      sd.y=~state2(intercept),
                      corr.xy=~1))
  
  # simulate "high resolution" tracks; this can take a while...
  origData <- tryCatch(stop(),error=function(e) e)
  while(inherits(origData,"error")){
    initState <- sample.int(2,nbAnimals,replace=TRUE)
    initPos <- matrix(c(sample(ncell(rhabitToRaster(UD1)),nbAnimals,replace=FALSE,prob=exp(getValues(rhabitToRaster(UD1)))/sum(exp(getValues(rhabitToRaster(UD1))))),
                        sample(ncell(rhabitToRaster(UD2)),nbAnimals,replace=FALSE,prob=exp(getValues(rhabitToRaster(UD2)))/sum(exp(getValues(rhabitToRaster(UD2)))))),2,nbAnimals,byrow=TRUE)
    initialPosition <- mapply(function(x) xyFromCell(rhabitToRaster(UD1),initPos[initState[x],x]),1:nbAnimals,SIMPLIFY = FALSE)
    origData <- tryCatch(simCTHMM(nbAnimals=nbAnimals, obsPerAnimal=obsPerAnimal, ncores=ncores,
                                  nbStates=2,
                                  initialPosition=initialPosition,
                                  dist=dist["mu"],
                                  DM=DM2,
                                  Par=list(mu=c(1,beta1,1,beta2,1,beta1,1,beta2,log(sd_1*resol),log(sd_1*resol),log(sd_2*resol)-log(sd_1*resol),log(sd_1*resol),log(sd_1*resol),log(sd_2*resol)-log(sd_1*resol),0,0)),#,aux=c(10,2)),
                                  beta=beta0,
                                  formulaDelta=~0+ID,
                                  delta = matrix(ifelse(initState==1,-1.e+100,1.e+100),nbAnimals,1),
                                  spatialCovs = spatialCovs, 
                                  covs=data.frame(intercept=1),
                                  gradient=TRUE,
                                  mvnCoords = "mu",
                                  lambda = lambda, # observation rate (1/lambda is expected time between successive observations); state switches can only occur at times of observations
                                  states=TRUE,
                                  TMB=TRUE),error=function(e) e) 
  }
  
  DM1 <-  list(mu=list(mean.x=~0+mu.x_tm1+langevin(cov1.x)+langevin(cov2.x)+langevin(cov3.x)+langevin(d2c.x),
                       mean.y=~0+mu.y_tm1+langevin(cov1.y)+langevin(cov2.y)+langevin(cov3.y)+langevin(d2c.y),
                       sd.x=~1,
                       sd.y=~1,
                       corr.xy=~1))
  
  langData <- list()
  for(j in seq_along(trimSeq)){
    trackSamp <- numeric()
    for(i in 1:nbAnimals){
      trackSamp <- c(trackSamp,c(0,cumsum(table(origData$ID)))[i]+c(1,sort(sample.int(obsPerAnimal-1,trimSeq[j]-1,replace=FALSE)+1)))
    }
    langData[[j]] <- origData[trackSamp,]
  }
  
  nbStates <- list(2,1)
  DM <- list(DM2,DM1)
  Par0 <- list(list(mu=c(1,beta1,1,beta2,1,beta1,1,beta2,log(sd_1*resol),log(sd_1*resol),log(sd_2*resol)-log(sd_1*resol),log(sd_1*resol),log(sd_1*resol),log(sd_2*resol)-log(sd_1*resol),0,0)),
               list(mu=c(1,apply(rbind(beta1,beta2),2,mean),1,apply(rbind(beta1,beta2),2,mean),log((sd_1+sd_2)/2*resol),log((sd_1+sd_2)/2*resol),0)))
  fbeta0 <- list(beta0,NULL)
  fixPar <- list(fixPar2,fixPar1)
  
  modNames <- list("fit2","fit1")
  
  cat("Fitting models to subsampled data in parallel...")
  doFuture::registerDoFuture()
  future::plan(future::multisession,workers=ncores)
  results[[isim]] <- foreach(j=seq_along(trimSeq),lD=langData) %:% 
    foreach(i=1:2,.final = function(x) setNames(x, modNames)) %dopar% {
      fit <- suppressMessages(fitCTHMM(lD,nbStates=nbStates[[i]],
                                       dist=dist["mu"],
                                       DM=DM[[i]],
                                       Par0=Par0[[i]],
                                       beta0=fbeta0[[i]],
                                       fixPar=fixPar[[i]],
                                       mvnCoords = "mu",
                                       optMethod="TMB",
                                       control=list(silent=TRUE)))
      
      ret <- list()
      ret$CIbeta <- fit$CIbeta
      ret$CIreal <- fit$CIreal
      if(nbStates[[i]]>1){
        ret$states <- fit$data$states
        ret$viterbi <- viterbi(fit)
        ret$stateProbs <- stateProbs(fit)
      }
      return(ret)
    }
  future::plan(future::sequential)
  cat("Done\n\n")
}

# summarize output
estimates <- list()
for(j in names(results[[1]][[1]])){
  estimates[[j]] <- NULL
  if(j=="fit2") muInd <- c(1,2,3,4,5,7,8,9,10,21,23,27)
  else if(j=="fit1") muInd <- c(1,2,3,4,5,11,13)
  for(k in as.character(obsPerAnimal/trimSeq[1:7]*1/lambda)){
    tmp <- as.data.frame(do.call(rbind,lapply(results,function(x) x[[which(as.character(obsPerAnimal/trimSeq*1/lambda)==k)]][[j]]$CIbeta$mu$est[muInd])))
    tmp$seq <- rep((obsPerAnimal/trimSeq*1/lambda)[which(as.character(obsPerAnimal/trimSeq*1/lambda)==k)],nsims)
    if(grepl("fit2",j)) {
      tmp$sigma2_1 <- do.call(rbind,lapply(results,function(x) x[[which(as.character(obsPerAnimal/trimSeq*1/lambda)==k)]][[j]]$CIreal$mu$est[3,1]^2))
      tmp$sigma2_2 <- do.call(rbind,lapply(results,function(x) x[[which(as.character(obsPerAnimal/trimSeq*1/lambda)==k)]][[j]]$CIreal$mu$est[3,2]^2))
      tmp$q_12 <- do.call(rbind,lapply(results,function(x) x[[which(as.character(obsPerAnimal/trimSeq*1/lambda)==k)]][[j]]$CIbeta$beta$est[1,1]))
      tmp$q_21 <- do.call(rbind,lapply(results,function(x) x[[which(as.character(obsPerAnimal/trimSeq*1/lambda)==k)]][[j]]$CIbeta$beta$est[1,2]))
      tmp$viterbi <- do.call(rbind,lapply(results,function(x) mean(x[[which(as.character(obsPerAnimal/trimSeq*1/lambda)==k)]][[j]]$viterbi==x[[which(as.character(obsPerAnimal/trimSeq*1/lambda)==k)]][[j]]$states)))
      tmp$stateProbs <- do.call(rbind,lapply(results,function(x) {
        stateProbs <- x[[which(as.character(obsPerAnimal/trimSeq*1/lambda)==k)]][[j]]$stateProbs;
        states <- x[[which(as.character(obsPerAnimal/trimSeq*1/lambda)==k)]][[j]]$states;
        mean(mapply(function(y) stateProbs[y,states[y]],1:length(states))>0.5)
      }))
      colnames(tmp)[c(2:9,14:15,16:17,18:19)] <- c("beta[11]","beta[12]","beta[13]","beta[14]","beta[21]","beta[22]","beta[23]","beta[24]","sigma[1]^2","sigma[2]^2","log(q[12])","log(q[21])","Viterbi","p[s]>0.5")
    } else {
      tmp$sigma2 <- do.call(rbind,lapply(results,function(x) x[[which(as.character(obsPerAnimal/trimSeq*1/lambda)==k)]][[j]]$CIreal$mu$est[3,1]^2))
      colnames(tmp)[c(2:5,9)] <- c("beta[1]","beta[2]","beta[3]","beta[4]","sigma^2")
    }
    estimates[[j]] <- rbind(estimates[[j]],tmp)
    
  }
}

library(ggplot2)
# Prepare data set of estimates
parnames <- truepar <- zeroline <- list()
parnames[[1]] <- c("beta[11]","beta[12]","beta[13]","beta[14]","beta[21]","beta[22]","beta[23]","beta[24]","sigma[1]^2","sigma[2]^2","log(q[12])","log(q[21])")#,"viterbi","stateProbs")
parnames[[2]] <- c("beta[1]","beta[2]","beta[3]","beta[4]","sigma^2")
truepar[[1]] <- data.frame(par = parnames[[1]][1:12], value = c(beta1, beta2,sd_1^2,sd_2^2,beta0[1,1],beta0[1,2]), cols=c(rep(2,4),rep(4,4),2,4,2,4))
truepar[[2]] <- data.frame(par = rep(parnames[[2]],2), value = c(beta1, sd_1^2,beta2,sd_2^2), cols=c(rep(2,5),rep(4,5)))
zeroline[[1]] <- data.frame(yint = c(rep(0,8)), par=parnames[[1]][c(1:8)])
zeroline[[2]] <- data.frame(yint = 0, par=parnames[[2]][1:4])

p <- list()
for(i in 1:length(estimates)) {
  df <- estimates[[i]][,parnames[[i]]]
  colnames(df) <- parnames[[i]]
  df2 <- reshape2::melt(df)
  colnames(df2) <- c("par", "value")
  plotdat <- cbind(dt = estimates[[i]]$seq, df2)
  plotdat$dt <- as.factor(plotdat$dt)
  
  p[[i]] <- ggplot(plotdat, aes(dt, value)) + geom_hline(data=zeroline[[i]], aes(yintercept=yint), lty=2) + 
    geom_hline(data=truepar[[i]], mapping=aes(yintercept=value), lty=2, col=truepar[[i]]$cols) + 
    facet_wrap("par", scale = "free", labeller = label_parsed) + geom_boxplot(outlier.size = 0.6) +
    theme(strip.text = element_text(size = 12), axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)) + xlab("Time interval") + ylab("Estimated value")
  
}

pdf("sim2results.pdf",width=13,height=7)
p[[1]]
dev.off()

pdf("sim1results.pdf",width=10,height=5)
p[[2]]
dev.off()

parnames <- zeroline <- list()
parnames[[1]] <- c("Viterbi","p[s]>0.5")
zeroline[[1]] <-   data.frame(yint = c(0.95,0.95), par=parnames[[1]])

df <- estimates[[1]][,parnames[[1]]]
colnames(df) <- parnames[[1]]
df2 <- reshape2::melt(df)
colnames(df2) <- c("par", "value")
plotdat <- cbind(dt = estimates[[1]]$seq, df2)
plotdat$dt <- as.factor(plotdat$dt)

sp <- ggplot(plotdat, aes(dt, value)) + geom_hline(data=zeroline[[1]], aes(yintercept=yint), lty=2) + 
  facet_wrap("par", scale = "free", labeller = label_parsed) + geom_boxplot(outlier.size = 0.6) +
  theme(strip.text = element_text(size = 12), axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) + xlab("Time interval") + ylab("Proportion")

pdf("stateResults.pdf",width=8.5,height=4.25)
sp
dev.off()

covsim <- sample.int(nsims,1)

UD1 <- rhabitToRaster(getUD(covariates=covlist[[covsim]], beta=beta1,log=TRUE))
names(UD1) <- "UD1"
UD2 <- rhabitToRaster(getUD(covariates=covlist[[covsim]], beta=beta2,log=TRUE))
names(UD2) <- "UD2"
estUD <- list()
estUD[[1]] <- UD1
estUD[[2]] <- UD2
estUD[[3]] <- rhabitToRaster(getUD(covariates=covlist[[covsim]], beta=results[[covsim]][[4]]$fit2$CIbeta$mu$est[2:5],log=TRUE))
estUD[[4]] <- rhabitToRaster(getUD(covariates=covlist[[covsim]], beta=results[[covsim]][[4]]$fit2$CIbeta$mu$est[7:10],log=TRUE))
estUD[[5]] <- rhabitToRaster(getUD(covariates=covlist[[covsim]], beta=results[[covsim]][[1]]$fit2$CIbeta$mu$est[2:5],log=TRUE))
estUD[[6]] <- rhabitToRaster(getUD(covariates=covlist[[covsim]], beta=results[[covsim]][[1]]$fit2$CIbeta$mu$est[7:10],log=TRUE))
estUD[[7]] <- rhabitToRaster(getUD(covariates=covlist[[covsim]], beta=results[[covsim]][[5]]$fit2$CIbeta$mu$est[2:5],log=TRUE))
estUD[[8]] <- rhabitToRaster(getUD(covariates=covlist[[covsim]], beta=results[[covsim]][[5]]$fit2$CIbeta$mu$est[7:10],log=TRUE))
estUD[[9]] <- rhabitToRaster(getUD(covariates=covlist[[covsim]], beta=results[[covsim]][[2]]$fit2$CIbeta$mu$est[2:5],log=TRUE))
estUD[[10]] <- rhabitToRaster(getUD(covariates=covlist[[covsim]], beta=results[[covsim]][[2]]$fit2$CIbeta$mu$est[7:10],log=TRUE))
estUD[[11]] <- rhabitToRaster(getUD(covariates=covlist[[covsim]], beta=results[[covsim]][[6]]$fit2$CIbeta$mu$est[2:5],log=TRUE))
estUD[[12]] <- rhabitToRaster(getUD(covariates=covlist[[covsim]], beta=results[[covsim]][[6]]$fit2$CIbeta$mu$est[7:10],log=TRUE))
estUD[[13]] <- rhabitToRaster(getUD(covariates=covlist[[covsim]], beta=results[[covsim]][[3]]$fit2$CIbeta$mu$est[2:5],log=TRUE))
estUD[[14]] <- rhabitToRaster(getUD(covariates=covlist[[covsim]], beta=results[[covsim]][[3]]$fit2$CIbeta$mu$est[7:10],log=TRUE))
estUD[[15]] <- rhabitToRaster(getUD(covariates=covlist[[covsim]], beta=results[[covsim]][[7]]$fit2$CIbeta$mu$est[2:5],log=TRUE))
estUD[[16]] <- rhabitToRaster(getUD(covariates=covlist[[covsim]], beta=results[[covsim]][[7]]$fit2$CIbeta$mu$est[7:10],log=TRUE))
estUD <- stack(estUD)
UDnames <- c("True UD1","True UD2","UD1: \u0394\u0074 = 0.1","UD2: \u0394\u0074 = 0.1","UD1: \u0394\u0074 = 0.01","UD2: \u0394\u0074 = 0.01","UD1: \u0394\u0074 = 0.2","UD2: \u0394\u0074 = 0.2","UD1: \u0394\u0074 = 0.02","UD2: \u0394\u0074 = 0.02","UD1: \u0394\u0074 = 0.5","UD2: \u0394\u0074 = 0.5","UD1: \u0394\u0074 = 0.05","UD2: \u0394\u0074 = 0.05","UD1: \u0394\u0074 = 1","UD2: \u0394\u0074 = 1")
names(estUD) <- UDnames

library(Cairo)
Cairo(type="pdf",file="simUD2.pdf",width=11,height=11, units='in')
spplot(estUD,
       layout=c(4, 4), # create a 4x4 layout for the data
       col.regions=viridis_pal()(5000), # add a color ramp
       useRaster=TRUE,
       xlab="easting",ylab="northing",
       names.attr=UDnames,at=seq(min(getValues(estUD)),max(getValues(estUD)),length=5000),colorkey=list(title=expression(log(pi)),title.control=list(side="right")))
dev.off()


estUD1 <- list()
estUD1[[1]] <- rhabitToRaster(getUD(covariates=covlist[[covsim]], beta=results[[covsim]][[1]]$fit1$CIbeta$mu$est[2:5],log=TRUE))
estUD1[[2]] <- rhabitToRaster(getUD(covariates=covlist[[covsim]], beta=results[[covsim]][[2]]$fit1$CIbeta$mu$est[2:5],log=TRUE))
estUD1[[3]] <- rhabitToRaster(getUD(covariates=covlist[[covsim]], beta=results[[covsim]][[3]]$fit1$CIbeta$mu$est[2:5],log=TRUE))
estUD1[[4]] <- rhabitToRaster(getUD(covariates=covlist[[covsim]], beta=results[[covsim]][[4]]$fit1$CIbeta$mu$est[2:5],log=TRUE))
estUD1[[5]] <- rhabitToRaster(getUD(covariates=covlist[[covsim]], beta=results[[covsim]][[5]]$fit1$CIbeta$mu$est[2:5],log=TRUE))
estUD1[[6]] <- rhabitToRaster(getUD(covariates=covlist[[covsim]], beta=results[[covsim]][[6]]$fit1$CIbeta$mu$est[2:5],log=TRUE))
estUD1[[7]] <- rhabitToRaster(getUD(covariates=covlist[[covsim]], beta=results[[covsim]][[7]]$fit1$CIbeta$mu$est[2:5],log=TRUE))
estUD1 <- stack(estUD1)
UDnames <- c("\u0394\u0074 = 0.01","\u0394\u0074 = 0.02","\u0394\u0074 = 0.05","\u0394\u0074 = 0.1","\u0394\u0074 = 0.2","\u0394\u0074 = 0.5","\u0394\u0074 = 1")
names(estUD1) <- UDnames

Cairo(type="pdf",file="simUD1.pdf",width=11,height=7, units='in')
spplot(estUD1,
       layout=c(4, 2), # create a 4x4 layout for the data
       col.regions=viridis_pal()(5000), # add a color ramp
       useRaster=TRUE,
       xlab="easting",ylab="northing",
       names.attr=UDnames,at=seq(min(getValues(estUD1)),max(getValues(estUD1)),length=5000),colorkey=list(title=expression(log(pi)),title.control=list(side="right")))
dev.off()