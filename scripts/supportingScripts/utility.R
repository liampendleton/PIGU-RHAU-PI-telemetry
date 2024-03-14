stateCols <- list()
stateCols[[1]] <- "#E69F00"
stateCols[[2]] <- c("#E69F00","#009E73")
stateCols[[3]] <- c("#E69F00", "#D55E00", "#009E73")
stateCols[[4]] <- c("#E69F00", "#56B4E9", "#D55E00", "#009E73" )
stateCols[[5]] <- c("#E69F00", "#0072B2", "#56B4E9", "#D55E00", "#009E73" )

# function to plot UDs
plotUD <- function(nbUD,parIndex,covNames,model,UDnames,UDstates,sign=NULL,cropRast=TRUE,return=FALSE){
  stateNames <- model$stateNames
  nbStates <- length(stateNames)
  par(mfrow=c(nbUD,1))
  if(length(model$stateNames)>1){
    states <- viterbi(model)
  } else states <- rep(1,nrow(model$data))
  betas <- covs <- logUD <- list()
  for(i in 1:nbUD){
    betas[[i]] <- model$CIbeta$mu$est[parIndex[[i]]]
    if(!is.null(sign[[i]])) betas[[i]] <- betas[[i]] * sign[[i]]
    if(cropRast) covs[[i]] <- covlistCrop[covNames[[i]]]
    else covs[[i]] <- covlist0[covNames[[i]]]
    logUD[[UDnames[i]]] <- betas[[i]][1]*covs[[i]][[1]]
    for(j in 2:length(betas[[i]])){
      logUD[[UDnames[i]]] <- logUD[[UDnames[i]]] + betas[[i]][j]*covs[[i]][[j]]
    }
    values(logUD[[UDnames[i]]]) <- getValues(logUD[[UDnames[i]]]) - log(sum(Brobdingnag::brob(getValues(logUD[[UDnames[i]]])[!is.na(getValues(logUD[[UDnames[i]]]))])))
    names(logUD[[UDnames[i]]]) <- UDnames[i]
    raster::plot(logUD[[UDnames[i]]],col=viridis_pal()(100),main=paste0(model$modelName,": ",UDnames[i]))
    points(model$data[which(states %in% UDstates[[i]]),c("mu.x_tm1","mu.y_tm1")],col=stateCols[[nbStates]][i])
  }
  par(mfrow=c(1,1))
  if(return) return(logUD)
}

getWorkBounds <- function(DM,nbStates,data){
  workBounds <- list()
  for(i in names(DM)){
    if(i=="mu") parNames <- list(mu=c("mean.x","mean.y","sd.x","sd.y","corr.xy"))
    else parNames <- list(dry="prob")
    fullDM <- momentuHMM:::make_matrices(momentuHMM:::make_formulas(DM[i],i,parNames,nbStates),data %>% mutate(mu.x_tm1=mu.x,mu.y_tm1=mu.y))$X_fe
    parms <- colnames(fullDM)
    workBounds[[i]] <- matrix(rep(c(-Inf,Inf),ncol(fullDM)),ncol=2,byrow=TRUE,dimnames=list(parms,c("lower","upper")))
    if(i=="mu"){
      bInd <- which(grepl(".intercept",colnames(fullDM)))
      if(nbStates>2 & length(bInd)) workBounds[[i]][bInd,] <- matrix(c(-Inf,0),nrow=length(bInd),ncol=2,byrow=TRUE)
    } else {
      if(nbStates>1) {
        workBounds[[i]][1:(nbStates-1),] <- matrix(c(-Inf,0),nrow=nbStates-1,ncol=2,byrow=TRUE)
        workBounds[[i]][nbStates,] <- c(0,Inf)
      } else workBounds[[i]][nbStates,] <- c(-Inf,0)
    }
  }
  return(workBounds)
}

getFixPar <- function(DM,nbStates,data,meancons=NULL,sdcons=NULL,drycons=NULL){
  fixPar <- list()
  parNames <- list(mu=c("mean.x","mean.y","sd.x","sd.y","corr.xy"))
  fullDM <- momentuHMM:::make_matrices(momentuHMM:::make_formulas(DM["mu"],"mu",parNames,nbStates),data %>% mutate(mu.x_tm1=mu.x,mu.y_tm1=mu.y))$X_fe
  parms <- colnames(fullDM)
  fixPar$mu <- 1:length(parms)
  tm1Ind <- which(grepl("_tm1",parms))
  fixPar$mu[tm1Ind] <- NA
  if(nbStates>2){
    if(!is.null(meancons)){
      for(j in 1:length(meancons)){
        refInd <- which(grepl(paste0("mean.x.state",meancons[[j]][1]),parms))
        conInd <- NULL
        for(k in 2:length(meancons[[j]])){
          conInd <- c(conInd,which(grepl(paste0("mean.x.state",meancons[[j]][k]),parms)))
        }
        fixPar$mu[conInd] <- fixPar$mu[refInd]
      }
    }
  }
  meanInd <- which(grepl("mean",parms))
  fixPar$mu[which(duplicated(gsub(".y","",gsub(".x","",parms[meanInd]))))] <- fixPar$mu[which(!duplicated(gsub(".y","",gsub(".x","",parms[meanInd]))))]
  
  sdInd <- which(grepl("sd",parms))
  if(nbStates>1){
    fixPar$mu[which(grepl(paste0("mu.sd.x.state",nbStates,".(Intercept)"),parms,fixed=TRUE))] <- fixPar$mu[sdInd[1]]
    if(!is.null(sdcons)){
      for(j in 1:length(sdcons)){
        refInd <- which(grepl(paste0("sd.x.state",sdcons[[j]][1]),parms))
        conInd <- NULL
        for(k in 2:length(sdcons[[j]])){
          conInd <- c(conInd,which(grepl(paste0("sd.x.state",sdcons[[j]][k]),parms)))
        }
        fixPar$mu[conInd] <- fixPar$mu[refInd]
      }
    }
  }
  sdInd <- which(grepl("sd",parms))
  fixPar$mu[length(meanInd)+which(duplicated(gsub(".y.","",gsub(".x.","",parms[sdInd]))))] <- fixPar$mu[length(meanInd)+which(!duplicated(gsub(".y.","",gsub(".x.","",parms[sdInd]))))]
  
  fixPar$mu <- momentuHMM:::reorderFix(fixPar$mu)
  fixPar$mu[which(grepl("corr.xy",parms))] <- NA
  names(fixPar$mu) <- parms
  
  parNames <- list(dry="prob")
  fullDM <- momentuHMM:::make_matrices(momentuHMM:::make_formulas(DM["dry"],"dry",parNames,nbStates),data)$X_fe
  parms <- colnames(fullDM)
  fixPar$dry <- 1:ncol(fullDM)
  if(!is.null(drycons)){
    for(j in 1:length(drycons)){
      refInd <- which(grepl(paste0("prob.state",drycons[[j]][1]),parms))
      conInd <- NULL
      for(k in 2:length(drycons[[j]])){
        conInd <- c(conInd,which(grepl(paste0("prob.state",drycons[[j]][k]),parms)))
      }
      fixPar$dry[conInd] <- fixPar$dry[refInd]
    }
    fixPar$dry <- momentuHMM:::reorderFix(fixPar$dry)
  }
  
  if(nbStates>2){
    fixPar$beta <- matrix(1:(2*nbStates*(nbStates-1)),2,nbStates*(nbStates-1))
    fixPar$beta[,nbStates-1] <- NA
    fixPar$beta[,nbStates*(nbStates-1)] <- NA
    fixPar$beta[2,nbStates*(nbStates-1)-(nbStates-2):1] <- NA
    fixPar$beta <- momentuHMM:::reorderFix(fixPar$beta)
  }
  fixPar$delta <- c(rep(NA,nbStates-1),NA*(nbStates-1))
  return(fixPar)
}

getLangIndex <- function(fp,nbStates){
  estParNames <- names(fp$mu)[which(!is.na(fp$mu) & !duplicated(fp$mu))]
  parIndex <- list()
  if(any(grepl("state1.langevin",estParNames))) parIndex[[1]] <- match(estParNames[which(grepl("state1.langevin",estParNames))],names(fp$mu))
  if(nbStates>2){
    for(i in 2:(nbStates-1)){
      if(any(grepl(paste0("state",i,".langevin"),estParNames))) parIndex[[i]] <- match(estParNames[which(grepl(paste0("state",i,".langevin"),estParNames))],names(fp$mu))
    }
  }
  return(parIndex)
}

getPrior <- function(fixPar,sd=10){
  if(!is.null(fixPar$beta)){
    parInd <- which(is.na(fixPar$beta))
    nPar <- sum(unlist(lapply(fixPar[c("mu","dry")],function(x) length(x))))
    prior <- list(beta=matrix(c(0,sd),nrow=length(fixPar$beta),ncol=2,byrow=TRUE))
    prior$beta[parInd,] <- NA
    return(prior)
  } else return()
}

penAIC <- function(models,sampleSize=3,penalty=1/2*100^(-1)){
  paic <- aic <- nll <- logpenalty <- nbPar <- numeric(length(models))
  modelName <- unlist(lapply(models,function(x) x$modelName))
  d <- penalty/sampleSize
  for(i in 1:length(models)){
    nPar <- sum(unlist(lapply(models[[i]]$conditions$fullDM,function(x) ncol(x))))
    if(isTRUE(models[[i]]$conditions$optMethod=="TMB")) parInd <- nPar+which(!is.na(models[[i]]$conditions$fixPar$beta))
    else parInd <- nPar+which(is.na(models[[i]]$conditions$fixPar$beta))
    F <- matrix(0,length(models[[i]]$mod$estimate),length(models[[i]]$mod$estimate))
    F[parInd,parInd] <- diag(length(parInd))
    R <- F %*% t(F)
    q <- Matrix::rankMatrix(R)
    I <- models[[i]]$mod$Sigma
    Tmat <- t(F) %*% I %*% F
    T <- eigen(Tmat)$values[1:q]
    logpenalty[i] <- ifelse(is.null(models[[i]]$prior),0,models[[i]]$prior(models[[i]]$mod$estimate))
    nll[i] <- models[[i]]$mod$minimum + logpenalty[i]
    maxLogLike <- -nll[i]
    nbPar[i] <- length(models[[i]]$mod$wpar)
    paic[i] <- -2*maxLogLike + 2*(nbPar[i] - sum(2*d*T / (1+2*d*T)))
    aic[i] <- -2*maxLogLike + 2*nbPar[i]
  }
  ord <- order(aic)
  weight <- exp(-0.5*(aic[ord]-min(aic)))/sum(exp(-0.5*(aic[ord]-min(aic))))
  delta <- aic[ord]-min(aic)
  pweight <- exp(-0.5*(paic[ord]-min(paic)))/sum(exp(-0.5*(paic[ord]-min(paic))))
  pdelta <- paic[ord]-min(paic)
  return(data.frame(Model=modelName[ord],AIC=aic[ord],Weight=weight,deltaAIC=delta,pAIC=paic[ord],pWeight=pweight,deltapAIC=pdelta,NLL=nll[ord],Penalty=logpenalty[ord],K=nbPar[ord]))
}

# plot mu pseudo-residuals
plotRes <- function(model){
  pr <- pseudoRes(model)$muRes
  if(any(is.infinite(pr))){
    warning("some pseudo-residuals for ","mu"," are infinite; these will be set to NA for plotting")
    pr[which(is.infinite(pr))] <- NA
  }
  par(mfrow=c(1,3))
  ylim <- range(unlist(pr[which(is.finite(pr))]))
  qq <- stats::qqnorm(pr,plot=FALSE)
  limInf <- min(min(qq$x,na.rm=TRUE),min(qq$y,na.rm=TRUE))
  limSup <- max(max(qq$x,na.rm=TRUE),max(qq$y,na.rm=TRUE))
  plot(pr,type="l",xlab="Observation index",ylab=expression(paste(mu," pseudo-residuals")),main="",ylim=ylim)
  stats::qqnorm(pr,main=model$modelName,col="red",xlim=c(limInf,limSup),ylim=c(limInf,limSup))
  abline(0,1,lwd=2)
  ACF <- acf(pr,lag.max=40,na.action=na.pass,plot=FALSE)
  acf(pr,lag.max=40,na.action=na.pass,main="",ylim=c(0,1))
  par(mfrow=c(1,1))
}

# calculate and plot predicted and residuals
plotPredRes <- function(model,UDinputs){
  parmInd <- UDinputs$parIndex
  covInd <- UDinputs$covNames
  covInd[[length(covInd)+1]] <- c("depth","slope","depthslope","d2site")
  
  nbStates <- length(model$stateNames)
  nbAnimals <- length(unique(model$data$ID))
  
  if(nbStates>1) states <- stateProbs(model)
  else states <- matrix(1,nrow=nrow(model$data),ncol=)
  model$data$pred.x <- NA
  model$data$pred.y <- NA
  sigma2 <- matrix(NA,nbAnimals,nbStates)
  
  cat("Calculating predicted steps and residuals for model '",model$modelName,"'...",sep="")
  for(zoo in 1:nbAnimals){
    sigma2[zoo,] <-  (CIreal(model,covs=data.frame(ID=unique(tracks$ID)[zoo]),parms="mu")$mu$est[3,])^2
    for(i in which(model$data$ID==unique(tracks$ID)[zoo])){
      model$data$pred.x[i] <- model$data$mu.x_tm1[i]
      model$data$pred.y[i] <- model$data$mu.y_tm1[i]
      for(st in 1:nbStates){
        if(nbStates==1 || st < nbStates){
          beta <- model$CIbeta$mu$est[parmInd[[st]]]
        } else {
          beta <- rep(0,4)
        }
        model$data$pred.x[i] <- model$data$pred.x[i] + states[i,st]*sigma2[zoo,st]*model$data$dt[i]/2 * sum(beta * model$data[i,paste0(covInd[[st]],".x")])
        model$data$pred.y[i] <- model$data$pred.y[i] + states[i,st]*sigma2[zoo,st]*model$data$dt[i]/2 * sum(beta * model$data[i,paste0(covInd[[st]],".y")])
      }
    }
  }
  cat("Done\n\n")
  
  model$data$res.x <- model$data$mu.x - model$data$pred.x
  model$data$res.y <- model$data$mu.y - model$data$pred.y
  model$data$res2 <- sqrt(model$data$res.x^2 + model$data$res.y^2)
  
  my_x_zoom <- c(-1800, -1700)
  my_y_zoom <- c(550, 656)
  main_plot <- ggplot(model$data, aes(x = mu.x, y = mu.y)) +
    geom_path(aes(x=mu.x,y=mu.y,group=ID),size = 1, col = "red", alpha = 0.2) +
    geom_segment(aes(xend = mu.x,
                     yend = mu.y,
                     group=ID),
                 size = 1, col = "red", alpha = 0.1) + 
    geom_segment(aes(x=mu.x_tm1,
                     y=mu.y_tm1,
                     xend = pred.x,
                     yend = pred.y,
                     group=ID,
                     col = res2),
                 size = 1,
                 alpha = 0.5,
                 arrow = ggplot2::arrow(length = unit(0.3, "cm"))) +
    labs(col = "Residual") + 
    theme_bw() 
  
  p1 <- main_plot + annotate("rect", 
             xmin = my_x_zoom[1], xmax = my_x_zoom[2],
             ymin = my_y_zoom[1], ymax = my_y_zoom[2],
             alpha = 0.5, fill = "blue")
  
  
  p2 <- main_plot +
    coord_cartesian(xlim = my_x_zoom,
                    ylim = my_y_zoom)
  
  cat("Summary:\n")
  print(summary(model$data$res2))
  cat("\n")
  cat("Plotting histogram...")
  hist(model$data$res2,main=model$modelName,xlab="Residual norm",ylim=c(0,6000),xlim=c(0,140),breaks=seq(0,140,5))
  cat("Done\n")
  cat("Plotting predictions...")
  cat("Done\n")
  p1 + p2 + plot_layout(ncol = 2, guides = "collect") + plot_annotation(title = model$modelName)
}

tpmPlot <- function(model,distQuantile=0.99){
  
  Pmat <- list()
  d2site <- seq(min(model$data$d2site),quantile(model$data$d2site,distQuantile),length=200)
  cat("Calculating TPM and confidence intervals...")
  for(i in 1:length(d2site)){
    Pmat[[i]] <- CIreal(model,covs=data.frame(d2site=d2site[i],dt=median(model$data$dt)),parms="gamma")$gamma
  }
  cat("Done\n")
  
  d2mean <- mean(covlist$d2site[cellInd])
  d2sd <- sd(covlist$d2site[cellInd])
  d2site <- d2site * d2sd+d2mean
  nbStates <- length(model$stateNames)
  st <- viterbi(model)
  cols <- stateCols[[nbStates]]
  
  par(mfrow=c(nbStates,nbStates))
  par(mar=c(5,4.5,4,2)-c(0,0,1.5,1))
  for(i in 1:nbStates){
    for(j in 1:nbStates){
      std2site <- unlist(mapply(function(x) {
          ind <- which(model$data$ID==x);
          model$data$d2site[ind][which(st[ind][-length(ind)]==i & st[ind][-1]==j & model$data$d2site[ind[-length(ind)]]<=quantile(model$data$d2site,distQuantile))];
        },unique(model$data$ID))) * d2sd+d2mean
      counts <- hist(std2site,breaks=d2site,plot=FALSE,right=FALSE)$counts
      ind <- which(d2site<=max(std2site) & d2site>=min(std2site))
      stPmat <- unlist(lapply(Pmat,function(x) x$est[i,j]))[ind]
      
      plot(d2site,unlist(lapply(Pmat,function(x) x$est[i,j])),type="l",xlim=range(d2site),ylim=c(0,1),xlab="d2site",ylab=paste(i,"->",j),cex.lab=1.5,cex.axis=1.3)
      #lines(d2site[ind],unlist(lapply(Pmat,function(x) x$est[i,j]))[ind],type="l",col=cols[i])
      segments(d2site[ind[-length(ind)]],stPmat[-length(ind)],d2site[ind[-1]],stPmat[-1],xlim=range(d2site),col=cols[i],lwd=5*counts[ind]/max(counts[ind],na.rm=TRUE)+1)
      lines(d2site,unlist(lapply(Pmat,function(x) x$lower[i,j])),lty=3)
      lines(d2site,unlist(lapply(Pmat,function(x) x$upper[i,j])),lty=3)
    }
  }
}

statPlot <- function(model,distQuantile=0.99){
  nbStates <- length(model$stateNames)
  stat <- plotStationary(model,plotCI=TRUE,col=stateCols[[nbStates]],return=TRUE)
  stat$d2site <- do.call(rbind,stat$d2site)
  stat$d2site$state <- rep(model$stateNames,each=nrow(stat$d2site)/nbStates)
  stat$d2site <- stat$d2site[which(stat$d2site$cov<=quantile(model$data$d2site,distQuantile)),]
  d2mean <- mean(covlist$d2site[cellInd])
  d2sd <- sd(covlist$d2site[cellInd])
  stat$d2site$cov <- stat$d2site$cov * d2sd+d2mean
  size1 <- rel(1.15)
  size2 <- rel(1.25)
  size3 <- rel(1.35)
  ggplot(stat$d2site,aes(x=cov,y=est,color=factor(state,levels=model$stateNames),fill=factor(state,levels=model$stateNames)))+
    geom_line(aes(color=factor(state,levels=model$stateNames)))+
    scale_color_manual(values=stateCols[[nbStates]])+
    labs(fill = "state", col = "state",x="d2site", y = "Stationary state probabilities")+
    geom_ribbon(aes(ymin = lci, ymax = uci,group=factor(state,levels=model$stateNames)),alpha = 0.2) +
    scale_fill_manual(values=stateCols[[nbStates]])+
    theme(legend.text = element_text(size = size1),legend.title=element_text(size=size3),axis.title=element_text(size = size2),axis.text=element_text(size = size1))
}

sojournTimePlot <- function(model){
  
  par(mfrow=c(1,1))
  d2site <- seq(min(model$data$d2site),max(model$data$d2site),length=250)
  Qmat <- getTrProbs(data.frame(ID=1,time=1:length(d2site),d2site=d2site),nbStates=length(model$stateNames),Time.name="time",beta=model$CIbeta$beta$est,workBounds=model$conditions$workBounds["beta"],formula=model$conditions$formula,betaRef=model$conditions$betaRef,stateNames=model$stateNames,rateMatrix=TRUE,kappa=model$conditions$kappa)
  
  d2mean <- mean(covlist$d2site[cellInd])
  d2sd <- sd(covlist$d2site[cellInd])
  d2site <- d2site * d2sd+d2mean
  nbStates <- length(model$stateNames)
  st <- viterbi(model)
  cols <- stateCols[[nbStates]]
  
  std2site <- model$data$d2site[which(st==1)] * d2sd+d2mean
  counts <- hist(std2site,breaks=d2site,plot=FALSE,right=FALSE)$counts
  ind <- which(d2site<=max(std2site) & d2site>=min(std2site))
  ylim <- range(mapply(function(x) -1/apply(Qmat,3,function(y) y[x,x]),1:nbStates))
  plot(d2site[ind],-1/apply(Qmat,3,function(x) x[1,1])[ind],type="l",col=alpha(1,0),xlim=range(d2site),ylim=ylim,xlab="d2site",ylab=expression(1/q[ii]),main=model$modelName)
  for(i in 1:nbStates){
    message(paste0("State ",i," delta_t summary:"))
    print(summary(model$data$dt[which(st==i)]))
    
    std2site <- model$data$d2site[which(st==i)] * d2sd+d2mean
    counts <- hist(std2site,breaks=d2site,plot=FALSE,right=FALSE)$counts
    ind <- which(d2site<=max(std2site) & d2site>=min(std2site))
    stQmat <- -1/apply(Qmat,3,function(x) x[i,i])[ind]
    message(paste0("State ",i," d2site summary:"))
    print(summary(std2site))
    stdt <- model$data$dt[which(st==i)]
    message("\n")
    segments(d2site[ind[-length(ind)]],stQmat[-length(ind)],d2site[ind[-1]],stQmat[-1],xlim=range(d2site),col=cols[i],lwd=7.5*counts[ind]/max(counts[ind],na.rm=TRUE)+2.0*(counts[ind]>0))
    reg <- loess(stdt~std2site,family="symmetric")
    points(sort(std2site),predict(reg)[order(std2site)],col=cols[i],lty=2,lwd=1.5,type="l")
    #abline(h=median(stdt),lty=i+1,col=cols[i],lwd=1.5)
  }
  points(model$data$d2site * d2sd+d2mean,model$data$dt,col=alpha(cols[st],0.3),pch=20,cex=0.3)
  legend(250,ylim[2]*4/5,legend=model$stateNames,col=cols,lty=1,lwd=2)
  
}

# modification of Rhabit:::simMALA independently by state
simMALA <- function (beta, corr = 0, gamma2 = 1, times, loc0, cov_list = NULL, grad_fun = NULL, 
          silent = FALSE) 
{
  Rhabit:::checkCovGrad(cov_list, grad_fun)
  nb_obs <- length(times)
  J <- length(beta)
  xy <- matrix(NA, nb_obs, 2)
  xy[1, ] <- loc0
  dt <- diff(times)
  covgrad <- Rhabit:::gradLogUD(beta = beta, loc = xy[1, ], cov_list = cov_list, 
                       grad_fun = grad_fun, check = FALSE)
  g <- covgrad %*% beta
  logRSF <- Rhabit:::logRSFinterp(locs = xy[1, ], beta = beta, cov_list = cov_list)
  acc <- 0
  rej <- 0
  lag <- 0
  for (t in 2:nb_obs) {
    if (!silent) 
      cat("\rSimulating Langevin process...", round(100 * 
                                                      t/nb_obs), "%")
    rand_part <- stats::rnorm(2, 0, sqrt(gamma2 * dt[t - 1]))
    
    if(t>2) lag <- (xy[t - 1, ]-xy[t - 2, ]) * dt[t - 1]
    xyprime <- xy[t - 1, ] + lag * corr + 0.5 * g * gamma2 * dt[t - 1] + rand_part

    covgrad <- Rhabit:::gradLogUD(beta = beta, loc = xyprime, cov_list = cov_list, 
                         grad_fun = grad_fun, check = FALSE)
    gprime <- covgrad %*% beta
    logRSFprime <- Rhabit:::logRSFinterp(locs = xyprime, beta = beta, 
                                cov_list = cov_list)
    logProp <- sum(dnorm(xy[t - 1, ] + lag * corr, xyprime + gamma2 * 
                           dt[t - 1] * gprime/2, sqrt(gamma2 * dt[t - 1]), 
                         log = TRUE))
    logPropPrime <- sum(dnorm(xyprime, xy[t - 1, ] + lag * corr + gamma2 * 
                                dt[t - 1] * g/2, sqrt(gamma2 * dt[t - 1]), log = TRUE))
    logAR <- logRSFprime + logProp - logRSF - logPropPrime
    if (log(runif(1)) < logAR) {
      xy[t, ] <- xyprime
      g <- gprime
      logRSF <- logRSFprime
      acc <- acc + 1
    }
    else {
      xy[t, ] <- xy[t - 1, ]
      rej <- rej + 1
    }
  }
  if (!silent) 
    cat("\n")
  main_df <- data.frame(x = xy[, 1], y = xy[, 2], t = times)
  return(list(data = main_df, acc = acc, rej = rej))
}


# modification of Rhabit:::simMALA for state-dependent UDs
simStateMALA <- function (beta, corr = 0, gamma2 = 1, times, states, loc0, cov_list = NULL, grad_fun = NULL, 
                          silent = FALSE) 
{
  lapply(cov_list,function(x) Rhabit:::checkCovGrad(x, grad_fun))
  nb_obs <- length(times)
  J <- length(beta)
  xy <- matrix(NA, nb_obs, 2)
  xy[1, ] <- loc0
  dt <- diff(times)
  
  acc <- 0
  rej <- 0
  lag <- 0
  stateAcc <- stateRej <- numeric(length(gamma2))
  
  for (t in 2:nb_obs) {
    if (!silent) 
      cat("\rSimulating Langevin process...", round(100 * 
                                                      t/nb_obs), "%")
    rand_part <- stats::rnorm(2, 0, sqrt(gamma2[[states[t-1]]] * dt[t - 1]))
    
    covgrad <- Rhabit:::gradLogUD(beta = beta[[states[t-1]]], loc = xy[t-1, ], cov_list = cov_list[[states[t-1]]], 
                                  grad_fun = grad_fun, check = FALSE)
    g <- covgrad %*% beta[[states[t-1]]]
    logRSF <- Rhabit:::logRSFinterp(locs = xy[t-1, ], beta = beta[[states[t-1]]], cov_list = cov_list[[states[t-1]]])
    
    if(t>2) lag <- (xy[t - 1, ]-xy[t - 2, ]) * dt[t - 1]
    xyprime <- xy[t - 1, ] + lag * corr[[states[t-1]]] + 0.5 * g * gamma2[[states[t-1]]] * dt[t - 1] + rand_part
    
    covgrad <- Rhabit:::gradLogUD(beta = beta[[states[t-1]]], loc = xyprime, cov_list = cov_list[[states[t-1]]], 
                                  grad_fun = grad_fun, check = FALSE)
    gprime <- covgrad %*% beta[[states[t-1]]]
    logRSFprime <- Rhabit:::logRSFinterp(locs = xyprime, beta = beta[[states[t-1]]], 
                                         cov_list = cov_list[[states[t-1]]])
    logProp <- sum(dnorm(xy[t - 1, ] + lag * corr[[states[t-1]]], xyprime + gamma2[[states[t-1]]] * 
                           dt[t - 1] * gprime/2, sqrt(gamma2[[states[t-1]]] * dt[t - 1]), 
                         log = TRUE))
    logPropPrime <- sum(dnorm(xyprime, xy[t - 1, ] + lag * corr[[states[t-1]]] + gamma2[[states[t-1]]] * 
                                dt[t - 1] * g/2, sqrt(gamma2[[states[t-1]]] * dt[t - 1]), log = TRUE))
    logAR <- logRSFprime + logProp - logRSF - logPropPrime
    if (log(runif(1)) < logAR) {
      xy[t, ] <- xyprime
      g <- gprime
      logRSF <- logRSFprime
      acc <- acc + 1
      stateAcc[states[t-1]] <- stateAcc[states[t-1]] + 1
    }
    else {
      xy[t, ] <- xy[t - 1, ]
      rej <- rej + 1
      stateRej[states[t-1]] <- stateRej[states[t-1]] + 1
    }
  }
  if (!silent) 
    cat("\n")
  main_df <- data.frame(x = xy[, 1], y = xy[, 2], t = times)
  return(list(data = main_df, acc = acc, rej = rej, stateAcc = stateAcc, stateRej = stateRej))
}

# Simulate MALA from fitted model
malaSim <- function(model,UDinput,niter,globalStates=FALSE,sepStates=FALSE,ssl=TRUE,seed=114){

  set.seed(seed,kind="Mersenne-Twister",normal.kind = "Inversion")

  # parameter and covariate indicators
  parmInd <- UDinput$parIndex
  covInd <- UDinput$covNames
  
  nbAnimals <- length(unique(model$data$ID))
  nbStates <- length(model$stateNames)
  if(nbStates>1){
    if(globalStates) states <- viterbi(model)
    else states <- apply(stateProbs(model),1,which.max)
  } else states <- rep(1,nrow(model$data))
  
  if(nbStates>1) allrates <- matrix(NA, niter, nbAnimals + nbStates*nbAnimals,dimnames = list(NULL,c(paste0("ind",1:nbAnimals,":overall"),paste0("ind",rep(1:nbAnimals,nbStates),":",rep(model$stateNames,each=nbAnimals)))))
  else allrates <- matrix(NA,niter,nbAnimals,dimnames = list(NULL,paste0("ind",1:nbAnimals,":overall")))
  
  if(sepStates){
    acc <- rej <- matrix(0,niter,nbAnimals)
    # loop over states
    for(st in 1:nbStates){
      cat("\n######################")
      cat("\n######################\n")
      cat("## UD:",model$stateNames[st])
      cat("\n######################")
      cat("\n######################\n")
      # Loop over the tracks
      for(zoo in 1:nbAnimals) {
        cat("\n######################\n")
        cat("## Individual",zoo)
        cat("\n######################\n")
        ind <- which(tracks$ID==unique(tracks$ID)[zoo])
        stind <- which(tracks$ID==unique(tracks$ID)[zoo] & st==st)
        sti <- which(st[ind]==st)
        sigma2 <-  (CIreal(model,covs=data.frame(ID=unique(tracks$ID)[zoo]),parms="mu")$mu$est[3,st])^2
        if(nbStates==1 || st < (nbStates+ifelse(ssl,0,1))){
          beta <- model$CIbeta$mu$est[parmInd[[st]]]
          covs <- lapply(covlist0[covInd[[st]]],rasterToRhabit)
          corr <- rep(0,2)
          if(any(grepl(paste0("mean.x_",st,":crw(mu.x_tm1)"),colnames(model$CIbeta$mu$est),fixed=TRUE))){
            corr[1] <- model$CIbeta$mu$est[which(grepl(paste0("mean.x_",st,":crw(mu.x_tm1)"),colnames(model$CIbeta$mu$est),fixed=TRUE))]
          }
          if(any(grepl(paste0("mean.y_",st,":crw(mu.y_tm1)"),colnames(model$CIbeta$mu$est),fixed=TRUE))){
            corr[2] <- model$CIbeta$mu$est[which(grepl(paste0("mean.y_",st,":crw(mu.y_tm1)"),colnames(model$CIbeta$mu$est),fixed=TRUE))]
          }
        } else {
          beta <- rep(0,length(covInd[[1]]))
          covs <- lapply(covlist0[covInd[[1]]],rasterToRhabit)
          corr <- rep(0,2)
        }
        
        times <- as.numeric(tracks$date_time[ind])
        times <- (times-min(times))/3600
        
        # get state segments
        seqs <- R.utils::seqToIntervals(sti[-length(sti)])
        seqs[,2] <- seqs[,2]+1
        
        # Loop over the MALA simulations
        iter <- 1
        while(iter<=niter) {
          cat("Iteration",iter,"\n")
          tryCatch({
            # MALA simulation, based on estimated parameters and time grid of observation
            
            sim <- list()

            for(sq in 1:nrow(seqs)){
              sim[[sq]] <- simMALA(beta = beta, corr = corr, gamma2 = sigma2, 
                                            times = times[seq(seqs[sq,1],seqs[sq,2])], loc0 = c(tracks[seqs[sq,1],"mu.x"],tracks[seqs[sq,1],"mu.y"]), 
                                            cov_list = covs, silent=TRUE)
            }
            
            # Acceptance rate
            acc[iter,zoo] <- acc[iter,zoo] + sum(unlist(lapply(sim,function(x) x$acc)))
            rej[iter,zoo] <- rej[iter,zoo] + sum(unlist(lapply(sim,function(x) x$rej)))
            if(nbStates>1) {
              rate <- sum(unlist(lapply(sim,function(x) x$acc)))/(sum(unlist(lapply(sim,function(x) x$acc)))+sum(unlist(lapply(sim,function(x) x$rej))))
              allrates[iter,nbAnimals+seq(zoo,nbAnimals*nbStates,nbAnimals)] <- rate
            }
            iter <- iter + 1
          }, error = function(e) {
            # If the simulated track leaves the study region, start over
            cat("Error in the simulation -- starting over.\n")
          })
        }
      }
    }
    allrates[,1:nbAnimals] <- acc / (acc + rej)
  } else {
    
    # Loop over the tracks
    for(zoo in 1:nbAnimals) {
      cat("\n######################\n")
      cat("## Individual",zoo)
      cat("\n######################\n")
      ind <- which(tracks$ID==unique(tracks$ID)[zoo])
      beta <- sigma2 <- covs <- corr <- list()
      for(st in 1:nbStates){
        sigma2[[st]] <-  (CIreal(model,covs=data.frame(ID=unique(tracks$ID)[zoo]),parms="mu")$mu$est[3,st])^2
        if(nbStates==1 ||  st < (nbStates+ifelse(ssl,0,1))){
          beta[[st]] <- model$CIbeta$mu$est[parmInd[[st]]]
          covs[[st]] <- lapply(covlist0[covInd[[st]]],rasterToRhabit)
          corr[[st]] <- rep(0,2)
          if(any(grepl(paste0("mean.x_",st,":crw(mu.x_tm1)"),colnames(model$CIbeta$mu$est),fixed=TRUE))){
            corr[[st]][1] <- model$CIbeta$mu$est[which(grepl(paste0("mean.x_",st,":crw(mu.x_tm1)"),colnames(model$CIbeta$mu$est),fixed=TRUE))]
          }
          if(any(grepl(paste0("mean.y_",st,":crw(mu.y_tm1)"),colnames(model$CIbeta$mu$est),fixed=TRUE))){
            corr[[st]][2] <- model$CIbeta$mu$est[which(grepl(paste0("mean.y_",st,":crw(mu.y_tm1)"),colnames(model$CIbeta$mu$est),fixed=TRUE))]
          }
        } else {
          beta[[nbStates]] <- rep(0,length(covInd[[1]]))
          covs[[nbStates]] <- lapply(covlist0[covInd[[1]]],rasterToRhabit)
          corr[[nbStates]] <- rep(0,2)
        }
      }
      
      times <- as.numeric(tracks$date_time[ind])
      times <- (times-min(times))/3600
      st <- states[ind]
      
      # Loop over the MALA simulations
      iter <- 1
      while(iter<=niter) {
        cat("Iteration",iter,"\n")
        tryCatch({
          # MALA simulation, based on estimated parameters and time grid of observation
          
          sim  <- simStateMALA(beta = beta, corr = corr, gamma2 = sigma2, 
                                 times = times, states = st, loc0 = c(tracks[ind[1],"mu.x"],tracks[ind[1],"mu.y"]), 
                                 cov_list = covs)
          
          # Acceptance rate
          rate <- sim$acc/(sim$acc+sim$rej)
          stateRate <- sim$stateAcc / (sim$stateAcc + sim$stateRej)
          allrates[iter, zoo] <- rate
          if(nbStates>1) allrates[iter,nbAnimals+seq(zoo,nbAnimals*nbStates,nbAnimals)] <- stateRate
          
          iter <- iter + 1
        }, error = function(e) {
          # If the simulated track leaves the study region, start over
          cat("Error in the simulation -- starting over.\n")
        })
      }
    }
  }
  return(allrates)
}


statePlot <- function(model,UD){
  
  stateNames <- model$stateNames
  nbStates <- length(stateNames)
  st <- viterbi(model)
  
  land <- (covlistCrop$depth*sd(covlist$depth[cellInd])+mean(covlist$depth[cellInd])-attr(covlist$depth,"absmin"))>=0
  
  par(mfrow=c(nbStates-2,2), mar=c(0,0,3,3)+0.1)
  raster::plot(land, col=gray.colors(10, start = 1, end = 0.4, gamma = 2.2, alpha = NULL), axes=FALSE, legend=FALSE,box=FALSE,bty="n")
  points(model$data$mu.x_tm1,model$data$mu.y_tm1,col=alpha(stateCols[[nbStates]][st],0.3),pch=20,cex=2)
  points(model$data$mu.x_tm1[which(st==nbStates)],model$data$mu.y_tm1[which(st==nbStates)],col=alpha(stateCols[[nbStates]][nbStates],0.4),pch=20,cex=2)
  segments(model$data$mu.x_tm1[-nrow(model$data)],model$data$mu.y_tm1[-nrow(model$data)],model$data$mu.x_tm1[-1],model$data$mu.y_tm1[-1],col=alpha(stateCols[[nbStates]][st],0.4),cex=2)
  legend(-1600,900,stateNames,col=stateCols[[nbStates]],pch=20,cex=2)
  for(i in 1:(nbStates-1)){
    raster::plot(UD[[stateNames[i]]], col=viridis_pal()(100), axes = FALSE, main=stateNames[[i]], legend.args = list(text = parse(text=(paste0("log(pi[",i,"])"))),cex=1.2), legend.mar=8, box=FALSE, cex.main=2.5,legend.width=3,axis.args=list(cex.axis=1.4))
    points(model$data$mu.x_tm1[which(st==i)],model$data$mu.y_tm1[which(st==i)],col=alpha(stateCols[[nbStates]][i],0.4),pch=20,cex=2)
  }
}
