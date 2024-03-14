modelName <- nbStates <- DM <- formula <- workBounds <- stateNames <- inputUD <- fixPar <- kappa <- list()

##### 1 state models #####
# add additional habitat covariates
modelName[[1]] <- "S=1 sigma(.) beta(depth*slope+d2site^2)"
nbStates[[1]] <- 1
DM[[1]] <- list(mu = list(mean.x=~0+mu.x_tm1+langevin(depth.x)+langevin(slope.x)+langevin(depthslope.x)+langevin(d2site.x)+langevin(d2site2.x),
                          mean.y=~0+mu.y_tm1+langevin(depth.y)+langevin(slope.y)+langevin(depthslope.y)+langevin(d2site.y)+langevin(d2site2.y),
                          sd.x=~1,
                          sd.y=~1,
                          corr.xy=~1),
               dry = list(prob=~1))
formula[[1]] <- ~1
workBounds[[1]] <- getWorkBounds(DM[[1]],nbStates[[1]],tracks) 
fixPar[[1]] <- getFixPar(DM[[1]],nbStates[[1]],tracks)
kappa[[1]] <- Inf
stateNames[[1]] <- "outbound"
inputUD[[1]] <- list(nbUD=1,parIndex=getLangIndex(fixPar[[1]],nbStates[[1]]),covNames=list(c("depth","slope","depthslope","d2site","d2site2")),UDnames="outbound",UDstates=list(1),sign=list(1))


# add ID effects on speed parameter
modelName[[2]] <- "S=1 sigma(ID) beta(depth*slope+d2site^2)"
nbStates[[2]] <- 1
DM[[2]] <- list(mu = list(mean.x=~0+mu.x_tm1+langevin(depth.x)+langevin(slope.x)+langevin(depthslope.x)+langevin(d2site.x)+langevin(d2site2.x),
                          mean.y=~0+mu.y_tm1+langevin(depth.y)+langevin(slope.y)+langevin(depthslope.y)+langevin(d2site.y)+langevin(d2site2.y),
                          sd.x=~0+ID,
                          sd.y=~0+ID,
                          corr.xy=~1),
                dry = list(prob=~1))
formula[[2]] <- ~1
workBounds[[2]] <- getWorkBounds(DM[[2]],nbStates[[2]],tracks) 
fixPar[[2]] <- getFixPar(DM[[2]],nbStates[[2]],tracks)
kappa[[2]] <- Inf
stateNames[[2]] <- "outbound"
inputUD[[2]] <- list(nbUD=1,parIndex=getLangIndex(fixPar[[2]],nbStates[[2]]),covNames=list(c("depth","slope","depthslope","d2site","d2site2")),UDnames="outbound",UDstates=list(1),sign=list(1))

##### 2 state models #####

modelName[[3]] <- "S=2 sigma(.) beta(depth*slope+d2site^2) gamma(d2site)"
nbStates[[3]] <- 2
DM[[3]] <- list(mu = list(mean.x=~0+mu.x_tm1+state1(langevin(depth.x)+langevin(slope.x)+langevin(depthslope.x)+langevin(d2site.x)+langevin(d2site2.x)),
                          mean.y=~0+mu.y_tm1+state1(langevin(depth.y)+langevin(slope.y)+langevin(depthslope.y)+langevin(d2site.y)+langevin(d2site2.y)),
                          sd.x=~state2(intercept),
                          sd.y=~state2(intercept),
                          corr.xy=~1),
                dry = list(prob=~1))
formula[[3]] <- ~state1(d2site)
workBounds[[3]] <- getWorkBounds(DM[[3]],nbStates[[3]],tracks) 
fixPar[[3]] <- getFixPar(DM[[3]],nbStates[[3]],tracks)
kappa[[3]] <- 4
stateNames[[3]] <- c("outbound","haulout")
inputUD[[3]] <- list(nbUD=1,parIndex=getLangIndex(fixPar[[3]],nbStates[[3]]),covNames=list(c("depth","slope","depthslope","d2site","d2site2")),UDnames="outbound",UDstates=list(1),sign=list(1))

# add ID effects on speed parameter for "at sea" state
modelName[[4]] <- "S=2 sigma(ID) beta(depth*slope+d2site^2) gamma(d2site)"
nbStates[[4]] <- 2
DM[[4]] <- DM[[3]]
DM[[4]]$mu$sd.x <- DM[[4]]$mu$sd.y <- ~state1(0+ID)+state2(intercept)
formula[[4]] <- ~state1(d2site)
workBounds[[4]] <- getWorkBounds(DM[[4]],nbStates[[4]],tracks) 
fixPar[[4]] <- getFixPar(DM[[4]],nbStates[[4]],tracks)#list(mu=c(NA,1,2,3,4,5,NA,NA,1,2,3,4,5,NA,6,7,8,6,9,6,7,8,6,9,NA,NA),
                    #delta=rep(NA,nbStates[[4]]))
kappa[[4]] <- 4
stateNames[[4]] <- c("outbound","haulout")
inputUD[[4]] <- list(nbUD=1,parIndex=getLangIndex(fixPar[[4]],nbStates[[4]]),covNames=list(c("depth","slope","depthslope","d2site","d2site2")),UDnames="outbound",UDstates=list(1),sign=list(1))

##### 3 state models #####

modelName[[5]] <- "S=3 sigma(ID,1=2) beta(depth*slope+d2site^2) p(1=2) gamma(d2site)"
nbStates[[5]] <- 3
DM[[5]] <- list(mu = list(mean.x=~0+mu.x_tm1+state1(langevin(depth.x)+langevin(slope.x)+langevin(depthslope.x)+langevin(d2site.x)+langevin(d2site2.x))
                                          +state2(langevin(depth.x)+langevin(slope.x)+langevin(depthslope.x)+langevin(d2site.x)+langevin(d2site2.x)),
                          mean.y=~0+mu.y_tm1+state1(langevin(depth.y)+langevin(slope.y)+langevin(depthslope.y)+langevin(d2site.y)+langevin(d2site2.y))
                                          +state2(langevin(depth.y)+langevin(slope.y)+langevin(depthslope.y)+langevin(d2site.y)+langevin(d2site2.y)),
                          sd.x=~state1(0+ID)+state2(0+ID)+state3(intercept),
                          sd.y=~state1(0+ID)+state2(0+ID)+state3(intercept),
                          corr.xy=~1),
                dry = list(prob=~1))
formula[[5]] <- ~state1(d2site)+state2(d2site)
workBounds[[5]] <- getWorkBounds(DM[[5]],nbStates[[5]],tracks) 

fixPar[[5]] <- getFixPar(DM[[5]],nbStates[[5]],tracks,sdcons=list(c(1,2)),drycons=list(c(1,2))) #list(mu=c(NA,1,2,3,4,5,NA,6,7,8,9,10,NA,NA,1,2,3,4,5,NA,6,7,8,9,10,NA,11,12,13,11,12,13,11,14,11,12,13,11,12,13,11,14,NA,NA,NA),
                    #delta=rep(NA,nbStates[[5]]))
kappa[[5]] <- 4
stateNames[[5]] <- c("outbound","inbound","haulout")
inputUD[[5]] <- list(nbUD=2,parIndex=getLangIndex(fixPar[[5]],nbStates[[5]]),covNames=list(c("depth","slope","depthslope","d2site","d2site2"),c("depth","slope","depthslope","d2site","d2site2")),UDnames=c("outbound","inbound"),UDstates=list(1,2),sign=list(1,1))

# allow "inbound" state to have own speed and p
modelName[[6]] <- "S=3 sigma(ID) beta(depth*slope+d2site^2) gamma(d2site)"
nbStates[[6]] <- 3
DM[[6]] <- DM[[5]]
formula[[6]] <- ~state1(d2site)+state2(d2site)
workBounds[[6]] <- getWorkBounds(DM[[6]],nbStates[[6]],tracks) 

fixPar[[6]] <- getFixPar(DM[[6]],nbStates[[6]],tracks)
fixPar[[6]]$beta <- c(1:9,NA,10,NA)
kappa[[6]] <- 4
stateNames[[6]] <- c("outbound","inbound","haulout")
inputUD[[6]] <- list(nbUD=2,parIndex=getLangIndex(fixPar[[6]],nbStates[[6]]),covNames=list(c("depth","slope","depthslope","d2site","d2site2"),c("depth","slope","depthslope","d2site","d2site2")),UDnames=c("outbound","inbound"),UDstates=list(1,2),sign=list(1,1))

##### 4 state models #####

# "outbound", "foraging1", and "inbound" have same speed; "outbound" and "inbound" have opposite potential surfaces
modelName[[7]] <- "S=4 sigma(ID,1=2=3) beta(depth*slope+d2site^2,1=-3) p(1=2=3) gamma(d2site)"
nbStates[[7]] <- 4
DM[[7]] <- list(mu = list(mean.x=~0+mu.x_tm1+state1(langevin(depth.x)+langevin(slope.x)+langevin(depthslope.x)+langevin(d2site.x)+langevin(d2site2.x))
                          +state2(langevin(depth.x)+langevin(slope.x)+langevin(depthslope.x)+langevin(d2site.x))
                          +state3(langevin(-depth.x)+langevin(-slope.x)+langevin(-depthslope.x)+langevin(-d2site.x)+langevin(-d2site2.x)),
                          mean.y=~0+mu.y_tm1+state1(langevin(depth.y)+langevin(slope.y)+langevin(depthslope.y)+langevin(d2site.y)+langevin(d2site2.y))
                          +state2(langevin(depth.y)+langevin(slope.y)+langevin(depthslope.y)+langevin(d2site.y))
                          +state3(langevin(-depth.y)+langevin(-slope.y)+langevin(-depthslope.y)+langevin(-d2site.y)+langevin(-d2site2.y)),
                          sd.x=~state1(0+ID)+state2(0+ID)+state3(0+ID)+state4(intercept),
                          sd.y=~state1(0+ID)+state2(0+ID)+state3(0+ID)+state4(intercept),
                          corr.xy=~1),
                dry = list(prob=~1))
formula[[7]] <- ~state1(d2site)+state2(d2site)+state3(d2site)
workBounds[[7]] <- getWorkBounds(DM[[7]],nbStates[[7]],tracks) 
fixPar[[7]] <- getFixPar(DM[[7]],nbStates[[7]],tracks,meancons=list(c(1,3)),sdcons = list(c(1,2,3)),drycons=list(c(1,2,3)))
kappa[[7]] <- 4
stateNames[[7]] <- c("outbound","foraging","inbound","haulout")
inputUD[[7]] <- list(nbUD=3,parIndex=getLangIndex(fixPar[[7]],nbStates[[7]])[c(1,2,1)],covNames=list(c("depth","slope","depthslope","d2site","d2site2"),c("depth","slope","depthslope","d2site"),c("depth","slope","depthslope","d2site","d2site2")),UDnames=c("outbound","foraging","inbound"),UDstates=list(1,2,3),sign=list(1,1,-1))

# allow "foraging1" state to have own speed
modelName[[8]] <- "S=4 sigma(ID,1=3) beta(depth*slope+d2site^2,1=-3) p(1=3) gamma(d2site)"
nbStates[[8]] <- 4
DM[[8]] <- DM[[7]]
formula[[8]] <- ~state1(d2site)+state2(d2site)+state3(d2site)
workBounds[[8]] <- getWorkBounds(DM[[8]],nbStates[[8]],tracks) 
fixPar[[8]] <- getFixPar(DM[[8]],nbStates[[8]],tracks,meancons=list(c(1,3)),sdcons = list(c(1,3)),drycons=list(c(1,3)))
kappa[[8]] <- 4
stateNames[[8]] <- c("outbound","foraging","inbound","haulout")
inputUD[[8]] <- list(nbUD=3,parIndex=getLangIndex(fixPar[[8]],nbStates[[8]])[c(1,2,1)],covNames=list(c("depth","slope","depthslope","d2site","d2site2"),c("depth","slope","depthslope","d2site"),c("depth","slope","depthslope","d2site","d2site2")),UDnames=c("outbound","foraging","inbound"),UDstates=list(1,2,3),sign=list(1,1,-1))

# "outbound", "foraging1", and "inbound" have same speed but each state has own potential surface
modelName[[9]] <- "S=4 sigma(ID,1=2=3) beta(depth*slope+d2site^2) p(1=2=3) gamma(d2site)"
nbStates[[9]] <- 4
DM[[9]] <- list(mu = list(mean.x=~0+mu.x_tm1+state1(langevin(depth.x)+langevin(slope.x)+langevin(depthslope.x)+langevin(d2site.x)+langevin(d2site2.x))
                          +state2(langevin(depth.x)+langevin(slope.x)+langevin(depthslope.x)+langevin(d2site.x))
                          +state3(langevin(depth.x)+langevin(slope.x)+langevin(depthslope.x)+langevin(d2site.x)+langevin(d2site2.x)),
                          mean.y=~0+mu.y_tm1+state1(langevin(depth.y)+langevin(slope.y)+langevin(depthslope.y)+langevin(d2site.y)+langevin(d2site2.y))
                          +state2(langevin(depth.y)+langevin(slope.y)+langevin(depthslope.y)+langevin(d2site.y))
                          +state3(langevin(depth.y)+langevin(slope.y)+langevin(depthslope.y)+langevin(d2site.y)+langevin(d2site2.y)),
                          sd.x=~state1(0+ID)+state2(0+ID)+state3(0+ID)+state4(intercept),
                          sd.y=~state1(0+ID)+state2(0+ID)+state3(0+ID)+state4(intercept),
                          corr.xy=~1),
                dry = list(prob=~1))
formula[[9]] <- ~state1(d2site)+state2(d2site)+state3(d2site)
workBounds[[9]] <- getWorkBounds(DM[[9]],nbStates[[9]],tracks) 
fixPar[[9]] <- getFixPar(DM[[9]],nbStates[[9]],tracks,sdcons = list(c(1,2,3)),drycons=list(c(1,2,3)))
kappa[[9]] <- 4
stateNames[[9]] <- c("outbound","foraging","inbound","haulout")
inputUD[[9]] <- list(nbUD=3,parIndex=getLangIndex(fixPar[[9]],nbStates[[9]]),covNames=list(c("depth","slope","depthslope","d2site","d2site2"),c("depth","slope","depthslope","d2site"),c("depth","slope","depthslope","d2site","d2site2")),UDnames=c("outbound","foraging","inbound"),UDstates=list(1,2,3),sign=list(1,1,1))

# allow "foraging1" state to have own speed, "outbound" and "inbound" have own potential surfaces
modelName[[10]] <- "S=4 sigma(ID,1=3) beta(depth*slope+d2site^2) p(1=3) gamma(d2site)"
nbStates[[10]] <- 4
DM[[10]] <- DM[[9]]
formula[[10]] <- ~state1(d2site)+state2(d2site)+state3(d2site)
workBounds[[10]] <- getWorkBounds(DM[[10]],nbStates[[10]],tracks) 
fixPar[[10]] <- getFixPar(DM[[10]],nbStates[[10]],tracks,sdcons = list(c(1,3)),drycons=list(c(1,3)))
kappa[[10]] <- 4
stateNames[[10]] <- c("outbound","foraging","inbound","haulout")
inputUD[[10]] <- list(nbUD=3,parIndex=getLangIndex(fixPar[[10]],nbStates[[10]]),covNames=list(c("depth","slope","depthslope","d2site","d2site2"),c("depth","slope","depthslope","d2site"),c("depth","slope","depthslope","d2site","d2site2")),UDnames=c("outbound","foraging","inbound"),UDstates=list(1,2,3),sign=list(1,1,1))

# "outbound", "foraging1", and "inbound" have own speed, dry, and potential surface
modelName[[11]] <- "S=4 sigma(ID) beta(depth*slope+d2site^2) gamma(d2site)"
nbStates[[11]] <- 4
DM[[11]] <- DM[[9]]
formula[[11]] <- ~state1(d2site)+state2(d2site)+state3(d2site)
workBounds[[11]] <- getWorkBounds(DM[[11]],nbStates[[11]],tracks) 
fixPar[[11]] <- getFixPar(DM[[11]],nbStates[[11]],tracks)
fixPar[[11]]$beta <- c(1:19,NA,20,NA,21,NA)
kappa[[11]] <- 4
stateNames[[11]] <- c("outbound","foraging","inbound","haulout")
inputUD[[11]] <- list(nbUD=3,parIndex=getLangIndex(fixPar[[11]],nbStates[[11]]),covNames=list(c("depth","slope","depthslope","d2site","d2site2"),c("depth","slope","depthslope","d2site"),c("depth","slope","depthslope","d2site","d2site2")),UDnames=c("outbound","foraging","inbound"),UDstates=list(1,2,3),sign=list(1,1,1))

##### 5 state models #####

# 5 state model where two "foraging" states share same potential surface
modelName[[12]] <- "S=5 sigma(ID,1=3=4) beta(depth*slope+d2site^2,2=3) p(1=3=4) gamma(d2site)"
nbStates[[12]] <- 5
DM[[12]] <- list(mu = list(mean.x=~0+mu.x_tm1+state1(langevin(depth.x)+langevin(slope.x)+langevin(depthslope.x)+langevin(d2site.x)+langevin(d2site2.x))
                           +state2(langevin(depth.x)+langevin(slope.x)+langevin(depthslope.x)+langevin(d2site.x))
                           +state3(langevin(depth.x)+langevin(slope.x)+langevin(depthslope.x)+langevin(d2site.x))
                           +state4(langevin(depth.x)+langevin(slope.x)+langevin(depthslope.x)+langevin(d2site.x)+langevin(d2site2.x)),
                           mean.y=~0+mu.y_tm1+state1(langevin(depth.y)+langevin(slope.y)+langevin(depthslope.y)+langevin(d2site.y)+langevin(d2site2.y))
                           +state2(langevin(depth.y)+langevin(slope.y)+langevin(depthslope.y)+langevin(d2site.y))
                           +state3(langevin(depth.y)+langevin(slope.y)+langevin(depthslope.y)+langevin(d2site.y))
                           +state4(langevin(depth.y)+langevin(slope.y)+langevin(depthslope.y)+langevin(d2site.y)+langevin(d2site2.y)),
                           sd.x=~state1(0+ID)+state2(0+ID)+state3(0+ID)+state4(0+ID)+state5(intercept),
                           sd.y=~state1(0+ID)+state2(0+ID)+state3(0+ID)+state4(0+ID)+state5(intercept),
                           corr.xy=~1),
                 dry = list(prob=~1))
formula[[12]] <- ~state1(d2site)+state2(d2site)+state3(d2site)+state4(d2site)
workBounds[[12]] <- getWorkBounds(DM[[12]],nbStates[[12]],tracks) 
fixPar[[12]] <- getFixPar(DM[[12]],nbStates[[12]],tracks,meancons=list(c(2,3)),sdcons=list(c(1,3,4)),drycons=list(c(1,3,4)))
kappa[[12]] <- 4
stateNames[[12]] <- c("outbound","foraging1","foraging2","inbound","haulout")
inputUD[[12]] <- list(nbUD=4,parIndex=getLangIndex(fixPar[[12]],nbStates[[12]])[c(1,2,2,4)],covNames=list(c("depth","slope","depthslope","d2site","d2site2"),c("depth","slope","depthslope","d2site"),c("depth","slope","depthslope","d2site"),c("depth","slope","depthslope","d2site","d2site2")),UDnames=c("outbound","foraging1","foraging2","inbound"),UDstates=list(1,2,3,4),sign=c(1,1,1,1))

# two "foraging" states have own speeds
modelName[[13]] <- "S=5 sigma(ID,1=4) beta(depth*slope+d2site^2,2=3) p(1=4) gamma(d2site)"
nbStates[[13]] <- 5
DM[[13]] <- DM[[12]] 
formula[[13]] <- ~state1(d2site)+state2(d2site)+state3(d2site)+state4(d2site)
workBounds[[13]] <- getWorkBounds(DM[[13]],nbStates[[13]],tracks) 
fixPar[[13]] <- getFixPar(DM[[13]],nbStates[[13]],tracks,meancons=list(c(2,3)),sdcons=list(c(1,4)),drycons=list(c(1,4)))
kappa[[13]] <- 4
stateNames[[13]] <- c("outbound","foraging1","foraging2","inbound","haulout")
inputUD[[13]] <- list(nbUD=4,parIndex=getLangIndex(fixPar[[13]],nbStates[[13]])[c(1,2,2,4)],covNames=list(c("depth","slope","depthslope","d2site","d2site2"),c("depth","slope","depthslope","d2site"),c("depth","slope","depthslope","d2site"),c("depth","slope","depthslope","d2site","d2site2")),UDnames=c("outbound","foraging1","foraging2","inbound"),UDstates=list(1,2,3,4),sign=c(1,1,1,1))

# two "foraging" states have own potential surface
modelName[[14]] <- "S=5 sigma(ID,1=3=4) beta(depth*slope+d2site^2) p(1=3=4) gamma(d2site)"
nbStates[[14]] <- 5
DM[[14]] <- DM[[12]]
formula[[14]] <- ~state1(d2site)+state2(d2site)+state3(d2site)+state4(d2site)
workBounds[[14]] <- getWorkBounds(DM[[14]],nbStates[[14]],tracks) 
fixPar[[14]] <- getFixPar(DM[[14]],nbStates[[14]],tracks,sdcons=list(c(1,3,4)),drycons=list(c(1,3,4)))
kappa[[14]] <- 4
stateNames[[14]] <- c("outbound","foraging1","foraging2","inbound","haulout")
inputUD[[14]] <- list(nbUD=4,parIndex=getLangIndex(fixPar[[14]],nbStates[[14]]),covNames=list(c("depth","slope","depthslope","d2site","d2site2"),c("depth","slope","depthslope","d2site"),c("depth","slope","depthslope","d2site"),c("depth","slope","depthslope","d2site","d2site2")),UDnames=c("outbound","foraging1","foraging2","inbound"),UDstates=list(1,2,3,4),sign=c(1,1,1,1))

# kitchen sink model where all 5 states have own speed, potential surface, and dry parameters
modelName[[15]] <- "S=5 sigma(ID) beta(depth*slope+d2site^2) gamma(d2site)"
nbStates[[15]] <- 5
DM[[15]] <- DM[[12]]
formula[[15]] <- ~state1(d2site)+state2(d2site)+state3(d2site)+state4(d2site)
workBounds[[15]] <- getWorkBounds(DM[[15]],nbStates[[15]],tracks) 
fixPar[[15]] <- getFixPar(DM[[15]],nbStates[[15]],tracks)
kappa[[15]] <- 4
stateNames[[15]] <- c("outbound","foraging1","foraging2","inbound","haulout")
inputUD[[15]] <- list(nbUD=4,parIndex=getLangIndex(fixPar[[15]],nbStates[[15]]),covNames=list(c("depth","slope","depthslope","d2site","d2site2"),c("depth","slope","depthslope","d2site"),c("depth","slope","depthslope","d2site"),c("depth","slope","depthslope","d2site","d2site2")),UDnames=c("outbound","foraging1","foraging2","inbound"),UDstates=list(1,2,3,4),sign=c(1,1,1,1))

prior <- mapply(function(x) getPrior(fixPar[[x]],sd=10),1:length(DM))
