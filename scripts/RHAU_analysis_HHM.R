library(here)
library(ggplot2)
library(momentuHMM)
library(sf)
library(dplyr)
#########################################################
## RHAU DATA
# Load RHAU data
RHAU_data <- read.csv(here("data", "RHAU_data", "RHAU_data_UTM.csv"))
RHAU_data$time <- as.POSIXct(RHAU_data$time,tz="UTC") #convert times to POSIX 

RHAU_data <- RHAU_data %>%
  group_by(ID) %>%
  arrange(time, .by_group = TRUE)

# Multiple imputation to address temporal irregularity
# fit crawl model
crwOut <- crawlWrap(obsData = RHAU_data, timeStep = "15 mins",
                    theta=list(c(1.993,1.272),
                               c(1.705,2.053),
                               c(1.795,1.745),
                               c(1.449,1.893),
                               c(1.046,2.291),
                               c(1.509,2.251),
                               c(1.237,1.646)), #theta values obtained from iterating through different values and monitoring output; see below
                    fixPar=c(NA,NA))

crwOutFits <- crwOut$crwFits #observe parameter estimates and AIC/LL
print(crwOutFits) #iterate through theta values above (crawlWrap) using sigma and beta

RHAU_data_crw <- crwOut$crwPredict[,c(3,6,9:17)] #isolate interpolated data

RHAU_data_crw <- RHAU_data_crw %>%
  group_by(ID) %>%
  arrange(time, .by_group = TRUE)

#Create time of day covariate
RHAU_data_crw$tod <- as.numeric(format(RHAU_data_crw$time, "%H")) + as.numeric(format(RHAU_data_crw$time, "%M"))/60 #decimal hours

#Finish data
RHAU_data <- RHAU_data_crw

#########################################################
## MODEL SETUP

# Prepare data for momentuHMM and calculate gradients
tracks <- prepData(data = RHAU_data,
                   type = "UTM",
                   coordNames = c("mu.x", "mu.y"))

RHAU_data$hour <- as.integer(strftime(RHAU_data$time, format = "%H", tz="GMT")) #add cosinor covariate base on hour of day


tracks <- tracks %>%
  filter(if_all(everything(), ~ !is.na(.)))


# Get an idea of dist parameter specifications by plotting step length and turn angle
# Step length
hist(tracks$step, breaks=50)
summary(tracks$step) #look at mean and median
sd(tracks$step, na.rm = TRUE)
# Turn angle
hist(tracks$angle, breaks = 25)
summary(tracks$angle)
sd(tracks$angle, na.rm = TRUE)



stateNames <- c("resting", "exploring", "foraging")

stepMu0 <- c(10, 1000, 3000) #refer to stateNames. low resting movement, high exploring movement, intermediate foraging movement
stepSig0 <- c(0.5, 500, 200)
stepPar0 <- c(stepMu0, stepSig0)

anglePar0 <- c(10, 8, 0.4) #refer to stateNames. initial concentration parameters; little variation in angle, intermediate, high

# MODEL 1
# Fit first model without TOD covariate
RHAU_HMM_1 <- fitHMM(data = tracks, nbStates = 3, 
                     dist = list(step="gamma", angle="vm"), 
                     Par0 = list(step=stepPar0, angle=anglePar0),
                     estAngleMean = list(angle=FALSE),
                     stateNames = stateNames)

# MODEL 2
# Formula for transition probabilities; add TOD
formula <- ~cosinor(tod, period = 24)
# Get parameter estimates to initialize next model
Par0_2 <- getPar0(model = RHAU_HMM_1,
                  formula = formula)

# Use parameter estimates and coefficient from previous model; include TOD covariate
RHAU_HMM_2 <- fitHMM(data = tracks, nbStates = 3, 
                     dist = list(step="gamma", angle="vm"), 
                     Par0 = list(step=Par0_2$Par$step, angle=Par0_2$Par$angle),
                     estAngleMean = list(angle=FALSE),
                     stateNames = stateNames,
                     formula = formula)


# MODEL VALIDATION
states <- viterbi(RHAU_HMM_1) #decode most likely state sequence; change out HMM models
table(states)/nrow(tracks) #derive percentage of time spent in each state

# Compare models
AIC(RHAU_HMM_1, RHAU_HMM_2)

# Psuedo-residuals for steps and angles
pr <- pseudoRes(RHAU_HMM_1)

# Plot ACF of step pseudo-residuals
acf(pr$stepRes[!is.na(pr$stepRes)],lag.max = 300)

# Save the output of best model!
save(RHAU_HMM_1, file = here("results", "RHAU_HMM_1.RData"))
################################################################################
