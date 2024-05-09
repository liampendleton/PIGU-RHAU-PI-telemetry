library(here)
library(momentuHMM)
library(sf)
library(dplyr)
#########################################################
## PIGU DATA
# Load PIGU data
PIGU_data <- read.csv(here("data", "PIGU_data", "PIGU_data_UTM.csv"))
PIGU_data$time <- as.POSIXct(PIGU_data$time,tz="UTC") #convert times to POSIX 

PIGU_data <- PIGU_data %>%
  group_by(ID) %>%
  arrange(time, .by_group = TRUE)

PIGU_data <- PIGU_data %>% filter(y > 5327000) #removing erroneous points

# Multiple imputation to address temporal irregularity
# fit crawl model
crwOut <- crawlWrap(obsData = PIGU_data, timeStep = "15 mins",
                    theta=list(c(7.566,2.715),
                               c(7.552,2.890),
                               c(8.316,14.981),
                               c(6.736,13.340),
                               c(7.157,2.991),
                               c(7.175,2.647)), #theta values obtained from iterating through different values and monitoring output; see below
                    fixPar=c(NA,NA))

crwOutFits <- crwOut$crwFits #observe parameter estimates and AIC/LL
print(crwOutFits) #iterate through theta values above (crawlWrap) using sigma and beta

PIGU_data_crw <- crwOut$crwPredict[,c(3,6,9:17)] #isolate interpolated data

PIGU_data_crw <- PIGU_data_crw %>%
  group_by(ID) %>%
  arrange(time, .by_group = TRUE)

#Create time of day covariate
PIGU_data_crw$tod <- as.numeric(format(PIGU_data_crw$time, "%H")) + as.numeric(format(PIGU_data_crw$time, "%M"))/60 #decimal hours

#Finish data
PIGU_data <- PIGU_data_crw

#########################################################
## MODEL SETUP

# Prepare data for momentuHMM and calculate gradients
tracks <- prepData(data = PIGU_data,
                   type = "UTM",
                   coordNames = c("mu.x", "mu.y"))

tracks$hour <- as.integer(strftime(tracks$time, format = "%H", tz="GMT")) #add cosinor covariate base on hour of day

# Remove NAs
tracks <- tracks %>%
  filter(if_all(everything(), ~ !is.na(.)))

# Get an idea of dist parameter specifications by plotting step length and turn angle
# Step length
hist(tracks$step, breaks = 50, xlim=c(0,4000))
summary(tracks$step) #look at mean and median
sd(tracks$step, na.rm = TRUE)
# Turn angle
hist(tracks$angle, breaks = 100)
summary(tracks$angle)
sd(tracks$angle, na.rm = TRUE)

stateNames <- c("resting", "exploring", "foraging")

stepMu0 <- c(5, 500, 100) #refer to stateNames. low resting movement, high exploring movement, intermediate foraging movement
stepSig0 <- c(2, 100, 50) #deviations
stepPar0 <- c(stepMu0, stepSig0)

anglePar0 <- c(9, 6, 0.5) #refer to stateNames. initial concentration parameters; little variation in angle, intermediate, high

# MODEL 1
# Fit first model without covariates
PIGU_HMM_1 <- fitHMM(data = tracks, nbStates = 3, 
                     dist = list(step="gamma", angle="vm"), 
                     Par0 = list(step=stepPar0, angle=anglePar0),
                     estAngleMean = list(angle=FALSE),
                     stateNames = stateNames)

# MODEL 2
# Formula for transition probabilities; add TOD
formula1 <- ~cosinor(tod, period = 24)
# Get parameter estimates from previous model to initialize next model
Par0_2 <- getPar0(model = PIGU_HMM_1, 
                  formula = formula1)

# Use parameter estimates and coefficient from previous model; include TOD covariate
PIGU_HMM_2 <- fitHMM(data = tracks, nbStates = 3, dist = list(step="gamma", angle="vm"),
                     Par0 = list(step=Par0_2$Par$step, angle=Par0_2$Par$angle), 
                     beta0 = Par0_2$beta,
                     estAngleMean = list(angle=FALSE),
                     stateNames = stateNames,
                     formula = formula1)


# MODEL VALIDATION
states <- viterbi(PIGU_HMM_1) #decode most likely state sequence; change out HMM models
table(states)/nrow(tracks) #derive percentage of time spent in each state

# Compare models
AIC(PIGU_HMM_1, PIGU_HMM_2)

# Psuedo-residuals for steps and angles
pr <- pseudoRes(PIGU_HMM_1)

# Plot ACF of step pseudo-residuals
acf(pr$stepRes[!is.na(pr$stepRes)],lag.max = 300)

# Save the output of best model!
save(PIGU_HMM_1, file = here("results", "PIGU_HMM_1.RData"))
#############################################
