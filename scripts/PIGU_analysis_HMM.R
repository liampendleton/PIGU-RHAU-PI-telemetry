library(here)
library(momentuHMM)
library(sf)
library(dplyr)
library(tidyr)
#########################################################
## PIGU DATA
# Load PIGU data
PIGU_data <- read.csv(here("data", "PIGU_data", "PIGU_data_UTM_TEST.csv"))
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

PIGU_data_crw <- crwOut$crwPredict[,c(3,6,7:17)] #isolate interpolated data

PIGU_data_crw <- PIGU_data_crw %>%
  group_by(ID) %>%
  arrange(time, .by_group = TRUE)

# Create time of day covariate
PIGU_data_crw$tod <- as.numeric(format(PIGU_data_crw$time, "%H")) + as.numeric(format(PIGU_data_crw$time, "%M"))/60 #decimal hours

# Save interpolated data
# write.table(tracks, file = here("data", "PIGU_data", "PIGU_tracks.csv"), col.names = TRUE, row.names = FALSE, sep = ",", quote = FALSE)

# Finish data
PIGU_data <- PIGU_data_crw

# Fill in missing nest coordinates for each individual post-imputation
PIGU_data <- PIGU_data %>%
  group_by(ID) %>%
  fill(nest_x, nest_y, .direction = "downup") %>%
  ungroup()

# Calculate dist2nest
PIGU_data$dist2nest <- sqrt((PIGU_data$mu.x - PIGU_data$nest_x)^2 + (PIGU_data$mu.y - PIGU_data$nest_y)^2)


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

# tracks <- tracks %>% filter(ID == "44067")

# Get an idea of dist parameter specifications by plotting step length and turn angle
# Step length
hist(tracks$step, breaks = 50, xlim=c(0,1500), xlab="Step length (m)")
summary(tracks$step) #look at mean and median
sd(tracks$step, na.rm = TRUE)
# Turn angle
hist(tracks$angle, breaks = 25, xlab="Turn angle (radians)")
summary(tracks$angle)
sd(tracks$angle, na.rm = TRUE)


#################
# set up states #
stateNames <- c("resting", "exploring", "foraging")


############################
# initial parameter values #
stepMu0 <- c(5, 500, 100) #refer to stateNames. low resting movement, high exploring movement, intermediate foraging movement
stepSig0 <- c(2, 100, 50) #deviations
stepPar0 <- c(stepMu0, stepSig0)

anglePar0 <- c(9, 6, 0.5) #refer to stateNames. initial concentration parameters; little variation in angle, intermediate, high


###################################################
# use dist2nest to determine forced resting state #
tracks$rest <- as.double(tracks$dist2nest < 25)
tracks$rest <- as.factor(tracks$dist2nest < 25)
tracks$rest <- factor(tracks$dist2nest < 25, levels = c(FALSE, TRUE), labels = c("0", "1")) #convert rest to 1/0
knownStates <- ifelse(tracks$rest == "1", 1, NA)


###########
# MODEL 1 #
# Fit first model without covariates, without trying to specify knownStates. If you check object PIGU_HMM_1, it'll give you output of the model
PIGU_HMM_1 <- fitHMM(data = tracks, nbStates = 3, 
                     dist = list(step="gamma", angle="vm"), 
                     Par0 = list(step=stepPar0, angle=anglePar0),
                     estAngleMean = list(angle=FALSE),
                     stateNames = stateNames,
                     # knownStates = knownStates
                     )


###########
# MODEL 2 #
# Formulae from McClintock code
formula <- ~ 0 + rest
formulaDelta <- ~ 0 + rest

# Get parameter estimates from previous model to initialize next model
Par0_2 <- getPar0(model = PIGU_HMM_1,
                  formula = formula)

betaRef <- c(3,2,3) # make state 3 (foraging) the tpm reference for state 1 (resting)

# set up initial regression coefficients; rest0 is "if i'm not resting, how likely am I to do this action?". Rest1 is "if i AM resting..."
beta0 <- matrix(c(-100, 25, -25, NA, -25, NA,
                  100, -100, 25, NA, 25, NA),2,3*(3-1),byrow=TRUE,dimnames=list(c("rest0","rest1"),c("1 -> 1", "1 -> 3", "2 -> 1", "2 -> 3", "3 -> 1", "3 -> 2")))

# Use parameter estimates, regression coefficients, TPM, and other info to make a more informed model
PIGU_HMM_2 <- fitHMM(data = tracks, nbStates = 3, dist = list(step="gamma", angle="vm"),
                     Par0 = list(step=Par0_2$Par$step, angle=Par0_2$Par$angle), 
                     betaRef = betaRef,
                     beta0 = beta0,
                     estAngleMean = list(angle=FALSE),
                     stateNames = stateNames,
                     formula = formula,
                     formulaDelta = formulaDelta,
                     knownStates = knownStates)

###########
# MODEL 3 #

# Tweak regression coefficients based on output of previous model
beta0 <- matrix(c(-75, 25, -25, NA, -15, NA,
                  75, -75, 25, NA, 15, NA),2,3*(3-1),byrow=TRUE,dimnames=list(c("rest0","rest1"),c("1 -> 1", "1 -> 3", "2 -> 1", "2 -> 3", "3 -> 1", "3 -> 2")))

Par0_3 <- getPar0(model = PIGU_HMM_2,
                  formula = formula)
# Things start to break down here. Can't figure out how to hone in on answer.
PIGU_HMM_3 <- fitHMM(data = tracks, nbStates = 3, dist = list(step="gamma", angle="vm"),
                     Par0 = list(step=Par0_3$Par$step, angle=Par0_2$Par$angle),
                     betaRef = betaRef,
                     beta0 = beta0,
                     estANgleMean = list(angle=FALSE),
                     stateNames = stateNames,
                     formula = formula,
                     formulaDelta = formulaDelta,
                     knownStates = knownStates)

############### I suggest just running the full code from McColintock's answer to that guy in the thread and comparing his output to what we get from each model here to see how "rest is handled, and what its output looks like


# MODEL VALIDATION
states <- viterbi(PIGU_HMM_2) #decode most likely state sequence; change out HMM models
# table(states)/nrow(tracks) #derive percentage of time spent in each state
plot(PIGU_HMM_2,plotCI=TRUE,ask=FALSE)
unique(viterbi(PIGU_HMM_2)[which(tracks$rest==0)])
unique(viterbi(PIGU_HMM_2)[which(tracks$rest==1)])

# Compare models
AIC(PIGU_HMM_1, PIGU_HMM_2)

# Psuedo-residuals for steps and angles
pr <- pseudoRes(PIGU_HMM_1)

# Plot ACF of step pseudo-residuals
acf(pr$stepRes[!is.na(pr$stepRes)],lag.max = 300)

# Save the output of best model!
save(PIGU_HMM_1, file = here("results", "PIGU_HMM_1.RData"))
#############################################
