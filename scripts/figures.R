library(dplyr)
library(momentuHMM)
library(here)
library(ggplot2)
library(cowplot)
library(patchwork)

#################
### PIGU DATA ###
#################
# CREATE SAT MAPS OF ORIGINAL TRACKING DATA
PIGU_data_LL <- read.csv(here("data", "PIGU_data", "PIGU_data_latlon.csv"))
colnames(PIGU_data_LL) <- c("ID", "y", "x", "time")

PIGU44067 <- PIGU_data_LL %>% filter(ID == 44067)
PIGU44072 <- PIGU_data_LL %>% filter(ID == 44072)
PIGU44372 <- PIGU_data_LL %>% filter(ID == 44372 & y > 48.0) #removing erroneous points that were removed in analysis 
PIGU44505 <- PIGU_data_LL %>% filter(ID == 44505)
PIGU45657 <- PIGU_data_LL %>% filter(ID == 45657)
PIGU45658 <- PIGU_data_LL %>% filter(ID == 45658)

plotSat(PIGU44067, zoom=12, col = "white", ask=FALSE) #run through each ID


# CREATE MAPS OF HMM OUTPUT
load(file=here("results", "PIGU_HMM_1.RData"))

plot(x=PIGU_HMM_1, animals="44067", col=c(rgb(0.5, 0.75, 0.93, 0.2), rgb(1, 0.5, 0, 0.2), "red2"), cex=1.2, lwd=2) #change final digit of first two colors to 0.2 for map, 1 for histograms
plot(x=PIGU_HMM_1, animals="44072", col=c(rgb(0.5, 0.75, 0.93, 0.2), rgb(1, 0.5, 0, 0.2), "red2"), cex=1.2, lwd=2) #change final digit of first two colors to 0.2 for map, 1 for histograms
plot(x=PIGU_HMM_1, animals="44372", col=c(rgb(0.5, 0.75, 0.93, 0.2), rgb(1, 0.5, 0, 0.2), "red2"), cex=1.2, lwd=2) #change final digit of first two colors to 0.2 for map, 1 for histograms
plot(x=PIGU_HMM_1, animals="44505", col=c(rgb(0.5, 0.75, 0.93, 0.2), rgb(1, 0.5, 0, 0.2), "red2"), cex=1.2, lwd=2) #change final digit of first two colors to 0.2 for map, 1 for histograms
plot(x=PIGU_HMM_1, animals="45657", col=c(rgb(0.5, 0.75, 0.93, 0.2), rgb(1, 0.5, 0, 0.2), "red2"), cex=1.2, lwd=2) #change final digit of first two colors to 0.2 for map, 1 for histograms
plot(x=PIGU_HMM_1, animals="45658", col=c(rgb(0.5, 0.75, 0.93, 0.2), rgb(1, 0.5, 0, 0.2), "red2"), cex=1.2, lwd=2) #change final digit of first two colors to 0.2 for map, 1 for histograms

# Retrieve RGB components of "yellow4"
col_info <- col2rgb("skyblue2") / 255

# Create composite step and angle plots
plot(PIGU_HMM_1, plotCI=TRUE, breaks=25)

################################################################################
#################
### RHAU DATA ###
#################
# CREATE SAT MAPS OF ORIGINAL TRACKING DATA
RHAU_data_LL <- read.csv(here("data", "RHAU_data", "RHAU_data_latlon.csv"))
colnames(RHAU_data_LL) <- c("ID", "y", "x", "time")

RHAU44149 <- RHAU_data_LL %>% filter(ID == 44149)
RHAU45653 <- RHAU_data_LL %>% filter(ID == 45653)
RHAU45659 <- RHAU_data_LL %>% filter(ID == 45659)
RHAU45663 <- RHAU_data_LL %>% filter(ID == 45663)
RHAU45672 <- RHAU_data_LL %>% filter(ID == 45672)
RHAU45695 <- RHAU_data_LL %>% filter(ID == 45695)
RHAU45701 <- RHAU_data_LL %>% filter(ID == 45701)

plotSat(RHAU45701, zoom=11, col = "white", ask=FALSE) #run through each ID

# CREATE MAPS OF HMM OUTPUT
load(file=here("results", "RHAU_HMM_1.RData"))

# CREATE MAPS OF HMM OUTPUT
load(file=here("results", "RHAU_HMM_1.RData"))
plot(x=RHAU_HMM_1, animals="44149", col=c(rgb(0.5, 0.75, 0.93, 0.2), rgb(1, 0.5, 0, 0.2), "red2"), cex=1.2, lwd=2) #change final digit of first two colors to 0.2 for map, 1 for histograms
plot(x=RHAU_HMM_1, animals="45663", col=c(rgb(0.5, 0.75, 0.93, 0.2), rgb(1, 0.5, 0, 0.2), "red2"), cex=1.2, lwd=2) #change final digit of first two colors to 0.2 for map, 1 for histograms
plot(x=RHAU_HMM_1, animals="45672", col=c(rgb(0.5, 0.75, 0.93, 0.2), rgb(1, 0.5, 0, 0.2), "red2"), cex=1.2, lwd=2, legend.pos="bottomleft") #change final digit of first two colors to 0.2 for map, 1 for histograms
plot(x=RHAU_HMM_1, animals="45695", col=c(rgb(0.5, 0.75, 0.93, 0.2), rgb(1, 0.5, 0, 0.2), "red2"), cex=1.2, lwd=2) #change final digit of first two colors to 0.2 for map, 1 for histograms
plot(x=RHAU_HMM_1, animals="45701", col=c(rgb(0.5, 0.75, 0.93, 0.2), rgb(1, 0.5, 0, 0.2), "red2"), cex=1.2, lwd=2) #change final digit of first two colors to 0.2 for map, 1 for histograms

# Create composite step and angle plots
plot(RHAU_HMM_1, plotCI=TRUE, breaks=25)

################################################################################

track44067 <- tracks %>%
  filter(ID == 44067) %>%
  mutate(state = states[1:length(time)])

ggplot(track44067, aes(x = time)) +
  geom_point(aes(y = step, color = "red"), size = 2) +
  geom_point(aes(y = angle, color = "blue"), size = 2, shape = 17) +
  geom_line(aes(y = as.numeric(state), color = "State", lwd = 1)) +
  


p_step <- ggplot(track44067, aes(x = time, y = step)) +
  geom_point(aes(y = step, color = "red"), size = 2) +
  labs(y = "Step Length (meters)") +
  theme(axis.line = element_line())

p_angle <- ggplot(track44067, aes(x = time)) +
  geom_point(aes(y = angle, color = "blue"), size = 2, shape = 17) +
  labs(y = "Angle concentration (k)")

p_states <- ggplot(track44067, aes(x = time)) +
  geom_line(aes(y = as.numeric(state), color = "State", lwd = 0.5)) +
  labs(y = "State Sequence")


wrap_elements(get_plot_component(p_step, "ylab-l")) +
  wrap_elements(get_y_axis(p_step)) +
  wrap_elements(get_plot_component(p_angle, "ylab-l")) +
  wrap_elements(get_y_axis(p_angle)) +
  p_states +
  plot_layout(widths = c(3,1,3,1,40))










library(plotly)
plot_ly(data = track44067) %>% 
  add_lines(x = ~time, y = ~state, 
            color = I("red")) %>% 
  add_markers(x = ~time, y = ~step, yaxis = "y2", 
              color = I("blue")) %>%
  add_markers(x = ~time, y = ~angle, yaxis = "y3", 
           color = I("purple")) %>%
  layout(
    yaxis = list(showline = TRUE, side = "right", color = "red",  title = "State Sequence"),
    yaxis2 = list(showline = TRUE, side = "left", color = "blue", title = "Step Length (meters)"), 
    yaxis3 = list(showline = TRUE, side = "left", anchor = "free", overlaying = "y", color = "purple", title = "Angle Concentration (k)"),
    xaxis = list(showline = FALSE, zeroline = FALSE, dtick = 1, title = "Time"), 
    showlegend = FALSE, 
    margin = list(pad = 50, b = 90, l = 150, r = 90),
    legend = list(orientation = "h"))




# Create the plot with color adjustments
plot_ly(data = track44067) %>%
  add_lines(x = ~time, y = ~state, 
            color = I("black"), name = "State Sequence", line = list(width = 3)) %>%
  add_markers(x = ~time, y = ~step, yaxis = "y2", 
              color = I("red"), name = "Step Length") %>%
  add_markers(x = ~time, y = ~angle, yaxis = "y3", 
              color = I("purple"), name = "Turn Angle") %>%
  layout(
    yaxis = list(showline = TRUE, side = "right", color = "black", title = "State Sequence",
                 tickmode = "array", tickvals = c(1, 2, 3), ticktext = c("1", "2", "3")),
    yaxis2 = list(showline = TRUE, side = "left", color = "red", title = "Step Length (meters)"), 
    yaxis3 = list(showline = TRUE, side = "left", anchor = "free", overlaying = "y2", color = "purple", 
                  title = "Angle (degrees)", position = 0.1),  # Adjust the position if necessary
    xaxis = list(showline = FALSE, zeroline = FALSE, dtick = 1, title = "Time"), 
    showlegend = TRUE, 
    margin = list(pad = 50, b = 90, l = 150, r = 90),
    legend = list(orientation = "h", x = 0.1, y = 1.1))




