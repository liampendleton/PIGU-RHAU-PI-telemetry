
library(tidyverse)
library(sf)
library(here)

Tag44067 <- read_sf(here("data/Tag44067","Obs070722_105346_Tag44067.kml"))

plot(Tag44067)

