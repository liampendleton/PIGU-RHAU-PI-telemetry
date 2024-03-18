library(here)
data.22 <- read.csv(here("data", "22_tags_processed_3.10", "Obs070722_105346_Tag44067.csv"))


test <- do.call(rbind, lapply(list.files(path = here("data", "22_tags_processed_3.10"), pattern = '\\.csv'), read.csv))


myMergedData <- do.call(rbind, lapply(list.files(path = "N:/Ring data by cruise"), read.csv))
          
          