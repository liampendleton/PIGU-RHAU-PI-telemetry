library(here)
library(dplyr)

# Get the list of files in the directory
file_list <- list.files(path = here("data", "22_tags_processed_3.10", "CSV"))

# Initialize an empty list to store data frames
dfs <- list()

# Loop through each file and read it, skipping the first 5 lines
for (file in file_list) {
  df <- read.csv(file.path(here("data", "22_tags_processed_3.10", "CSV"), file), header = FALSE, skip = 5)
  names(df) <- NULL  #remove column names
  dfs[[file]] <- df
}




