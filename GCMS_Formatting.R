library(ggplot2)
library(tidyverse)

data <- read.csv("GCMS_Data.csv")
data <- as_tibble(data)
#create a long format table of time values to append later
time_col <- data
#define each column as a time point in new df
colnames(time_col) <- time_col[1, ]
#delete all columns and rows that aren't times
time_col <- time_col[-1, -1]
#convert to long format
time_col <- time_col |>
  pivot_longer(
    cols = everything(),
    names_to = "log2_hours",
    values_to = "values"
  )
#create a long format data table not yet including time
combined_data <- data[-1, ] |>
  pivot_longer(
    cols = -Compound,
    names_to = "yeast",
    values_to = "rel_abundance"
  )
#add time to the data table from col 1 of the time_col df
combined_data <- combined_data |>
  mutate(log2_hrs = pull(time_col, 1))
#if the data is from lcms data it will need this if to transform yeast names
if(str_length(combined_data$yeast[[1]]) > 4){
  combined_data <- combined_data |>
    mutate(yeast = tolower(substr(yeast, 19, 21)))
  #if the data is from gcms data it needs this else condition
} else {
  combined_data <- combined_data |>
    mutate(yeast = substr(yeast, 1, 3))
}

combined_data <- combined_data |>
  rename(hrs = log2_hrs) |>
  mutate(hrs = 2^as.numeric(hrs))



combined_data <- combined_data[, c("Compound", "hrs", "yeast", "rel_abundance")]
write.csv(combined_data, file = "GCMS_Data.csv", append = FALSE, quote = FALSE, row.names = FALSE)



