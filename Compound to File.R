library(ggplot2)
library(tidyverse)

data <- read.csv("GCMS_Data.csv")
#convert to tibble for ease of use
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

#eliminates all data triplicates containing points with differences greater than 1
combined_data <- combined_data |>
  group_by(yeast, log2_hrs, Compound) |>
  mutate(
    rel_abundance = sapply(rel_abundance, function(current_rel_abundance) {
      # Count how many other values in the group are within 1
      close_intensities <- rel_abundance[rel_abundance != current_rel_abundance & abs(rel_abundance - current_rel_abundance) <= 20]
      # If 2 or more values are close, keep it. Otherwise, change to 0.
      if (length(close_intensities) >= 2) {
        current_rel_abundance
      } else {
        0
      }
    })
  )
combined_data[combined_data == 0] <- NA
combined_data <- na.omit(combined_data)

#average the triplicate data and convert to non log rel_abundance
avg_data <- combined_data |> 
  group_by(Compound, yeast, log2_hrs) |>
  mutate(rel_abundance = 2^(mean(as.numeric(rel_abundance), na.rm = TRUE)))|>
  distinct()

#create subset of the averaged data that just includes your compound of interest
cmpd_avg_data <- avg_data |> 
  ungroup() |>
  filter(
    str_detect(Compound, "^benzeneacetal")
  )

#add a row for each yeast corresponding to the unique po2 row, 
#then remove that row
po2_data <- cmpd_avg_data |>
  filter(yeast == "po2") |>
  select(rel_abundance, log2_hrs, Compound)
yeast_list <- cmpd_avg_data |>
  filter(yeast != "po2") |>
  distinct(yeast)
new_rows <- yeast_list |>
  mutate(rel_abundance = po2_data$rel_abundance,
         log2_hrs = po2_data$log2_hrs,
         Compound = po2_data$Compound)
cmpd_avg_data <- cmpd_avg_data |>
  filter(yeast != "po2") |>
  bind_rows(new_rows)

#normalize by dividing every data point by that yeasts 
#corresponding first data point (aka the one we just added)
cmpd_avg_data <- cmpd_avg_data |>
  mutate(rel_abundance = rel_abundance/pull(po2_data, 1))

#convert from log2 hours to log hours
cmpd_avg_data <- cmpd_avg_data |>
  ungroup() |>
  rename(log_hrs = log2_hrs) |>
  mutate(log_hrs = log(2^as.numeric(log_hrs)))

write.csv(cmpd_avg_data[,-1], "phenylacetaldehyde_GCMS.csv", 
          quote = FALSE, row.names = FALSE)