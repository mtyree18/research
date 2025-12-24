library(tidyverse)

#for compounds from the normal phase LCMS data (neg_data and pos_data files) you have to find the mass beforehand in the files because they don't always have the exact mass to the 4th decimal
compound_string <- "phenylethyl alcohol"

data <- read.csv("Transformed for Line Graphs/GCMS_Data_Transformed.csv")
data[[4]][data[[4]] == 0] <- NA
data <- data |>
  group_by(Compound, yeast, hrs) |>
  filter(!all(is.na(across(rel_abundance)))) |>
  ungroup()

#average the triplicate data and convert to non log rel_abundance
if(class(data[[1]]) == "character") {
  avg_data <- data |> 
    group_by(Compound, yeast, hrs) |>
    mutate(rel_abundance = (mean(as.numeric(rel_abundance), na.rm = TRUE))) |>
    distinct() |>
    ungroup()
} else {
  avg_data <- data |> 
    group_by(Compound, yeast, hrs) |>
    mutate(rel_abundance = 2^(mean(as.numeric(rel_abundance)))) |>
    distinct() |>
    ungroup()
}

cmpd_data <- data |>
  filter(Compound == compound_string) |>
  filter(hrs == max(hrs))

tukey_formatted <- cmpd_data |>
  select(3,4)

anova <- aov(rel_abundance ~ yeast, data = tukey_formatted)

TukeyHSD(anova)