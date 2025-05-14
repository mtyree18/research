library(ggplot2)
library(tidyverse)
library(outliers)


combined_data <- read.csv("Amino_Acid_HILIC_Data.csv")
combined_data <- as_tibble(combined_data)
combined_data[is.na(combined_data)] <- 0
combined_data <- combined_data[, -c(8:10)]
combined_data <- combined_data[-c(1:3), ]
combined_data <- pivot_longer(combined_data, cols = -c(Compound, log2_hrs), names_to = "yeast", values_to = "rel_abundance")


#average the triplicate data and convert to non log rel_abundance
  if(max(combined_data[4], na.rm = TRUE) > 55) {
  avg_data <- combined_data |> 
    group_by(Compound, yeast, log2_hrs) |>
    mutate(rel_abundance = (mean(as.numeric(rel_abundance), na.rm = TRUE))) |>
    distinct()
  } else {
  avg_data <- combined_data |> 
    group_by(Compound, yeast, log2_hrs) |>
    mutate(rel_abundance = 2^(mean(as.numeric(rel_abundance), na.rm = TRUE))) |>
    distinct()
  }

#create subset of the averaged data that just includes your compound of interest
  cmpd_avg_data <- avg_data |> 
    ungroup() |>
    filter(
    str_detect(Compound, "^phenylalanine")
    )

# #add a row for each yeast corresponding to the unique po2 row, 
# #then remove that row
#   po2_data <- cmpd_avg_data |>
#     filter(yeast == "po2") |>
#     select(rel_abundance, log2_hrs, Compound)
#   yeast_list <- cmpd_avg_data |>
#     filter(yeast != "po2") |>
#     distinct(yeast)
#   new_rows <- yeast_list |>
#     mutate(rel_abundance = po2_data$rel_abundance,
#     log2_hrs = po2_data$log2_hrs,
#     Compound = po2_data$Compound)
#   cmpd_avg_data <- cmpd_avg_data |>
#     filter(yeast != "po2") |>
#     bind_rows(new_rows)
#   
# #normalize by dividing every data point by that yeasts 
# #corresponding first data point (aka the one we just added)
#   cmpd_avg_data <- cmpd_avg_data |>
#     mutate(rel_abundance = rel_abundance/pull(po2_data, 1))

#convert from log2 hours to log hours
  cmpd_avg_data <- cmpd_avg_data |>
    ungroup() |>
    rename(log_hrs = log2_hrs) |>
    mutate(log_hrs = log(as.numeric(log_hrs))) |>
    filter(log_hrs > 1.3)
  
ref_value <- cmpd_avg_data |>
  filter(log_hrs == 1.595339) |>
  pull(rel_abundance) 
cmpd_avg_data <- cmpd_avg_data |>
  mutate(rel_abundance = rel_abundance/ref_value)

#define label positions for series
xpos <- 5
ypos <- 200
offset <- 10

#plot data
ggplot(cmpd_avg_data, aes(x = log_hrs, y = rel_abundance, color = yeast, group = yeast)) +
    geom_line(linewidth = 1, show.legend = FALSE) +
    geom_point(show.legend = FALSE) +
    labs(title = "Phenylalanine") +
    xlab("Log(Time)") + 
    ylab("Relative Abundance") +
    theme_classic(base_size = 24) +
    theme(
      axis.text = element_text(color = "black"),
      plot.margin = margin(5.5, 60, 5.5, 5.5)
    ) +
    coord_cartesian(clip = "off")+
#label maker
  annotate("text", x = xpos, y = ypos+offset*4, label = "Belgian", color = "#F8766D", hjust = 0, size = 6) +
  annotate("text", x = xpos, y = ypos+offset*3, label = "California", color = "#A3A500", hjust = 0, size = 6) +
  annotate("text", x = xpos, y = ypos+offset*2, label = "English", color = "#00BF7D", hjust = 0, size = 6) +
  annotate("text", x = xpos, y = ypos+offset, label = "Pilsner", color = "#00B0F6", hjust = 0, size = 6) +
  annotate("text", x = xpos, y = ypos, label = "Sake", color = "#E76BF3", hjust = 0, size = 6)
  
  