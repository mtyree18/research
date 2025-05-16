library(ggplot2)
library(tidyverse)
library(outliers)

data <- read.csv("GCMS_Data.csv")
data[[4]][data[[4]] == 0] <- NA
data <- data |>
  group_by(Compound, yeast, hrs) |>
  filter(!any(is.na(across(rel_abundance)))) |>
  ungroup()

#average the triplicate data and convert to non log rel_abundance
  if(class(data[[1]]) == "character") {
  avg_data <- data |> 
    group_by(Compound, yeast, hrs) |>
    mutate(rel_abundance = (mean(as.numeric(rel_abundance)))) |>
    distinct() |>
    ungroup()
  } else {
  avg_data <- data |> 
    group_by(Compound, yeast, hrs) |>
    mutate(rel_abundance = 2^(mean(as.numeric(rel_abundance)))) |>
    distinct() |>
    ungroup()
  }

#create subset of the averaged data that just includes your compound of interest
  cmpd_avg_data <- avg_data |> 
    filter(
    str_detect(Compound, "^ethyl acetate")
    )
if(class(data[[1]]) == "character" & data[[2,1]] != "ethyl acetate") {
  ref_value <- cmpd_avg_data |>
  group_by(rel_abundance) |>
  mutate(dup_count = n()) |>
  filter(dup_count == 5) |>
  ungroup() |>
  slice_max(order_by = hrs, n = 1)
ref_value <- as.double(ref_value[[1,4]])
} else {
  ref_value <- cmpd_avg_data |>
    filter(yeast == "po2")
  ref_value <- as.double(ref_value[[1,4]])
}

cmpd_avg_data <- cmpd_avg_data |>
  mutate(rel_abundance = rel_abundance/ref_value)
  

#define label positions for series
xpos <- 6
ypos <- 2
offset <- 0.5

#plot data
ggplot(cmpd_avg_data, aes(x = log(hrs), y = rel_abundance, color = yeast, group = yeast)) +
    geom_line(linewidth = 1, show.legend = FALSE) +
    geom_point(show.legend = FALSE) +
    labs(title = "Phenylpyruvic acid") +
    xlab("time") + 
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
  
