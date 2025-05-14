library(ggplot2)
library(tidyverse)

  data <- read.csv("Phe_HILIC_Data.csv")
  data <- as_tibble(data, na.rm = TRUE)
  data[data == 0] <- 1
  data[, -1] <- log(data[,-1])
  long_data <- data |>
    pivot_longer(
    cols = -log2_hrs,
    names_to = "yeast",
    values_to = "rel_abundance")
  
  long_data[long_data == 0] <- NA
  
  long_data <- long_data |>
    group_by(yeast, log2_hrs) |>
    mutate(
      rel_abundance = sapply(rel_abundance, function(current_rel_abundance) {
        # Count how many other values in the group are within 1
        close_intensities <- rel_abundance[rel_abundance != current_rel_abundance & abs(rel_abundance - current_rel_abundance) <= 1]
        # If 2 or more values are close, keep it. Otherwise, change to 0.
        if (length(close_intensities) >= 2) {
          current_rel_abundance
        } else {
          0
        }
      })
    )
  # long_data[long_data == 0] <- NA
  # long_data <- na.omit(long_data)
  
  avg_data <- long_data |>
    group_by(log2_hrs, yeast) |>
    mutate(rel_abundance = mean(as.numeric(rel_abundance), na.rm = FALSE)) |>
    distinct()

  avg_data <- avg_data |>
    rename(log_hrs = log2_hrs) |>
    mutate(log_hrs = log(2^log_hrs))
  
  #normalize by dividing every data point by that yeasts 
  #corresponding first data point (aka the one we just added)
  avg_data <- avg_data |>
    mutate(rel_abundance = rel_abundance/rel_abundance[log_hrs == min_rank(log_hrs)])
  
xpos <- 6
ypos <- 1.05
offset <- 0.03
  
  ggplot(avg_data, aes(x = log_hrs, y = rel_abundance, color = yeast, group = yeast)) +
    geom_line(linewidth = 1, show.legend = FALSE) +
    geom_point(show.legend = FALSE) +
    # geom_text(
    #   data = tbl_labels, 
    #   aes(x = log2_hrs, y = rel_abundance, label = yeast, color = yeast),
    #   hjust = 0, 
    #   nudge_x = 0.1, 
    #   show.legend = FALSE
    # ) +
    labs(title = "Phenylalanine") +
    xlab("Log(Time)") + 
    ylab("rel_abundance") +
    theme_classic(base_size = 24) +
    theme(
      axis.text = element_text(color = "black"),
      plot.margin = margin(5.5, 60, 5.5, 5.5)
    ) +
    coord_cartesian(clip = "off")+

    annotate("text", x = xpos, y = ypos+offset*4, label = "Belgian", color = "#F8766D", hjust = 0, size = 6) +
    annotate("text", x = xpos, y = ypos+offset*3, label = "California", color = "#A3A500", hjust = 0, size = 6) +
    annotate("text", x = xpos, y = ypos+offset*2, label = "English", color = "#00BF7D", hjust = 0, size = 6) +
    annotate("text", x = xpos, y = ypos+offset, label = "Pilsner", color = "#00B0F6", hjust = 0, size = 6) +
    annotate("text", x = xpos, y = ypos, label = "Sake", color = "#E76BF3", hjust = 0, size = 6)
  