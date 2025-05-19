library(ggplot2)
library(tidyverse)
library(outliers)
`%||%` <- function(a, b) if (!is.null(a)) a else b

data <- read.csv("GCMS_Data_Transformed.csv")
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

#create subset of the averaged data that just includes your compound of interest
cmpd_avg_data <- avg_data |> 
  filter(
    #    Compound == "2-methylbutyraldehyde")
    str_detect(Compound, "^2 & 3-methyl-1-butanol"))
if(class(data[[1]]) == "character" & data[[2,1]] != "2-methyl-furan ") {
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
  if(nrow(ref_value) > 0) {
    ref_value <- as.double(ref_value[which.max(ref_value[[2]]), 4])
  } else {
    ref_value <- 1
  }
}

cmpd_avg_data <- cmpd_avg_data |>
  mutate(rel_abundance = rel_abundance/ref_value)


#define label positions for series
xpos <- 6
ypos <- 2
offset <- 0.3

cmpd_avg_data$yeast <- factor(cmpd_avg_data$yeast)

#plot data
my_plot <- ggplot(cmpd_avg_data, aes(x = log(as.numeric(hrs)), y = rel_abundance, color = yeast, group = yeast)) +
  geom_line(linewidth = 1, show.legend = FALSE) +
  geom_point(show.legend = FALSE) +
  labs(title = "2 & 3-methyl-1-butanol") +
  xlab("time") + 
  ylab("Relative Abundance") +
  theme_classic(base_size = 24) +
  theme(
    axis.text = element_text(color = "black"),
    plot.margin = margin(5.5, 60, 5.5, 5.5)
  ) +
  coord_cartesian(clip = "off")

plot_data <- ggplot_build(my_plot)
group_colors <- plot_data$data[[1]] |>
  select(group, colour) |>
  distinct()
yeast_levels <- levels(cmpd_avg_data$yeast)
yeast_colors <- setNames(group_colors$colour, yeast_levels)

my_plot +
  #label maker
  annotate("text", x = xpos, y = ypos+offset*5, label = "Pre O2", color = as.character(yeast_colors["po2"] %||% "#AAAAAA"), hjust = 0, size = 6) +
  annotate("text", x = xpos, y = ypos+offset*4, label = "California", color = as.character(yeast_colors["cal"] %||% "#AAAAAA"), hjust = 0, size = 6) +
  annotate("text", x = xpos, y = ypos+offset*3, label = "Belgium", color = as.character(yeast_colors["bel"] %||% "#AAAAAA"), hjust = 0, size = 6) +
  annotate("text", x = xpos, y = ypos+offset*2, label = "English", color = as.character(yeast_colors["eng"] %||% "#AAAAAA"), hjust = 0, size = 6) +
  annotate("text", x = xpos, y = ypos+offset, label = "Pilsner", color = as.character(yeast_colors["pil"] %||% "#AAAAAA"), hjust = 0, size = 6) +
  annotate("text", x = xpos, y = ypos, label = "Sake", color = as.character(yeast_colors["sak"] %||% "#AAAAAA"), hjust = 0, size = 6)

