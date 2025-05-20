library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(tidyverse)
library(dendsort)
library(seriation)
library(gridExtra)

compound_string <- "1-pentanol"

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

cmpd_avg_data <- avg_data |> 
  filter(str_detect(Compound, compound_string))

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

cmpd_avg_data <- cmpd_avg_data |>
  filter(hrs > 26)

#transform into matrix for heatmap
cmpd_avg_data <- cmpd_avg_data[, -1]
cmpd_avg_data <- cmpd_avg_data |>
  pivot_wider(
    names_from = yeast,
    values_from = rel_abundance
  )
cmpd_avg_data <- cmpd_avg_data |>
  column_to_rownames("hrs") |>
  as.matrix()
#everything before this point was just transforming the data to work for the heatmap package

#establishing the color function
color_function <- colorRamp2(seq(min(cmpd_avg_data, na.rm = TRUE), max(cmpd_avg_data, na.rm = TRUE)-1, length.out = 10), hcl_palette = "Purples 3", reverse = TRUE)
#defining the dendogram object
#this stores the distances between columns
col_dist <- dist(t(cmpd_avg_data), method = "euclidean", na.rm = TRUE)
#this then uses those distances to generate a dendrogram 
#which gets its branches flipped around to minimize the 
#distances between adjacent columns. This variable stores 
#the column order data as well as the hierarchical clustering 
#data which can then be used to generate a dendrogram which is 
#done on the subsequent line
ordering <- seriate(col_dist, method = "OLO")
dendrogram <- as.dendrogram(ordering[[1]])


Heatmap(cmpd_avg_data,
        name = "phenylalanine",
        col = color_function,
        na_col = "black",
        rect_gp = gpar(col = "black", lwd = 2),
        column_title = "Phenylalanine",
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 20, fontface = "bold"),
        row_title = "time(hrs)",
        row_title_side = "right",
        row_title_gp = gpar(fontsize = 20, fontface = "bold"),
        row_title_rot = 0,
        row_order = order(as.numeric(rownames(cmpd_avg_data))),
        cluster_columns = dendrogram,
        column_names_rot = 0,
        column_names_side = "bottom",
        row_names_gp = gpar(fontsize = 15),
        column_names_gp = gpar(fontsize = 15),
        column_names_centered = TRUE,
        row_names_centered = TRUE,
        column_split = 2,
        column_gap = unit(3, "mm"),
        column_dend_side = "top")

