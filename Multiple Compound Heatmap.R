library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(tidyverse)

#returns the matrix version of the csv data from excel
data_processing <- function(filename){
  raw_data <- read.csv(filename)
  wide_raw_data <- reshape(raw_data, idvar = "log_hrs", timevar = "yeast", direction = "wide")
  rownames(wide_raw_data) <- round(as.numeric(wide_raw_data[, 1]), 3)
  wide_raw_data <- wide_raw_data[, -1]
  colnames(wide_raw_data) <- gsub("rel_abundance.", "", colnames(wide_raw_data))
  return(as.matrix(wide_raw_data))
}


#populates a list with matrices from however many compound files you want
my_list <- list()
file_list <- c("phenylacetaldehyde_GCMS.csv")
for(file in file_list){
  my_list[[file]] <- data_processing(file)
}

#creates a list of heatmaps in vector form
pre_ht_list <- list()
for(i in 1:length(file_list)){
  col_seq <- seq(min(my_list[[i]], na.rm = TRUE), max(my_list[[i]], na.rm = TRUE), length.out = 10)
  col_fun <- colorRamp2(col_seq, hcl_palette = "Purples 3", reverse = TRUE)
  lgd_at <- col_seq[round(seq(1, length(col_seq), length.out = 3))]
  lgd = Legend(col_fun = col_fun, title = "Response", 
               at = lgd_at, 
               legend_height = unit(5, "cm"), 
               labels_gp = gpar(fontsize = 24), 
               title_gp = gpar(fontsize = 28))
               
  pre_ht_list[[i]] <- Heatmap(my_list[[i]], name = "Response\n", show_column_dend = TRUE, #column_order = sort(colnames(my_list[[i]])),
                    col = col_fun, #colorRamp2(seq(min(my_list[[i]]), max(my_list[[i]]), length = 10), hcl_palette = "Purples 3", reverse = TRUE), 
                    show_row_dend = FALSE, column_title = "  Phenylacetaldehyde", row_order = order(as.numeric(rownames(my_list[[i]]))),
                    row_names_gp = gpar(fontsize = 24), column_names_gp = gpar(fontsize = 24), column_title_gp = gpar(fontsize = 32),
                    heatmap_legend_param = list(legend_height = unit(5, "cm"), at = lgd_at, labels = round(lgd_at, 0), 
                    title_gp = gpar(fontsize = 28), labels_gp = gpar(fontsize = 24)))
}

#creates a heatmap list that is a type specific to the complex heatmap library
ht_list <- NULL
for(i in 1:length(pre_ht_list)){
ht_list <- ht_list + pre_ht_list[[i]]
}

ht_opt$HEATMAP_LEGEND_PADDING = unit(10, "mm")
#makes the plots
draw(ht_list, padding = unit(c(1,1,1,1), "cm"))
