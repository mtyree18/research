library(ggplot2)
library(tidyverse)

data <- read.csv("Amino_Acid_HILIC_Data.csv")
data <- as_tibble(data)
data[is.na(data)] <- 0
data <- data[, -c(8:10)]
data <- pivot_longer(data, cols = -c(Compound, log2_hrs), names_to = "yeast", values_to = "rel_abundance")
data <- data |>
  rename(hrs = log2_hrs)

write.csv(data, file = "Amino_Acid_HILIC_Data.csv", append = FALSE, quote = FALSE, row.names = FALSE)