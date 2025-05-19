library(tidyr)
library(ggplot2)

df <- read.csv("GCMS_Data_Unformatted.csv", na.strings = "", header = FALSE)
df <- as_tibble(df)
df1 <- df[-1, ]
df2 <- df[-2, ]
colnames(df1) <- df1[1, ]
colnames(df2) <- df2[1, ]
df1 <- df1[-1, ]
df2 <- df2[-1, ]
colnames(df1)[1] <- "Compound"
colnames(df2)[1] <- "Compound"
df1_longer <- df1 |>
  pivot_longer(cols = -Compound, names_to = "yeast", values_to = "rel_abundance")
df2_longer <- df2 |>
  pivot_longer(cols = -Compound, names_to = "hrs", values_to = "rel abundance")
combined_data <- df1_longer |>
  mutate(hrs = pull(df2_longer, 2))
combined_data <- as.data.frame(combined_data)
combined_data$yeast <- gsub("ExtrF", "sak", combined_data$yeast)
combined_data$yeast <- gsub("Eng", "eng", combined_data$yeast)
combined_data$yeast <- gsub("Cal", "cal", combined_data$yeast)
combined_data$yeast <- gsub("Bel", "bel", combined_data$yeast)
combined_data$yeast <- gsub("Pils", "pil", combined_data$yeast)
combined_data$yeast <- gsub("pre-ferm", "po2", combined_data$yeast)
combined_data <- combined_data[, c(1,4,2,3)]
write.csv(combined_data, file = "GCMS_Data_Transformed.csv", row.names = FALSE)
gcms_data <- read.csv("GCMS_Data_Transformed.csv", sep = ",")
