data <- read.csv()
blanks_rm_data <- data |>
  filter(is.na("blank 1") == TRUE)
blanks_rm_data <- blanks_rm_data[-"blank 1",]
avg_data <- blanks_rm_data |>
  pivot_longer(values_to = "intensity", cols = all, names_to = "file")
avg_data <- avg_data |>
  mutate("strain" = substr(start = 6, 9)) |>
  mutate("treatment" = substr(start = 10, 11))
avg_data <- avg_data |>
  group_by(strain, treatment) |>
  summarize(avg = mean(intensity))

FC_data <- avg_data |>
  group_by(strain) |>
  mutate("FC" = avg[, "S"]/avg[, "C"])
