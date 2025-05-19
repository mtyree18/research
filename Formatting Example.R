library(tidyr)
library(ggplot2)

df1 <- data.frame(Compound = c("x", "y", "z"),
                 one = c(1, 0, 3),
                 two = c(1, 1, 2),
                 three = c(2, 3, 1)
                 )
df2 <- data.frame(Compound = c("x", "y", "z"),
                  t1 = c(1, 0, 3),
                  t2 = c(1, 1, 2),
                  t3 = c(2, 3, 1)
)
df1_longer <- df1 |>
  pivot_longer(cols = -Compound, names_to = "yeast", values_to = "values")
df2_longer <- df2 |>
  pivot_longer(cols = -Compound, names_to = "hrs", values_to = "values")