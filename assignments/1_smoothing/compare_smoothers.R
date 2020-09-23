source("assignments/1_smoothing/smoother_loocv.R")
library(bench)
library(tidyverse)

Nuuk <- read_table2("data/nuuk.dat.txt", 
                    col_names = c("Year", 1:12), 
                    na = "-999", 
                    skip = 1) %>% 
  gather(key = "Month", value = "Temperature", -Year, convert = TRUE) %>% 
  mutate(Temperature = Temperature / 10) %>% 
  filter(Year > 1866)
Nuuk_year <- group_by(Nuuk, Year) %>% 
  summarise(
    Temperature = mean(Temperature),
    Median = median(Temperature),
    High = max(Temperature), 
    Low = min(Temperature)
  )
n <- nrow(Nuuk_year)

df <- tibble(x = Nuuk_year$Year, y = Nuuk_year$Temperature)

Y_decomp <- smooth_loocv(data = df, p = 50, decompose = TRUE)$y
Y_nodecomp <- smooth_loocv(data = df, p = 50, decompose = FALSE)$y
range(Y_decomp - Y_nodecomp)

benchmark <- bench::mark(
  "Decomp, LOOCV" = {smooth_loocv(data = df, p = 50, decompose = TRUE); 1},
  "No decomp, LOOCV" = {smooth_loocv(data = df, p = 50, decompose = FALSE); 1},
  "Decomp, λ given" = {smooth_loocv(data = df, p = 50, lambda = 3, decompose = TRUE); 1},
  "No decomp, λ given" = {smooth_loocv(data = df, p = 50, lambda = 3, decompose = FALSE); 1},
iterations = 100)
benchmark
autoplot(benchmark)

mark_lambda_given
autoplot(mark_lambda_given)


profvis(smooth_loocv(data = df, p = 50, decompose = TRUE))
profvis(smooth_loocv(data = df, p = 50, decompose = FALSE))
