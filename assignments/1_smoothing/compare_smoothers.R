library(bench)
library(tidyverse)
library(cowplot)

source("assignments/1_smoothing/smoother_loocv.R")

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

# Examples, use of mySmoothing
mySmoothing(data = df, decompose = TRUE) # Use SVD, number of splines as default
mySmoothing(data = df, decompose = FALSE) # Calculate S, number of splines as default
mySmoothing(data = df, p = 50, decompose = TRUE) # Use SVD, number of splines given
mySmoothing(data = df, lambda = 50, decompose = TRUE) # Use SVD, lambda given

## Compare our smoother to smooth.spline ----
Nuuk_year$smooth_R <- with(Nuuk_year, smooth.spline(Year, Temperature, cv = TRUE)$y)
Nuuk_year$smooth_LOOCV <- with(Nuuk_year, mySmoothing(data = tibble(x = Year, y = Temperature))$y)
  
plot1 <- Nuuk_year %>%
  ggplot(aes(x = Year, y = Temperature)) +
  theme(legend.position = "bottom") +
  geom_point() +
  geom_smooth(method = "mySmoothing", se = FALSE, aes(color = "mySmoothing-method")) +
  geom_line(aes(Year, smooth_R, color = "smooth.spline-method")) +
  scale_color_manual("", values = c("blue", "red"))
plot2 <- Nuuk_year %>% 
  ggplot(aes(x = Year, y = smooth_LOOCV - smooth_R)) +
  geom_line() +
  ylab("Numeric difference")
  
plot_grid(plot1, plot2, rel_widths = c(1.5, 1.5))




## COMPARE SPEED BASED ON SIMULATED DATA ----
sim <- function(n, a, b, trend_fun, rand_fun, ...) {
  x <- runif(n, a, b)
  y <- trend_fun(x) + rand_fun(n, ...)
  tibble(x = x, y = y)
}
f <- function(x) cos(x)

# Example of simulation
sim(1000, 0, 10, f, rnorm) %>% ggplot(aes(x, y)) +
  geom_point() +
  geom_smooth(method = "mySmoothing", se = FALSE)

n_list <- seq(10, 500, 50)
median_df <- tibble(n = NA, method = NA, median = NA, ms = NA)
for(n in n_list) {
  df <- sim(n, 0, 10, f, rnorm)
  benchmark <- bench::mark(
    "Decomp, LOOCV" = {mySmoothing(data = df, decompose = TRUE); 1},
    # "No decomp, LOOCV" = {mySmoothing(data = df, decompose = FALSE); 1},
    "Decomp, lambda given" = {mySmoothing(data = df, lambda = 50, decompose = TRUE); 1},
    # "No decomp, lambda given" = {mySmoothing(data = df, lambda = 50, decompose = FALSE); 1},
    "smooth.spline" = {smooth.spline(df$x, df$y, cv = TRUE); 1},
    iterations = 100
  )
  median = as.numeric(benchmark$median)*1000
  ms <- map(1:length(benchmark$min), ~ tibble(ms = as.numeric(benchmark$time[[.]])*1000))
  labs <- attr(benchmark$expression, "description")
  tmp_df <- tibble(n = n, method = labs, median = median, ms = ms)
  median_df <- bind_rows(median_df, tmp_df)
}
median_df <- median_df[-1,]

median_df %>% 
  filter(method %in% c("smooth.spline", "Decomp, LOOCV")) %>% 
  unnest(ms) %>% 
  group_by(n, method) %>% 
  filter(ms <= quantile(ms, 0.95)) %>%  # Filter out 5% slowest for each n
  ungroup() %>% 
  ggplot(aes(x = n, y = ms, col = method)) + 
  geom_jitter(aes(fill = method), size = 1, shape = 21, stroke = 0, height = 0, width = 7) +
  geom_line(aes(x = n, y = median)) +
  scale_colour_manual(values=c("black", "black"))
