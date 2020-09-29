source("assignments/2_unisamp/main.R")
library(profvis)
library(microbenchmark)
library(bench)
library(tidyverse)
library(patchwork)
library(plotly)
options(scipen=999)

n <- 100                  # Length of each random walk (same as n in assignment text)
S0 <- 30                  # Initial wealth (= 30 in assignment)
X_min <- -1.9; X_max <- 2 # Support of X ~ Unif(X_min, X_max)
theta <- -0.3             # Theta in Importance Sampling


N_sim <- 1000 # N_sim redefineres i løbet af dokumentet

# First check: Do R and Rcpp implementations return same values?
# For Monte Carlo
set.seed(100)
A <- ruin_MC(N_sim, n, S0, X_min, X_max)
set.seed(100)
B <- ruin_MC_rcpp(N_sim, n, S0, X_min, X_max)

all.equal(A, B) # True

# For Importance Sampling
set.seed(100)
A <- ruin_IS(N_sim, n, S0, X_min, X_max, theta)
set.seed(100)
B <- ruin_IS_rcpp(N_sim, n, S0, X_min, X_max, theta)

all.equal(A, B) # True

### The probability of ruin is approximately...
ruin_IS_rcpp(10^6, n, S0, X_min, X_max, theta)[10^6]

# For MC, confidence intervals given by +- 1.96 * sqrt(p(1-p) / 1:N_sim)
N_sim <- 10^5
tibble(N = 1:N_sim, cummean = ruin_MC_rcpp(N_sim, n, S0, X_min, X_max)) %>% 
  filter(N %% 10 == 0) %>% 
  mutate(p = cummean) %>% 
  ggplot(aes(x = N, y = cummean)) +
  geom_ribbon(
    mapping = aes(
      ymin = pmax(0, cummean - 1.96 * sqrt(p*(1-p) / N)),
      ymax = pmin(1, cummean + 1.96 * sqrt(p*(1-p) / N))
    ), fill = "gray") +
  coord_cartesian(ylim = c(0, 0.005)) +
  geom_line()

# Benchmark: 
N_sim <- 10^3
bench::mark(
  "MC R" = {ruin_MC(N_sim, n, S0, X_min, X_max); 1},
  "MC Rcpp" = {ruin_MC_rcpp(N_sim, n, S0, X_min, X_max); 1},
  "IS R" = {ruin_IS(N_sim, n, S0, X_min, X_max, theta); 1},
  "IS Rcpp" = {ruin_IS_rcpp(N_sim, n, S0, X_min, X_max, theta); 1},
)

pressed <- bench::press(
  N_sim = 10^(1:4),
  {
    bench::mark(
      "MC R" = {ruin_MC(N_sim, n, S0, X_min, X_max); 1},
      "MC Rcpp" = {ruin_MC_rcpp(N_sim, n, S0, X_min, X_max); 1},
      "IS R" = {ruin_IS(N_sim, n, S0, X_min, X_max, theta); 1},
      "IS Rcpp" = {ruin_IS_rcpp(N_sim, n, S0, X_min, X_max, theta); 1}
    )
  }
)

## Autoplot
autoplot(pressed)

## Log-log plot
pressed %>% 
  mutate(Method = attr(expression, "description")) %>% 
  ggplot(aes(x = N_sim, y = median, color = Method)) +
  geom_line() +
  geom_point() +
  scale_x_log10()

pressed %>% 
  mutate(Method = attr(expression, "description")) %>% 
  ggplot(aes(x = N_sim, y = median, color = Method)) +
  geom_line() +
  geom_point() +
  scale_y_continuous()

## Plot of sample paths for cummulative mean
k <- 10         # Number of sample paths simulated
N_sim <- 10^5   # Number of simulations per sample path

cummean <- tibble()
for (i in 1:k) {
  MC <- ruin_MC_rcpp(N_sim, n, S0, X_min, X_max)
  IS <- ruin_IS_rcpp(N_sim, n, S0, X_min, X_max, theta)
  cummean <- cummean %>% bind_rows(tibble(type = "MC", i = i, N = 1:N_sim, cummean = MC))
  cummean <- cummean %>% bind_rows(tibble(type = "IS", i = i, N = 1:N_sim, cummean = IS))
}
cummean %>% 
  filter(N %% 100 == 0) %>%
  ggplot(aes(x = N, y = cummean, col = as.factor(i))) +
  geom_line() +
  facet_wrap(~ type) +
  coord_cartesian(ylim = c(0.001, 0.0025)) +
  theme(legend.position = "none")



#### How do we choose theta? ----
# Warning: Takes a while to run
N_sim <- 50000; k <- 50
cummean <- tibble()
for (i in 1:k) {
  for (theta in seq(-0.3, -0.16, 0.02)) {
    IS <- ruin_IS_rcpp(N_sim, n, S0, X_min, X_max, theta)
    cummean <- cummean %>% bind_rows(tibble(theta = theta, i = i, N = 1:N_sim, cummean = IS))
  }
}

p1 <- cummean %>% 
  filter(N %% 10000 == 0) %>% 
  group_by(theta, N) %>% 
  summarise(sd = sd(cummean)) %>% 
  mutate(text = paste("N =", N, "\n θ =", theta, "\n log sd =", round(log(sd), 3))) %>% 
  ggplot(aes(x = theta, y = sd, col = factor(N), group = N, text = text)) +
  geom_line() +
  geom_point() +
  scale_y_log10()

p2 <- cummean %>% 
  filter(N %% 5000 == 0) %>% 
  # filter(round(100*theta) %% 4 == 0) %>% 
  group_by(theta, N) %>% 
  summarise(sd = sd(cummean)) %>% 
  mutate(text = paste("N =", N, "\n θ =", theta, "\n log sd =", round(log(sd), 3))) %>% 
  ggplot(aes(x = N, y = sd, col = factor(theta), group = theta, text = text)) +
  geom_line() +
  geom_point() +
  scale_y_log10()

p3 <- cummean %>% 
  filter(N %% 10000 == 0) %>% 
  group_by(theta, N) %>% 
  summarise(sd = sd(cummean)) %>% 
  mutate(text = paste("N =", N, "\n θ =", theta, "\n log sd =", round(log(sd), 3))) %>%  
  ggplot(aes(x = N, y = theta, fill = sd, text = text)) +
  geom_tile()

p1 | p2

(p1 / p3) | p2

# Plotly versions
ggplotly(p1, tooltip = "text")
ggplotly(p2, tooltip = "text")
ggplotly(p3, tooltip = "text")
