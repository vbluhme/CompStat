source("assignments/2_unisamp/main.R")
library(profvis)
library(microbenchmark)
library(bench)
library(tidyverse)

N_sim <- 10^4             # Total number of simulations
n <- 100                  # Length of each random walk (same as n in assignment text)
S0 <- 30                  # Initial wealth (= 30 in assignment)
X_min = -1.9; X_max = 2   # Support of X ~ Unif(X_min, X_max)

# Do they return same ruin probabilities?
set.seed(100)
A <- ruin_MC(N_sim, n, S0, runif, X_min, X_max)
set.seed(100)
B <- ruin_MC_rcpp(N_sim, n, S0, runif, X_min, X_max)

all.equal(A, B) # yes

### The probability of ruin is approximately...
p <- B[N_sim]
p                         # Mean
sqrt(p*(1-p) / N_sim)     # Standard deviation

qplot(1:N_sim, B) + 
  geom_ribbon(
    mapping = aes(
      ymin = pmax(0, B - 1.96 * sqrt(p*(1-p) / 1:N_sim)),
      ymax = B + 1.96 * sqrt(p*(1-p) / 1:N_sim)
    ), fill = "gray") +
  coord_cartesian(ylim = c(0, 0.01)) +
  geom_line() + 
  geom_point() 

# Benchmark: 
# - Rcpp faster than R implementation. Advantage grows with N_sim.
# - Last line shows just simulation of N_sim * n uniform RVs. This is a lower bound for both functions.
# - Rcpp adds only 1-2 ms after simulation.
bench::mark(
  "Base R" = {ruin_MC(N_sim, n, S0, runif, X_min, X_max); 1},
  "Rcpp" = {ruin_MC_rcpp(N_sim, n, S0, runif, X_min, X_max); 1},
  "runif" = {runif(N_sim * n, -1.9, 2); 1}
)

pressed <- bench::press(
  N_sim = 10^(1:4),
  {
    bench::mark(
      "BaseR" = {ruin_MC(N_sim, n, S0, runif, X_min, X_max); 1},
      "Rcpp" = {ruin_MC_rcpp(N_sim, n, S0, runif, X_min, X_max); 1},
      "runif" = {runif(N_sim * n, -1.9, 2); 1}
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
  