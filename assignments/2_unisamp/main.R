## Importance sampling
## Ruin probabilities
library(Rcpp)
Rcpp::sourceCpp("assignments/2_unisamp/ruin.cpp")

ruin_MC <- function(N_sim, n, S0, X_min, X_max) {
  ruin <- numeric(N_sim)
  X <- runif(N_sim * n, X_min, X_max)
  for (i in 1:N_sim) {
    S <- S0 + cumsum(X[((i-1)*n + 1):(i*n)])
    ruin[i] <- min(S) <= 0
  }
  cumsum(ruin) / 1:N_sim
}

ruin_MC_rcpp <- function(N_sim, n, S0, X_min, X_max) {
  ruin <- cpp_MC(S0, n, N_sim, X_min, X_max)
  cumsum(ruin) / 1:N_sim
}



## Implementation of Importance Sampling
ruin_IS <- function(N_sim, n, S0, X_min, X_max, theta) {
  # Simulation from density g
  rg <- function(n, theta) {
    Y <- runif(n)
    phi <- 1/theta*(exp(theta*X_max) - exp(theta*X_min))
    1/theta*log(Y*phi*theta + exp(X_min*theta))
  }
  g_multi <- function(X, theta) {
    phi <- 1/theta*(exp(theta*X_max) - exp(theta*X_min))
    n <- length(X)
    exp(theta*sum(X))/phi^n
  }
  # f <- function(X) 
  #   all(X >= -1.9 & X <= 2)*1/3.9^(length(X))
  wstar <- function(X, theta)
    (X_max - X_min)^(-length(X)) / g_multi(X, theta) 
  h <- function(X, S0)
    min(cumsum(X)) <= -S0
  
  X <- rg(N_sim * n, theta)
  ruin <- numeric(N_sim)
  for (i in 1:N_sim) {
    Xi <- X[((i-1)*n + 1):(i*n)]
    ruin[i] <- h(Xi, S0) * wstar(Xi, theta)
  }
  cumsum(ruin) / 1:N_sim
}

ruin_IS_rcpp <- function(N_sim, n, S0, X_min, X_max, theta) {
  ruin <- cpp_IS(S0, n, N_sim, theta, X_min, X_max)
  cumsum(ruin) / 1:N_sim
}
