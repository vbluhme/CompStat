## Importance sampling
## Ruin probabilities
library(Rcpp)
Rcpp::sourceCpp("assignments/2_unisamp/ruin.cpp")

ruin_MC <- function(N_sim, n, S0, X_gen, ...) {
  ruin <- numeric(N_sim)
  X <- X_gen(N_sim * n, ...)
  for (i in 1:N_sim) {
    S <- S0 + cumsum(X[((i-1)*n + 1):(i*n)])
    ruin[i] <- min(S) <= 0
  }
  cumsum(ruin) / 1:N_sim
}

ruin_MC_rcpp <- function(N_sim, n, S0, X_gen, ...) {
  X <- X_gen(N_sim * n, ...)
  ruin <- ruin_cpp(X, S0, n, N_sim)
  cumsum(ruin) / 1:N_sim
}
