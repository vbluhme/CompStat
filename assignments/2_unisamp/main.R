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



# Leg med importance sampling
# g er tætheden, vi simulerer fra
# f er tætheden, vi vil undgå at simulere fra
# w* = f / g
# rg simulerer fra tætheden g. Udregnet ved at invertere fordelingsfunktionen.
# h er funktionen, vi vil kende middelværdi af
# 
# Det lader faktisk til at virke.
g <- function(X, theta) {
  phi <- 1/theta*(exp(theta*2) - exp(-theta*1.9))
  exp(theta*sum(X))/phi^(length(X))
}

f <- function(X) prod(dunif(X, -1.9, 2))
wstar <- function(X, theta) f(X) / g(X, theta) 

rg <- function(n, theta) {
  Y <- runif(n)
  phi <- 1/theta*(exp(theta*2) - exp(-theta*1.9))
  1/theta*log(Y*phi*theta + exp(-1.9*theta))
}  

h <- function(X, S0) {
  min(cumsum(X)) <= -S0
}

X <- lapply(rep(100, 10^4), rg, theta)
apply_f <- function(X) h(X, 30)*wstar(X, theta)
mean(sapply(X, apply_f))