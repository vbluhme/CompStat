
# Leg med importance sampling
# g er tætheden, vi simulerer fra
# g_multi er den simultane tæthed af n uafh. g-fordelte variable (g_{θ, n}(x) fra opgaven)
# f er tætheden for den uniforme fordeling på (-1.9, 2)^n
# w* = f / g
# rg simulerer fra tætheden g. Udregnet ved at invertere fordelingsfunktionen.
# h er funktionen, vi vil kende middelværdi af
g <- function(x, theta) {
  phi <- 1/theta*(exp(theta*2) - exp(-theta*1.9))
  exp(theta*x)/phi
}
g_multi <- function(X, theta) {
  phi <- 1/theta*(exp(theta*2) - exp(-theta*1.9))
  n <- length(X)
  exp(theta*sum(X))/phi^n
}
rg <- function(n, theta) {
  Y <- runif(n)
  phi <- 1/theta*(exp(theta*2) - exp(-theta*1.9))
  1/theta*log(Y*phi*theta + exp(-1.9*theta))
}

f <- function(X) prod(dunif(X, -1.9, 2))
wstar <- function(X, theta) f(X) / g_multi(X, theta)

h <- function(X, S0) min(cumsum(X)) <= -S0

# We show that we can approximate the probability
theta <- -1/3; S0 <- 30; N_sim <- 10^5; n <- 100
X <- lapply(rep(n, N_sim), rg, theta)
apply_f <- function(X) h(X, S0)*wstar(X, theta)
cummean <- cumsum(sapply(X, apply_f))/ (1:N_sim)
plot(cummean, type = "l"); abline(h = cummean[N_sim], col = "red")

# The simulation function rg() works :)
par(mfrow = c(2,2))
for(theta in c(-1.5, -0.5, 0.5, 1.5)) {
  tmp_g <- function(x) g(x, theta)
  hist(rg(10^5, theta), prob = T, breaks = 80, main = paste("θ =", theta), xlab = "")
  curve(tmp_g, col = "red", add = T)
}
par(mfrow = c(1,1))
