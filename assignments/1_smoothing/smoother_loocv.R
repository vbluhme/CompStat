library(tidyverse)
library(splines)
library(bench)
library(profvis)

pen_mat <- function(inner_knots) {
  knots <- sort(c(rep(range(inner_knots), 3), inner_knots))
  d <- diff(inner_knots)  # The vector of knot differences; b - a 
  g_ab <- splineDesign(knots, inner_knots, derivs = 2) 
  knots_mid <- inner_knots[-length(inner_knots)] + d / 2
  g_ab_mid <- splineDesign(knots, knots_mid, derivs = 2)
  g_a <- g_ab[-nrow(g_ab), ]
  g_b <- g_ab[-1, ]
  (crossprod(d * g_a,  g_a) + 
      4 * crossprod(d * g_ab_mid, g_ab_mid) + 
      crossprod(d * g_b, g_b)) / 6 
}
decompose_S <- function(Phi, Omega) {
  Phi_svd <- svd(Phi)
  Omega_tilde <- t(crossprod(Phi_svd$v, Omega %*% Phi_svd$v)) / Phi_svd$d
  Omega_tilde <- t(Omega_tilde) / Phi_svd$d
  Omega_tilde_svd <- svd(Omega_tilde)  
  U_tilde <- Phi_svd$u %*% Omega_tilde_svd$u
  
  list(U_tilde = U_tilde, Gamma = Omega_tilde_svd$d)
}

mySmoothing <- function(..., data, p = NULL, lambda = NULL, m = 300, decompose = TRUE) {
  ord <- order(data$x)
  x <- data$x[ord]
  y <- data$y[ord]
  nx <- length(x)
  if (is.null(p)) {
    p <- .nknots.smspl(nx)
  }
  
  inner_knots <- seq(min(x), max(x), length.out = p - 2)
  Phi <- splineDesign(c(rep(range(inner_knots), 3), inner_knots), x)
  Omega <- pen_mat(inner_knots)
  
  if (decompose) {
    decomp <- decompose_S(Phi, Omega)
    decomp$U_tilde_y <- t(decomp$U_tilde) %*% y
    
    loocv <- function(lambda, decomp) {
      f_hat <- decomp$U_tilde_y / (1 + lambda * decomp$Gamma)
      f_hat <- decomp$U_tilde %*% f_hat
      diag <- decomp$U_tilde^2 %*% as.matrix(1 / (1 + lambda * decomp$Gamma))
      sum((y - f_hat)^2 / (1 - diag)^2)
    }
    
    if (is.null(lambda)) {
      lambda <- optimize(loocv, c(0, m), decomp = decomp)$minimum
    }
    
    f_hat <- decomp$U_tilde_y / (1 + lambda * decomp$Gamma)
    f_hat <- decomp$U_tilde %*% f_hat
  } else {
    compute_S <- function(lambda, Phi, Omega)
      Phi %*% solve(crossprod(Phi) + lambda*Omega, t(Phi))
    
    loocv <- function(lambda, Phi, Omega) {
      S <- compute_S(lambda, Phi, Omega)
      f_hat <- S %*% y
      sum((y - f_hat)^2 / (1-diag(S))^2)
    }
    if (is.null(lambda)) {
      lambda <- optimize(loocv, c(0, m), Phi = Phi, Omega = Omega)$minimum
    }
    
    f_hat <- compute_S(lambda, Phi, Omega) %*% y
  }
  
  structure(list(x = x, y = f_hat), class = "mySmoothing")
}
predict.mySmoothing <- function(object, newdata, ...) {
  approx(object$x, object$y, newdata$x)$y
}

