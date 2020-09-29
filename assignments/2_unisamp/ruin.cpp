#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector cpp_MC(double S0, int n, int N_sim, double X_min, double X_max) {
  NumericVector ruin (N_sim);
  NumericVector Xi (n);
  double S;
  for(int i = 0; i < N_sim; i++) {
    Xi = runif(n, X_min, X_max);
    S = S0;
    for (int j = 0; j < n; j++) {
      S += Xi[j];
      if (S <= 0) {
        ruin[i] = 1;
        break;
      }
    }
  }
  return ruin;
}

// [[Rcpp::export]]
NumericVector cpp_IS(double S0, int n, int N_sim, double theta, double X_min, double X_max) {
  NumericVector Xi (n);
  NumericVector ruin (N_sim);
  double S;
  double phi = 1/theta * (exp(X_max*theta) - exp(X_min*theta));
  double phi_div_xmax_xmin = phi / (X_max - X_min);
  
  for(int i = 0; i < N_sim; i++) {
    Xi = runif(n);
    for (int j = 0; j < n; j++) {
      Xi[j] = 1/theta*log(Xi[j]*phi*theta + exp(X_min*theta));
    }
    
    S = S0;
    for (int j = 0; j < n; j++) {
      S += Xi[j];
      if (S <= 0) {
        ruin[i] = 1;
        break;
      }
    }
    if (ruin[i] == 1) {
      for (int j = 0; j < n; j++) {
        ruin[i] *= phi_div_xmax_xmin / exp(theta*Xi[j]);
      }
    }
  }
  return ruin;
}