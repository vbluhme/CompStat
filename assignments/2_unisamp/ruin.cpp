#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector ruin_cpp(NumericVector X, double S0, int n, int N_sim) {
  NumericVector ruin (N_sim);
  double S;
  for(int i = 0; i < N_sim; i++) {
    S = S0;
    for (int j = 0; j < n; j++) {
      S += X[i*n + j];
      if (S <= 0) {
        ruin[i] = 1;
        break;
      }
    }
  }
  return ruin;
}