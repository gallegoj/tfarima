#include "RcppArmadillo.h"
#include "diff.h"
using namespace arma;


// [[Rcpp::export]]
arma::colvec diffC(const arma::colvec &z, const arma::colvec &nabla, const bool &bc) {  
  int N, n, d, t, i;
  double x;
  N = z.n_elem;
  d = nabla.n_elem-1;
  n = N - d;
  vec y(N), w(n);
  
  if (bc) {
    for (t=0; t<N; t++) {
      if (z(t) <= 0.0)
        Rcpp::stop("Invalid Box-Cox transformation");
      y(t) = log(z(t));
    }
  } else { 
    y = z;
  }
  
  if (d>0) {
    for (t=0; t<n; t++) {
      x = 0;
      for (i=0; i<=d; i++) 
        x += nabla(i)*y(t+d-i);
      w(t) = x;
    }
    return w;
  } else {
    return y;
  }
}

