#include "RcppArmadillo.h"
#include "acovtoma.h"
using namespace arma;

// MA parameters .
//
// \code{acovtomaC} computes the MA parameters from a vector
// of autocovariances.
//
// @param g Numeric vector, c(g_0, g_1, ..., g_q).
// 
// @return \code{ma} returns a numeric vector 
// c(sig2, theta_1, ..., theta_q). 
// 
// [[Rcpp::export]]
arma::colvec acovtomaC(const arma::colvec &g) {
  int i, j, q;
  q = g.n_elem;
  
  vec f(q), t0(q, fill::zeros), t1(q);
  mat T(q, q);
  
  t0(0) = sqrt(g(0));
  while(true) {
    for (i = 0; i < q; i++) {
      for (j = 0; j < q-i; j++) T(i, j) = t0(i+j);
      for (j = q-i; j < q; j++) T(i, j) = 0;
      for (j=i; j < q; j++) T(i, j) += t0(j-i);
      f(i) = -g(i);
      for (j = 0; j < q-i; j++)
        f(i) += t0(j)*t0(i+j);
    }
    t1 = t0 - solve(T, f);
    if ( max(abs(t1-t0)) < 1e-6 ){
      t1 = t0;
      break;
    }
    t0 = t1;
  }
  
  for (j = 1; j < q; j++)
    t1(j) /= -t1(0);

  return t1;
  
}
