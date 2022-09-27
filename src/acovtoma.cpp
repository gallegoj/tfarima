#include "RcppArmadillo.h"
#include "pol.h"
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
// c(sig, theta_1, ..., theta_q). 
// 
// [[Rcpp::export]]
arma::colvec acovtomaC(const arma::colvec &g, int &code) {
  int i, j, q, iter, niter;
  bool again = true;
  code = 0;
  q = g.n_elem;
  if (q == 1) return sqrt(g);
  
  vec f(q), t0(q, fill::zeros), t1(q), t2(q);
  mat T(q, q);
  
  t0(0) = sqrt(g(0));
  t2 = - g;
  iter = 0;
  niter = 100;
  while(iter < niter) {
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
    } else {
      if (max(abs(t1)) < max(abs(t2)))
        t2 = t1;
    }
    t0 = t1;
    ++iter;
    if (iter == niter && again) {
      iter = 0;
      checkmaC(t1);
      again = false;
    }
  }

  if (iter == niter) {
    code = 1;
    t1 = t2;
  }
    
  for (j = 1; j < q; j++)
    t1(j) /= -t1(0);

  return t1;
  
}

void checkmaC(arma::colvec &ma) {
  double s;
  mat A;
  s = ma(0);
  ma *= -1;
  ma(0) = 1;
  A = polyrootsC(ma);
  ma = roots2polC(A, true);
  ma *= -1;
  ma(0) = s;  
}
