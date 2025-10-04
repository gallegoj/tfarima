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
// c(theta_0, theta_1, ..., theta_q). 
// 
// [[Rcpp::export]]
arma::colvec acovtomaC(const arma::colvec &g, int &code, double tol, 
                       int max_iter) {
  int i, j, q, iter;
  code = 0;
  q = g.n_elem;
  if (q == 1) return arma::colvec({sqrt(g(0))});
  
  vec f(q), t0(q, fill::zeros), t1(q), t2(q);
  mat T(q, q);
  
  t0(0) = sqrt(g(0));
  iter = 0;
  while(iter < max_iter) {
    for (i = 0; i < q; i++) {
      for (j = 0; j < q-i; j++) T(i, j) = t0(i+j);
      for (j = q-i; j < q; j++) T(i, j) = 0;
      for (j=i; j < q; j++) T(i, j) += t0(j-i);
      f(i) = -g(i);
      for (j = 0; j < q-i; j++)
        f(i) += t0(j)*t0(i+j);
    }
    if ( all(abs(f) < tol) ) {
      t1 = t0;
      break;
    }
    
    try {
      t1 = t0 - solve(T, f, solve_opts::no_approx);
    } catch (...) {
      code = -1;
      return t0;
    }      
    t0 = t1;
    ++iter;
  }
  
  if (iter == max_iter)
    code = -1;

  //checkmaC(t1);
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
