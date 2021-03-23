#include "RcppArmadillo.h"
#include "filter.h"
using namespace arma;

// Filter input 
//
// \code{filterC} filters an input with a rational polynomial.
//
//
// @param x input, numeric vector.
// @param omega polynomial numerator, numeric vector.
// @param delta polynomial denominator, numeric vector.
// @param b delay.
// 
// @return \code{filterC} returns the filtered input. 
// 
// [[Rcpp::export]]
arma::colvec filterC(const arma::colvec &x, const arma::colvec &omega,
                     const arma::colvec &delta, int b){

  int t, j, tlag, n, r, s;
  double d;
  
  n = x.n_elem;
  r = delta.n_elem-1;
  s = omega.n_elem-1;
  b = abs(b);
  
  colvec a(n);
  
  if (r>0||s>0) {
    for (t = 0; t < n; t++) {
      d = 0;
      for (j = 0; j <= s; j++) {
        tlag = t-b-j;
        if (tlag>-1) d += omega(j)*x(tlag);
      }
      for (j = 1; j <= r; j++) {
        tlag = t-j;
        if(tlag>-1) d -= delta(j)*a(tlag);
      }
      a(t) = d;
    }
  }
  else {
    for (t = 0; t < b; t++) 
      a(t) = 0;
    for (t = b; t < n; t++)
      a(t) = omega(0)*x(t-b);
  } 
  
  return a;
  
}
