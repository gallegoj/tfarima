#include "RcppArmadillo.h"
#include "gls.h"
#include "res.h"
using namespace arma;

// Generalized Least Squares
//
// \code{glsC} GLS estimation.
//
// @param y Dependent variable, numeric vector.
// @param X Design matrix, numeric vector.
// @param phi AR polynomial, numeric vector.
// @param theta MA polynomial, numeric vector.
// 
// @return \code{glsC} returns the GLS estimates. 
// 
// [[Rcpp::export]]
const arma::colvec glsC(const arma::colvec &y, const arma::mat &X , const arma::colvec &phi, const arma::colvec &theta) {
  
  int i, j, nc, nr;
  nr = y.n_rows;
  nc = X.n_cols;

  mat XX(nc, nc);
  colvec Xy(nc);
  mat Z(nr, nc); 

  colvec z = condresC(y, phi, theta, true);
  for (i = 0; i < nc; i++) {
    Z.col(i) = condresC(X.col(i), phi, theta, true);
  }
  
  for (i = 0; i < nc; i++) {
    XX(i, i) = dot(Z.col(i), Z.col(i));
    Xy(i) = dot(Z.col(i), z);
    for (j = 0; j < i; j++)
      XX(i, j) = XX(j, i) = dot(X.col(i), X.col(j));
  }
  
  return solve(XX, Xy);

}
