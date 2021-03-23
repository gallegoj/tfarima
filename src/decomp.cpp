#include "RcppArmadillo.h"
#include "forecast.h"
#include "diff.h"
#include "res.h"
#include "pol.h"
using namespace arma;

// [[Rcpp::export]]
arma::mat decompHC(const arma::mat &T, const double mu) {
  
  int h, i, j, k, l, r, nr, order;
  double freq, coef;
  
  nr = T.n_rows;
  r = 0;
  for (i=0; i<nr; i++) r += (int)T(i, 5);
  
  if (mu != 0) r += 1;
  mat H(4, r, fill::zeros);
  
  k = 0;
  for (j=0; j<nr; j++) {
    freq = T(j, 3);
    order = (int)T(j, 5);
    if (simeqC(freq, 0)) {
      coef = T(j, 0);
      if (simeqC(coef, 1.0)) {
        l = 0;
        if (mu != 0) order++;
      } else l = 2;
      for (h=0; h<order; h++)
        H(l, k+h) = 1;
      k += order;
    } else {
      coef = T(j, 2);
      if (simeqC(coef, 1.0)) l = 1;
      else l = 3;
      for (h=0; h<order; h++)
        H(l, k+h) = 1;
      k += order;
    } 
  }  
  
  return H;
  
}

// [[Rcpp::export]]
arma::mat decompFC(const arma::mat &T, const double mu) {

  int h, i, j, k, r, nr, order;
  double arg, freq, coef;
  
  nr = T.n_rows;
  r = 0;
  for (i=0; i<nr; i++) r += (int)T(i, 5);
  
  if (mu != 0) r += 1;
  mat F(r+1, r);
  
  k = 0;
  for (j=0; j<nr; j++) {
    freq = T(j, 3);
    order = (int)T(j, 5);
    if (simeqC(freq, 0)) {
      coef = T(j, 0);
      if (simeqC(coef, 1.0)) {
        if (mu != 0) order++;
      }
      for (h=0; h<order; h++) {
        for (i=0; i <= r; i++) 
          F(i, k) = pow(coef, (i+1))*pow(i+1.0, h); 
        k++;
      }
    } else if (simeqC(freq, 0.5)) {
      coef = T(j, 0);
      for (h=0; h<order; h++) {
        for (i=0; i <= r; i++) 
          F(i, k) = pow(coef, (i+1))*pow(i+1.0, h); 
        k++;
      }
    } else if(T(j, 1) > 0.0){
      coef = T(j, 2);
      arg = 2*datum::pi*freq;
      for (h=0; h<order; h++) {
        for (i=0; i <= r; i++) 
          F(i, k) = pow(coef, i+1)*cos(arg*(i+1))*pow(i+1.0, h); 
        k++;
      }
      
      for (h=0; h<order; h++) {
        for (i=0; i <= r; i++) 
          F(i, k) = pow(coef, i+1)*sin(arg*(i+1))*pow(i+1.0, h); 
        k++;
      }
    }
  }  
  
  return F;
  
}

// [[Rcpp::export]]
arma::mat deceffBC(const arma::colvec &y, const bool &bc, const double &mu,
                  const arma::colvec &phi, const arma::colvec &nabla,
                  const arma::colvec &theta, double &sig2, const arma::mat &F,
                  int type) {

  mat B;
  int i, j, p, d, r, N;
  N = y.n_elem;
  p = phi.n_elem - 1;
  d = nabla.n_elem - 1;
  r = p + d;
  
  if (type == 2) {
    vec yc = flipud(y);
    B = deceffBC(yc, bc, mu, phi, nabla, theta, sig2, F, 1);
    B = flipud(B);
    return(B);  
  } else if (type > 2) {
    B = (deceffBC(y, bc, mu, phi, nabla, theta, sig2, F, 1) + 
            deceffBC(y, bc, mu, phi, nabla, theta, sig2, F, 2))/2;
    for (i = 0; i < r; i++) {
      for (j = 0; j < p; j++) {
        B(i, j) *= 2;
        B(N - 1 - i, j) *= 2;
      }
    }
    return B;
  }  
  
  if (mu != 0) r += 1;
  
  if (r != (int)F.n_cols || r+1 != (int)F.n_rows) Rcpp::stop("Wrong matrix F");
  
  mat F1(r, r), F2(r, r), F3(r, r);
  B.zeros(N, r+1);
  mat Y;
  if (bc)
    Y = forecastC(log(y), false, mu, phi, nabla, theta, sig2, r, r);
  else
    Y = forecastC(y, false, mu, phi, nabla, theta, sig2, r, r);

  vec a = exactresC(diffC(y, nabla, bc) - mu, phi, theta);
  if ((int)a.n_elem != N) {
    if ((int)a.n_elem < N) a.insert_rows(0, N - a.n_elem);
    else a.shed_rows(0, a.n_elem - N - 1);
  }
  
  vec psi = polyratioC(theta, polymultC(phi, nabla), r);
  psi.shed_row(0);
  vec yh(r), b(r);
  for(j = 0; j < r; j++) yh(j) = Y(r + j, 0);

  for (j = 0; j < r; j++) {
    for (i = 0; i < r; i++) {
      F1(i, j) = F(i, j);
      F2(i, j) = F(i + 1, j);
    }
  }
  
  b = solve(F1, yh);
  F3 = solve(F1, F2);
  psi = solve(F1, psi);
  for (i = r; i < N; i++) {
    B(i, r) = a(i);
    for (j = 0; j < r; j++)
      B(i, j) = b(j);
    b = F3*b + psi*a(i);
  }

  // Initial values
  
  F3 = solve(F2, F1);
  psi = solve(F2, F1*psi);
  
  for(j = 0; j < r; j++) 
    b(j) = B(r, j);
  
  for (i = r - 1; i > -1; i--) {
    B(i, r) = a(i);
    b = F3*b - psi*a(i);
    for (j = p; j < r; j++)
      B(i, j) = b(j);
  }
  
  return B;
  
}


// [[Rcpp::export]]
arma::vec seasadjC(const arma::colvec &y, const bool &bc, const double &mu,
                    const arma::colvec &phi, const arma::colvec &nabla,
                    const arma::colvec &theta, double &sig2, 
                    const arma::cx_colvec &ariroots, int method) {

  int i, j, p, d, r, N;
  double sum;
  N = y.n_elem;
  p = phi.n_elem - 1;
  d = nabla.n_elem - 1;
  
  if (p+d != (int)ariroots.n_elem) Rcpp::stop("Wrong number of roots");
  
  vec Y(N);
  mat R = sortrootsC(1/ariroots);
  mat H = decompHC(R, mu);
  mat F = decompFC(R, mu);
  mat B = deceffBC(y, FALSE, mu, phi, nabla, theta, sig2, F, method);
  r = B.n_cols - 1;
  vec f(r);

  for (j=0; j<r; j++) f(j) = H(1, j)*F(0, j);  
  if (bc) {
    for (i=0; i<N; i++) {
      sum = 0;
      for (j=0; j<r; j++)
        sum += B(i, j)*f(j);
      Y(i) = exp(log(y(i) - sum));
    }
  } else {
    for (i=0; i<N; i++) {
      sum = 0;
      for (j=0; j<r; j++)
        sum += B(i, j)*f(j);
      Y(i) = y(i) - sum;
    }
  }

  return Y;
}
