#include "RcppArmadillo.h"
#include "res.h"
#include "tacov.h"
#include "pol.h"
using namespace arma;

// Conditional residuals of an ARMA(p, q) model.
//
// \code{condresC} computes the conditional residuals of 
// an ARMA(p, q) model.
//
// The conditional residuals are computed fixing the initial
// conditions  at zero. The AR and MA polynomials are coded as 
// numeric vectors c(1, a_1, ..., a_p).
//
// @param w Stationary time series, numeric vector.
// @param phi AR polynomial, numeric vector.
// @param theta MA polynomial, numeric vector.
// @param forward logical. If TRUE/FALSE, conditional residuals are computed 
// for the forward/backward representation.
// 
// @return \code{condresC} returns a column vector containing the conditional
//  residuals. 
// 
// [[Rcpp::export]]
arma::colvec condresC(const arma::colvec &w, const arma::colvec &phi, 
                      const arma::colvec &theta, const bool forward) {

  int t, j, tlag, n, p, q;
  double x;
  
  n = w.n_elem;
  p = phi.n_elem-1;
  q = theta.n_elem-1;
  
  colvec res(n);
  
  if (forward) {
    if (p>0||q>0) {
      for (t = 0; t < n; t++) {
        x = 0;
        for (j = 0; j<= p; j++) {
          tlag = t-j;
          if (tlag>-1) x +=  phi(j)*w(tlag);
          else break;
        }
        for (j = 1; j<= q; j++) {
          tlag = t-j;
          if (tlag>-1) x -=  theta(j)*res(tlag);
          else break;
        }
        res(t) = x;
      }
    } else
      res = w;
  } else {
    if (p > 0 || q>0) {
      for (t = n-1; t>-1; t--) {
        x = 0;
        for (j = 0; j <= p; j++) {
          tlag = t+j;
          if (tlag<n) x +=  phi(j)*w(tlag);
          else break;
        }
        for (j = 1; j <= q; j++) {
          tlag = t+j;
          if (tlag<n) x -=  theta(j)*res(tlag);
          else break;
        }
        res(t) = x;
      }
    } else
      res = w;
  }
  
  return res;
  
}

// Presample values or initial conditions of an ARMA(p, q) model.
//
// \code{inicondC} computes the initial conditions of an ARMA(p, q) model. 
//
// @inheritParams condresC
// 
// @return \code{inicondC} returns a column vector containing
// the p+q presample values  
// c(w_{-p+1}, w_{-p+1}, ..., w_0, a_{-q+1}, ..., a_0). 
// 
// [[Rcpp::export]]
arma::colvec inicondC(const arma::colvec &w, const arma::colvec &phi, const arma::colvec &theta) {

  int i, j, k, n, p, q, r, s, t;
  double sum;
  
  n = w.n_elem;
  p = phi.n_elem-1;
  q = theta.n_elem-1;
  s = (p+q);
  
  r = p;
  if (r<q) r = q;
  
  vec a = condresC(w, phi, theta);     // Conditional residuals  
  vec g;                               // g0, ..., gp-1
  mat Pu(s, s, fill::zeros);           // Cholesky factor of V(u)
  vec psi = polyratioC(theta, phi, q);
  mat FPu(r, s, fill::zeros);          // 1) F; 2) FPu
  mat aux(s, s, fill::zeros);

  vec xt(r, fill::zeros);              // t-th row matrix X
  mat XX(r, r, fill::zeros);    
  vec Xy(r, 1, fill::zeros);
  vec u(s);                            // initial conditions
  
  // Pu = V(u)  
  if (p > 0) { 
    g = tacovC(phi, theta, 1, p-1);
    for (i = 0; i < p; i++) {            // Block (1, 1) = toeplitz(g)
      Pu(i, i) = g(0);
      for (j = 0; j < i; j++) {
        Pu(i, j) = Pu(j, i) = g(i - j);
      }
    }
  }
  
  if (q > 0) { // Block (2, 2)
    for (i = 0; i < q; i++) {
      Pu(p+i, p+i) = 1;
    }
  }
  
  if (p > 0 && q > 0) { // Blocks (1, 2) y (2, 1)
    k = p;
    if (k > q) k = q;
    for (i = 0; i < k; i++) {
      for (j = i; j < k; j++) {
        Pu(p+q-1-j, p-1-j+i) = Pu(p-1-j+i, p+q-1-j) = psi(i);
      }
    }
  }

  // F
  for (i = 0; i<p; i++)
    for (j = i; j<p; j++)
      FPu(i, j) = -phi(p-j+i);

  for (i = 0; i<q; i++)
    for (j = i; j<q; j++)
      FPu(i, p+j) = theta(q-j+i);
  Pu = chol(Pu, "lower");
  FPu = FPu*Pu;

  if (q>0)
  {
    xt(0) = 1;    
    for (t = 0; t<n; t++) {
      for (i = 0; i<r; i++) {
        for (j = 0; j<r; j++)
          XX(i, j) +=  xt(i)*xt(j);
        Xy(i) +=  xt(i)*a(t);
      }
      sum = 0;
      for (i = 0; i<q; i++)
        sum -=  theta(i+1)*xt(i);
      for (j = r-1; j >= 1; j--)
        xt(j) = xt(j-1);
      xt(0) = sum;
    }
  } else {
    for (t = 0; t<p; t++) {
      Xy(t) +=  a(t);
      XX(t, t) +=  1;
    }
  }
  
  aux = FPu.t()*XX*FPu;
  for (i = 0; i<s; i++)
    aux(i, i) +=  1.0;
  aux = inv_sympd(aux);
  //Rcpp::Rcout << aux << endl;
  Xy = FPu.t()*Xy;
  u = Pu*aux*Xy;

  return u;  
  
}

// Exact residuals of an ARMA(p, q) model.
//
// \code{exactresC} computes the exact residuals of 
// an ARMA(p, q) model.
//
// The exact residuals are computed replacing the presample
// values  by the conditional expectations to the stationary 
// time series. The AR and MA polynomials are coded as 
// numeric vectors c(1, a_1, ..., a_p).
//
// @inheritParams condresC
// 
// @return \code{exactresC} returns a column vector containing the exact
//  residuals. 
// 
// [[Rcpp::export]]
arma::colvec exactresC(const arma::colvec &w, const arma::colvec &phi, const arma::colvec &theta) {
  int t, j, tlag, n, p, q;

  n = w.n_rows;
  p = phi.n_elem-1;
  q = theta.n_elem-1;
  
  vec res(n+q, fill::zeros);
  vec u = inicondC(w, phi, theta);

  j = p;
  for (t = 0; t<q; t++) {
      res(t) = u(j++);
  }

  if (p>0||q>0) {
    for (t = 0; t<n; t++) {
      for (j = 0; j<= p; j++) {
        tlag = t - j;
        if (tlag>-1)
          res(t+q) +=  phi(j)*w(tlag);
        else
          res(t+q) +=  phi(j)*u(p+tlag);
      }
      
      for (j = 1; j<= q; j++) {
        tlag = t-j+q;
        res(t+q) -=  theta(j)*res(tlag);
      }
    }
  } else
    res = w;
  
  return res;
  
}

// [[Rcpp::export]]
double cssrC(const arma::colvec &w, const arma::colvec &phi, const arma::colvec &theta) {
  vec res = condresC(w, phi, theta);
  int n = res.n_elem;
  double ssr = 0;
  for (int t = 0; t <n; t++)
    ssr +=  pow(res[t], 2);
  return(ssr);
}


