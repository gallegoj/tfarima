#include "RcppArmadillo.h"
#include "loglik.h"
#include "res.h"
#include "tacov.h"
#include "pol.h"
using namespace arma;

// [[Rcpp::export]]
double ellarmaC(const arma::colvec &w, const arma::colvec &phi, const arma::colvec &theta) {
  int n, p, q, pq, r, s, t, i, j, k, h;
  double sum, ssr, detXX, llf;
  
  n = w.n_elem;
  p = phi.n_elem-1;
  q = theta.n_elem-1;
  r = p;
  if(r<q) r = q;
  s = p;
  if (s>q) s = q;
  
  pq = (p+q);
  
  vec a = condresC(w, phi, theta);     // Conditional residuals  
  vec g;                         // g_0, ..., g_p-1
  vec psi = polyratioC(theta, phi, q);

  mat F(r, pq, fill::zeros);   //
  mat FVF(r, r, fill::zeros);    
  mat aux(pq, pq, fill::zeros);

  vec Xt(r, fill::zeros);     // t-th row matrix X
  mat XX(r, r, fill::zeros);    
  vec Xy(r, fill::zeros);
  mat ssr0(1, 1);
  
  // aux = V(u)  
  if (p>0) { // Block (1, 1)
    g = tacovC(phi, theta, 1, p-1); // 1) g0, ..., g_p-1;    
    for (h=0; h<p; h++) {
      aux(h, h) = g(0);
      for (k=0; k<h; k++) {
        aux(k, h) = aux(h, k) = g(h-k);
      }
    }
  }
  
  if (q>0) { // Block (2, 2)
    for (h=0; h<q; h++) {
      aux(p+h, p+h) = 1.0;
    }
  }
  
  if (p>0 && q>0) { // Blocks (1, 2) y (2, 1)
    for (k=0; k<s; k++) {
      for (h=k; h<s ; h++) {
        i = p-1-h+k;
        j = p+q-1-h;
        aux(j, i) = aux(i, j) = psi(k);
      }
    }
  }
  
  // F
  for (h=0; h<p; h++)
    for (k=h; k<p; k++)
      F(h, k) = -phi(p-k+h);

  for (h=0; h<q; h++)
    for (k=h; k<q; k++)
          F(h, p+k) = theta(q-k+h);

  FVF = F*aux*F.t(); 
  if (q>0) {
    ssr = 0;
    Xt(0) = 1.0;
    for (t=0; t<n; t++) {
      ssr += pow(a(t), 2);
      for (i=0; i<r; i++) {
        for (j=0; j<r; j++) {
          XX(i, j) += Xt(i)*Xt(j);
        }
        Xy(i) += Xt(i)*a(t);
      }
      
      sum = 0;
      for (h=0; h<q; h++)
        sum -= theta(h+1)*Xt(h);

      for (j=r-1; j>=1; j--)
        Xt(j) = Xt(j-1);
      Xt(0) = sum;
    }
  } else {
    for (t=0; t<p; t++) {
      Xy(t) = a(t);
      XX(t, t) = 1.0;
    }
    
    ssr = 0;
    for (t=0; t<n; t++) {
      ssr += a(t)*a(t);
    }
  }
  
  aux = FVF*XX;
  for (i=0; i<r; i++)
    aux(i, i) += 1.0;

  detXX = det(aux);
  aux = inv(aux);

  ssr0 = Xy.t()*aux*FVF*Xy;
  ssr -= ssr0(0, 0);
  
  llf = -0.5*n*( 1.0 + log(2*datum::pi) + log(ssr/n) + log(detXX)/n );
      
  return llf;
  
}

// Generalized residuals of an ARMA model
//
// \code{gresC} Compute the generalized residuals of an ARMA model.
//
// These residuals allow the evaluation and maximization of the exact 
// likelihood function of an ARMA process in terms of the residual sum 
// of squares. Besides, they allow compute the variance of the estimations
// from the Jacobian matrix. See, e.g., Nicholls and Hall (1979) and
// Gallego (2009).
//
// @param w Stationary time series, numeric vector.
// @param phi AR polynomial, numeric vector.
// @param theta MA polynomial, numeric vector.
// 
// @return \code{gresC} Generalized residual, column vector.  
// 
// [[Rcpp::export]]
arma::colvec gresC(const arma::colvec &w, const arma::colvec &phi, const arma::colvec &theta) {
  int n, p, q, pq, r, s, t, i, j, k, h;
  double sum;
  
  n = w.n_elem;
  p = phi.n_elem-1;
  q = theta.n_elem-1;
  r = p;
  if(r<q) r = q;
  s = p;
  if (s>q) s = q;
  
  pq = (p+q);
  
  vec a = condresC(w, phi, theta);     // Conditional residuals  
  vec g;                         // g_0, ..., g_p-1
  vec psi = polyratioC(theta, phi, q);
  
  mat F(r, pq, fill::zeros);
  mat V(pq, pq, fill::zeros); 
  mat iP(pq, pq, fill::zeros); // Cholesky factor of inv(V)  
  vec y(n+pq, fill::zeros);
  vec Xt(r, fill::zeros);     // t-th row matrix X
  mat Z(n+pq, pq, fill::zeros);
  mat iZZ(pq, pq, fill::zeros);
  vec b(pq);
  
  if (p>0) { // Block (1, 1)
    g = tacovC(phi, theta, 1, p-1); // 1) g0, ..., g_p-1;    
    for (h=0; h<p; h++) {
      V(h, h) = g(0);
      for (k=0; k<h; k++) {
        V(k, h) = V(h, k) = g(h-k);
      }
    }
  }
  
  if (q>0) { // Block (2, 2)
    for (h=0; h<q; h++) {
      V(p+h, p+h) = 1.0;
    }
  }
  
  if (p>0 && q>0) { // Blocks (1, 2) y (2, 1)
    for (k=0; k<s; k++) {
      for (h=k; h<s ; h++) {
        i = p-1-h+k;
        j = p+q-1-h;
        V(j, i) = V(i, j) = psi(k);
      }
    }
  }
  
  // F
  for (h=0; h<p; h++)
    for (k=h; k<p; k++)
      F(h, k) = -phi(p-k+h);
  
  for (h=0; h<q; h++)
    for (k=h; k<q; k++)
      F(h, p+k) = theta(q-k+h);

  if (p>0)
    iP = inv(chol(V, "lower"));
  else
    iP = V;
  
  if (q>0) {
    Xt(0) = 1.0;
    for (i = 0; i< pq; i++) 
      Z(0, i) = F(0, i);  
    for (t=1; t<n; t++) {
      sum = 0;
      for (h=0; h<q; h++)
        sum -= theta(h+1)*Xt(h);
      for (j=r-1; j>=1; j--)
        Xt(j) = Xt(j-1);
      Xt(0) = sum;
      for (i = 0; i< pq; i++) {
        sum = 0;
        for (j = 0; j < r; j++) 
          sum += Xt(j) * F(j, i);
        Z(t, i) = sum;
      }
    }
  } else {
    for (t=0; t<p; t++) 
      for (j=0; j< p; j++)
      Z(t, t) = F(t, j);
  }
  
  for (t = 0; t < n; t++) 
    y(t) = a(t);

  for (i = 0; i < pq; i++)
    for (j = 0; j <= i; j++)
      Z(n+i, j) = iP(i, j);
  iZZ = Z.t()*Z;
  sum = pow(sqrt(det(iZZ)*det(V)), 1.0/n);
  iZZ = inv(iZZ);
  b = iZZ*Z.t()*y;
  y -= Z*b;
  y *= sum;
  
  return y;
  
}

// [[Rcpp::export]]
double ssrC(const arma::colvec &w, const arma::colvec &phi, const arma::colvec &theta) {
  int n, p, q, pq, r, s, t, i, j, k, h;
  double sum, ssr;
  
  n = w.n_elem;
  p = phi.n_elem-1;
  q = theta.n_elem-1;
  r = p;
  if(r<q) r = q;
  s = p;
  if (s>q) s = q;
  
  pq = (p+q);
  
  vec a = condresC(w, phi, theta);     // Conditional residuals  
  vec g;                         // g_0, ..., g_p-1
  vec psi = polyratioC(theta, phi, q);
  
  mat F(r, pq, fill::zeros);   //
  mat FVF(r, r, fill::zeros);    
  mat aux(pq, pq, fill::zeros);
  
  vec Xt(r, fill::zeros);     // t-th row matrix X
  mat XX(r, r, fill::zeros);    
  vec Xy(r, fill::zeros);
  mat ssr0(1, 1);
  
  // aux = V(u)  
  if (p>0) { // Block (1, 1)
    g = tacovC(phi, theta, 1, p-1); // 1) g0, ..., g_p-1;    
    for (h=0; h<p; h++) {
      aux(h, h) = g(0);
      for (k=0; k<h; k++) {
        aux(k, h) = aux(h, k) = g(h-k);
      }
    }
  }
  
  if (q>0) { // Block (2, 2)
    for (h=0; h<q; h++) {
      aux(p+h, p+h) = 1.0;
    }
  }
  
  if (p>0 && q>0) { // Blocks (1, 2) y (2, 1)
    for (k=0; k<s; k++) {
      for (h=k; h<s ; h++) {
        i = p-1-h+k;
        j = p+q-1-h;
        aux(j, i) = aux(i, j) = psi(k);
      }
    }
  }
  
  // F
  for (h=0; h<p; h++)
    for (k=h; k<p; k++)
      F(h, k) = -phi(p-k+h);
  
  for (h=0; h<q; h++)
    for (k=h; k<q; k++)
      F(h, p+k) = theta(q-k+h);
  
  FVF = F*aux*F.t(); 
  
  if (q>0) {
    ssr = 0;
    Xt(0) = 1.0;
    for (t=0; t<n; t++) {
      ssr += pow(a(t), 2);
      for (i=0; i<r; i++) {
        for (j=0; j<r; j++) {
          XX(i, j) += Xt(i)*Xt(j);
        }
        Xy(i) += Xt(i)*a(t);
      }
      
      sum = 0;
      for (h=0; h<q; h++)
        sum -= theta(h+1)*Xt(h);
      
      for (j=r-1; j>=1; j--)
        Xt(j) = Xt(j-1);
      Xt(0) = sum;
    }
  } else {
    for (t=0; t<p; t++) {
      Xy(t) = a(t);
      XX(t, t) = 1.0;
    }
    
    ssr = 0;
    for (t=0; t<n; t++) {
      ssr += a(t)*a(t);
    }
  }
  
  aux = FVF*XX;
  for (i=0; i<r; i++)
    aux(i, i) += 1.0;
  
  aux = inv(aux);
  
  ssr0 = Xy.t()*aux*FVF*Xy;
  ssr -= ssr0(0, 0);
  
  return ssr;
  
}

// [[Rcpp::export]]
double cllarmaC(const arma::colvec &w, const arma::colvec &phi, const arma::colvec &theta) {
  int n, t;
  double ssr, llf;
  
  n = w.n_elem;
  vec a = condresC(w, phi, theta);     // Conditional residuals  

  ssr = 0;
  for (t=0; t<n; t++)
    ssr += a(t)*a(t);

  llf = -0.5*n*( 1.0 + log(2*datum::pi) + log(ssr/n) );
  
  return llf;
  
}


