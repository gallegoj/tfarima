#include "RcppArmadillo.h"
#include "sim.h"
#include "tacov.h"
#include "pol.h"

using namespace arma;

// [[Rcpp::export]]
arma::colvec simC(arma::colvec a, const bool bc, const double mu,
                  const arma::colvec &phi, const arma::colvec &nabla,
                  const arma::colvec &theta, const arma::colvec &y0) {
  int N, n, p, d, q, pq, r, s, t, tlag, i, j, k, h;
  double sum;
  
  p = phi.n_elem-1;
  d = nabla.n_elem-1;
  q = theta.n_elem-1;
  N = a.n_elem;
  n = N - d;
  
  r = p;
  if(r<q) r = q;
  s = p;
  if(s>q) s = q;
  pq = (p+q);
  // vec a = randn(n);

  vec w(n);
  vec y(N, fill::zeros);
  
  vec g;            // g_0, ..., g_p-1
  vec psi = polyratioC(theta, phi, q);
  
  mat F(r, pq, fill::zeros);   //
  mat V(pq, pq, fill::zeros);
  vec Fu(r);
  
  // V(u)  
  if (p>0) { // Block (1,1)
    g = tacovC(phi, theta, 1, p-1); // 1) g0, ..., g_p-1;    
    for(h=0; h<p; h++){
      V(h, h) = g(0);
      for (k=0; k<h; k++) {
        V(k, h) = V(h, k) = g(h-k);
      }
    }
  }
  
  if (q>0) { // Block (2,2)
    for(h=0; h<q; h++){
      V(p+h, p+h) = 1.0;
    }
  }
  
  if (p>0 && q>0) { // Blocks (1,2) y (2,1)
    for (k=0; k<s; k++){
      for (h=k; h<s ;h++) {
        i = p-1-h+k;
        j = p+q-1-h;
        V(j, i) = V(i, j) = psi(k);
      }
    }
  }
  
  V = chol(V);

  // F
  for (h=0; h<p; h++)
    for (k=h; k<p; k++)
      F(h, k) = -phi(p-k+h);
  
  for (h=0; h<q; h++)
    for (k=h; k<q; k++)
      F(h, p+k) = theta(q-k+h);
  
  Fu = F*V.t()*randn(pq);

  for (t = n-1; t > -1; t--) {
    sum = a(t);
    for (j=1; j<= q; j++) {
      tlag = t - j;
      if (tlag > -1)
        sum += theta(j)*a(tlag);
      else break;
    }
    a(t) = sum;
  }
  
  for (t = 0; t < r; t++)
    a(t) += Fu(t);
  
  for (t = 0; t <n; t++) {
    sum = 0;
    for (j=1; j <= p; j++) {
      tlag = t-j;
      if (tlag > -1)
        sum -= phi(j)*w(tlag); 
      else break;  
    }
    w(t) = sum + a(t);
  }
  
  for (t = 0; t <n; t++) w(t) += mu;
    
  if ((int)y0.n_elem >= d) {
    if (bc) {
      for (t=0; t < d; t++)
        y(t) = log(y0(t));
    } else {
      for (t=0; t < d; t++)
        y(t) = y0(t);
    }
  }
  
  for (t=d; t<N ; t++) {
    sum = 0;
    for (j=1; j<=d; j++) 
      sum -= nabla(j)*y(t-j);
    y(t) = sum + w(t-d);
  }
  
  if (bc) {
    if ((int)y0.n_elem >= d) d = 0;
    for (t=d; t<N ; t++)
      y(t) = exp(y(t));
  }

  return y;
  
}

