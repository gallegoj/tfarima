#include "RcppArmadillo.h"
#include "presample.h"
#include "tacov.h"
#include "pol.h"
using namespace arma;

// [[Rcpp::export]]
arma::mat presampleCovC(const arma::colvec &phi, const arma::colvec &theta, 
                        bool fvf) {
  int p, q, pq, r, s, t, i, j, k, h;

  p = phi.n_elem-1;
  q = theta.n_elem-1;
  r = p;
  if(r < q) r = q;
  s = p;
  if (s>q) s = q;
  pq = (p+q);
  vec g;
  vec psi = polyratioC(theta, phi, q);
  
  mat F(r, pq, fill::zeros);
  mat FVF(r, r, fill::zeros);    
  mat aux(pq, pq, fill::zeros);
  
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
  if (fvf) {
    FVF = F*aux*F.t(); 
    return FVF;  
  } else {
    return aux;
  }
}

