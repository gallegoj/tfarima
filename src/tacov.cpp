#include "RcppArmadillo.h"
#include "pol.h"
#include "tacov.h"
using namespace arma;

// Theoretical autocovariancies.
//
// \code{taucovC} computes some theoretical autocovariances
// for an ARMA(p, q) process.
//
// @param phi Numeric vector, c(1, phi_1, ..., phi_p).
// @param theta Numeric vector, c(1, theta_1, ..., theta_q).
// @param sigma2 Double, variance of the error.
// @param nlags Integer, number of autocovariances.
// 
// @return \code{taucov} returns a numeric vector 
// c(gamma_0, gamma_1, ..., gamma_nlags). 
// 
// [[Rcpp::export]]
arma::colvec tacovC(const arma::colvec &phi, const arma::colvec &theta, 
                   double sigma2, int nlags) {
  int i, j;
  int p, q;
  double x;
  
  p = phi.n_elem-1;
  q = theta.n_elem-1;

  if (nlags<0) nlags *= -1;

  // Cross-covariances G_0, G_1, ..., G_nlags
  // are obtained by solving the SEL Ax=b
  
  vec g(nlags+1, fill::zeros); // g = [g_0, g_1, ..., g_nlags];
  mat A(p, p, fill::zeros); // 
  vec c(q+1, fill::zeros); // E[w_t-k, e_t']
  vec b(p, 1, fill::zeros); // subvector of vec(B) or [vec(B) | zeros()] 
  vec psi, aux;

  // Creating vector c
  psi = polyratioC(theta, phi, q);
  for (i=0; i<=q; i++) {
    x = 0;
    for (j=0; j<=(q-i); j++)
      x += psi(j)*theta(i+j);
    c(i) = x;  
  }
  
  if (p == 0) { // c contains the autocovariances
    if (q<nlags) j = q+1;
    else j = nlags+1;
    for (i=0; i<j; i++)
      g(i) = c(i);
    return g*sigma2;
  }

  // Creating matrix A
  for (i=0; i<p; i++) {
    for (j=1; j<p-i; j++) // Fill diagonal
      A(i+j, j) = phi(i);
    for (j=0; j<=i; j++) // Fill diagonal
      A(i-j, j) += phi(i);
  }
  
  if (p>1) {
    for (j=1; j<p; j++) // Fill diagonal
      A(p-j, j) += phi(p);
  }
  
  for (i=0; i<p; i++) // (Phi_p o Phi_p-i)
    A(0, i) -= phi(p)*phi(p-i);

  if (q>=p) c(0) -= c(p)*phi(p);

  if (q<p) j = q+1;
  else j = p;
  for (i=0; i<j; i++)
    b(i) = c(i);

  aux = solve(A, b); 
  for (i=0; i<p; i++) {
    if (nlags<i)
      break;
    g(i) = aux(i);
  }
  
  for (i=p; i<=nlags; i++) {
    x = 0;
    for (j=1; j<=p; j++)
      x += g(i-j)*phi(j);
    g(i) -= x;
    if (i<=q) g(i) += c(i);
  }
  
  return g*sigma2;

}

// [[Rcpp::export]]
const arma::colvec pacorrC(const arma::colvec &phi, const arma::colvec &theta, 
                int nlags) {
  int i, j;
  double sum1, sum2;
  if (nlags < 1) nlags = 1;
  vec g = tacovC(phi, theta, 1, nlags);
  mat G(nlags, nlags, fill::zeros);
  g /= g(0);
  G(0, 0) = g(1);
  for (i=1; i<nlags; i++) {
    sum1 = sum2 = 0;
    for (j=0; j<i; j++) {
      sum1 += G(i-1, j)*g(i-j);
      sum2 += G(i-1, j)*g(j+1);
    }
    G(i, i) = (g(i+1) - sum1)/(1 - sum2);
    for (j=0; j<i; j++) 
      G(i, j) = G(i-1, j) - G(i, i)*G(i-1, i - j - 1);
  }
  //Rcpp::Rcout << G;
  return diagvec(G);
}


