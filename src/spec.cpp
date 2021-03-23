#include "RcppArmadillo.h"
#include "spec.h"
using namespace arma;

// Theoretical autocovariancies.
//
// \code{spectrumC} Computes the spectrum of an ARMA(p,q) process.
//
// @param phi Numeric vector, c(1, phi_1, ..., phi_p).
// @param theta Numeric vector, c(1, theta_1, ..., theta_q).
// @param sigma2 Double, variance of the error.
// @param nfreq Integer, number of autocovariances.
// 
// @return \code{spectrumC} returns a numeric vector 
// c(gamma_0, gamma_1, ..., gamma_nlags). 
// 
// [[Rcpp::export]]
arma::mat spectrumC(const arma::colvec &phi, const arma::colvec &theta,
                    double sigma2, int nfreq) {
  int i, j, p, q, N;
  double d, k, w, f, c, s;
  
  p = phi.n_elem-1;
  q = theta.n_elem-1;
  
  if (nfreq < 1) nfreq = 501;

  mat A(nfreq, 2, fill::zeros); // A = (k, f)
  N = 2*(nfreq-1);
  d = 2.0*datum::pi;
  for (i = 0; i < nfreq; i++) {
    k = 1.0*i/N; // frequency from 0 to 0.5
    w = d*k; // from 0 to pi
    c = s = 0;
    for (j = 1; j <= q; j++) {
      c += cos(w*j)*theta(j);
      s += sin(w*j)*theta(j);
    }
    f = pow(1+c,2) + pow(s, 2);
    c = s = 0;
    for (j = 1; j <= p; j++) {
      c += cos(w*j)*phi(j);
      s += sin(w*j)*phi(j);
    }
    f /= (pow(1+c, 2) + pow(s, 2));
    A(i, 0) = k;
    A(i, 1) = f*sigma2/d;
  }

  return A;
 
}

// [[Rcpp::export]]
arma::mat pgramC(const arma::colvec &y, bool cpgram) {
  int t, k, q, n;
  double ak, bk, sum, pi2;
  n = y.n_elem;
  
  if (n%2) q = (n-1)/2;
  else q = n/2;
  mat A(q, 2); // A = (k, p)
  
  sum = 0;
  pi2 = 2*datum::pi;
  for (k = 1; k <= q; k++) {
    ak = 0; bk = 0;
    for (t = 0; t < n; t++) {
      ak += cos(pi2*k*t/n)*y(t);
      bk += sin(pi2*k*t/n)*y(t);
    }
    A(k-1, 0) = 1.0*k/n;
    sum += A(k-1, 1) = (pow(ak, 2) + pow(bk, 2))*2.0/n;
  }
  
  if (n%2 == 0) {
    A(q-1, 1) *= 0.5;
    sum += A(q-1, 1);
  }
  
  if(cpgram) {
    for (k = 1; k < q; k++) 
      A(k, 1) += A(k-1, 1);
  }

  for (k = 0; k < q; k++) 
    A(k, 1) /= sum;
  
  return A;
  
}
