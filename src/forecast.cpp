#include "RcppArmadillo.h"
#include "forecast.h"
#include "diff.h"
#include "res.h"
#include "pol.h"
using namespace arma;

// [[Rcpp::export]]
arma::mat forecastC(const arma::colvec &y, const bool bc, const double &mu,
                    const arma::colvec &phi, const arma::colvec &nabla,
                    const arma::colvec &theta, double sig2,
                    int ori, const int hor) {
  int t, h, j, p, d, q, N, nres;
  double sum;

  vec w = diffC(y, nabla, bc) - mu;
  vec a = exactresC(w, phi, theta);
  
  p = phi.n_elem - 1;
  d = nabla.n_elem - 1;
  q = theta.n_elem - 1;
  N = y.n_elem;
  nres = a.n_elem;
  if (ori < 1 || ori < p + d || ori < q  || ori > N)
    Rcpp::stop("Invalid forecast origin");
  
  mat W(ori+hor, 4, fill::zeros); // y, w, a, v(f)
  vec pol = polymultC(phi, nabla);
  vec psi = polyratioC(theta, pol, hor-1);
  psi = cumsum(pow(psi, 2))*sig2;
  
  for (t = 0; t < ori; t++) W(t, 0) = y(t);
  for (t = d; t < ori; t++) W(t, 1) = w(t - d);
  if (nres < N) {
    j = N - nres;
    for (t = j; t < ori; t++) W(t, 2) = a(t-j);
  } else {
    j = nres - N;
    for (t = 0; t < ori; t++) W(t, 2) = a(t+j);
  }
  
  for (h=0; h<hor; h++) {
    sum = 0;
    for (j = 1; j <= p; j++) {
      if(ori+h-j>-1) {
        sum -= phi(j)*W(ori+h-j, 1);
      } else break;
    }
    for (j=0; j<=q; j++) {
      if(ori+h-j>-1) {
        sum += theta(j)*W(ori+h-j, 2);
      } else break;
    }
    W(ori+h, 1) = sum;
  }

  for (h=0; h<hor; h++) {
    sum = W(ori+h, 1) + mu;
    for (j=1; j<=d; j++) {
      if(ori+h-j>-1) {
        if (bc) sum -= nabla(j)*log(W(ori+h-j, 0));
        else sum -= nabla(j)*W(ori+h-j, 0);
      } else break;
    }
    if (bc) W(ori+h, 0) = exp(sum);
    else W(ori+h, 0) = sum;
    W(ori+h, 3) = psi(h);
  }
  
  return(W);
}


// [[Rcpp::export]]
arma::colvec backcastC(const arma::colvec &y, const bool bc, const double &mu,
                    const arma::colvec &phi, const arma::colvec &nabla,
                    const arma::colvec &theta, double sig2,
                    int ori, const int hor) {
  int N = y.n_elem;
  vec z(N);
  for(int i=0; i<N; i++) z(i) = y(N-1-i);
  mat A = forecastC(z, bc, mu, phi, nabla, theta, sig2, N-ori+1, hor); 
  return flipud(A.col(0));
}
