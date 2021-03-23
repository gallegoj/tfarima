#include "RcppArmadillo.h"
#include "outliers.h"
#include "res.h"
#include "diff.h"
#include "pol.h"
#include "gls.h"

using namespace arma;

void tratioC(int tau, arma::colvec &x, const arma::colvec &y, double &b, double &t) {
  int i, N;
  double Sxx, Sxy;
  N = y.n_elem;
  Sxx = Sxy = 0;
  for (i = tau; i < N; i++) {
    Sxy += y(i)*x(i - tau);
    Sxx += pow(x(i - tau), 2);
  }
  b = Sxy/Sxx;
  Sxy = 0;
  for (i = 0; i < tau; i++) Sxy += pow(y(i), 2);
  for (i = tau; i < N; i++) Sxy += pow(y(i) - b*x(i-tau), 2);
  Sxy /= (N - 1); 
  t = b/sqrt(Sxy/Sxx);

}
  

void tempchangeC(int tau, arma::colvec &x, const arma::colvec &y, double &b, 
                 double &t, double &w) {
// t.ratio for temporary change
  int i, n, iter;
  double h, Sxx, Sxy, t0, b1, w1;
  
  t0 = t;
  h = sqrt(2.2e-16);
  n = x.n_elem;
  vec x1(n, fill::zeros), u(n), u1(n);

  w = 0.7;
  iter = 0;  
  while(iter++ < 25) {
    x1(tau) = 1.0;
    Sxx = 1; Sxy = y(tau);
    for (i = tau + 1; i < n; i++){
      x1(i) = w*x1(i-1) + x(i - tau);
      Sxx += pow(x1(i), 2);
      Sxy += y(i)*x1(i);
    }
    b = Sxy/Sxx;
    u = x1*b - y;

    w1 = w + h;
    Sxx = 1; Sxy = y(tau);
    for (i = tau + 1; i < n; i++){
      x1(i) = w1*x1(i-1) + x(i - tau);
      Sxx += pow(x1(i), 2);
      Sxy += y(i)*x1(i);
    }
    b1 = Sxy/Sxx;
    u1 = x1*b1 - y;
    x1 = (u1 - u)/(w1 - w);    
    Sxx = dot(x1, x1);
    Sxy = dot(x1, u);
    if (fabs(Sxy) < 0.0001) break;
    w1 = w - Sxy/Sxx;
    w = w1;
  }

  if (w > 0.3 && w < 1 && iter < 25) {
    x1(tau) = 1.0;
    Sxx = 1; Sxy = y(tau);
    for (i = tau + 1; i < n; i++) {
      x1(i) = w*x1(i-1) + x(i - tau);
      Sxx += pow(x1(i), 2);
      Sxy += y(i)*x1(i);
    }
    x = x1;
    b = Sxy/Sxx;
    w1 = dot(y, y) - pow(b, 2)*Sxx;
    w1 /= (x.n_elem-1);
    t = b/sqrt(w1/Sxx);
    if (fabs(t) < fabs(t0))
      t = w = 0;
  } else t = w = 0;

}

// Automatic outlier detection
//
// \code{outliersC} C function to detect outliers.
//
// @param z the time series.
// @param bc logical. If \code{TRUE} logs are taken.
// @param mu the mean of the stationary series.
// @param phi,nabla,theta numeric vectors containing the coefficients
// of this polynomiall.
// @param timing a integer vector with the indices of outliers. 
// @param c cut point to classify an observation as outlier. 
// 
// @return \code{outliersC} returns a matrix with the detected outliers.  
// 
// @section Warning:
// This C function is called by the R functions \code{\link{outliers.um}} and
//  \code{\link{outliers.tfm}}. 
// 
// [[Rcpp::export]]
arma::mat outliersC(const arma::colvec &z, bool bc, double mu, const arma::colvec &phi,
                    const arma::colvec &nabla, const arma::colvec &theta, 
                    arma::ucolvec &timing, double c) {
  
  int t, tau, i, N, T, k, iter;
  double d1, d2, sa;
  
  N = z.n_elem;
  vec w = diffC(z, nabla, bc);

  vec a = exactresC(w, phi, theta);
  if ((int)a.n_elem > N) a.shed_rows(1, a.n_elem - N);
  else if ((int)a.n_elem < N) a.insert_rows(0, N - a.n_elem);
  vec a1 = a;

  vec ph = polymultC(phi, nabla);
  vec x1 = polyratioC(ph, theta, N-1); // pi-weights
  vec x2 = cumsum(x1); // cumulative pi-weights
  vec x3(N);
  
  if (timing(0) < 1) {
    timing.resize(N);
    for (t = 1; t <= N; t++)
      timing(t-1) = t;
  }
  T = timing.n_elem;
  
  mat A(T, 6, fill::zeros);
  sa = stddev(a);
  iter = 0;
  while(iter++ < 4) {
    k = 0;  
    for (t = 0; t < T; t++) {
      if (A(t, 0) < 0.5) { // outlier has not been handled
        tau = timing(t) - 1;
        d1 = a(tau);
        d2 = d1/sa;
        if (fabs(d2) > c||T < N) {
          A(t, 0) = timing(t);
          A(t, 1) = 1; 
          A(t, 2) = d1; 
          A(t, 3) = d2; 
          if (tau < N-1) {
            tratioC(tau, x1, a, d1, d2);
            if ( fabs(d2) > fabs(A(t, 3)) ) {
              A(t, 1) = 2; 
              A(t, 2) = d1; 
              A(t, 3) = d2; 
            }
            tratioC(tau, x2, a, d1, d2);
            if (fabs(d2) > fabs(A(t, 3)) ) {
              A(t, 1) = 3; 
              A(t, 2) = d1; 
              A(t, 3) = d2; 
            }
            if ((int)A(t, 1) < 3) {
              x3 = x1;
              d1 = A(t, 2); d2 = A(t, 3);
              tempchangeC(tau, x3, a, d1, d2, A(t, 4));
              if ( fabs(d2) > fabs(A(t, 3)) ) {
                A(t, 1) = 4; 
                A(t, 2) = d1; 
                A(t, 3) = d2;
              } else A(t, 4) = 0;
            }
          } else A(t, 1) = 2;
          if (A(t, 1) == 1) {
            a(tau) -= A(t, 2);
          } else if (A(t, 1) == 2) {
            for (i = tau; i < N; i++) 
              a(i) -= A(t, 2)*x1(i-tau);
          } else if(A(t, 1) == 3) {
            for (i = tau; i < N; i++) 
              a(i) -= A(t, 2)*x2(i-tau);
          } else if (A(t, 1) == 4) {
            for (i = tau; i < N; i++) 
              a(i) -= A(t, 2)*x3(i-tau);
          }
          sa = stddev(a);
          ++k;
        }
      }
    }
    if (k == 0) break;  
  }

  if (k==0 && iter==1) return mat(1, 1, fill::zeros);

  for (t = T - 1; t > - 1; t--) {
    if (A(t, 0) < 0.5) A.shed_row(t);
  }
  
  k = A.n_rows;
  if (k > 1) {
    vec x(k, fill::zeros);
    mat XX(k, k, fill::zeros);
    vec Xy(k, fill::zeros);
    vec b(k);
    mat vb(k, k);
    int j = 0;
    for (t = 0; t < N; t++) {
      for (j = 0; j < k; j++) {    
        switch((int)A(j, 1)) {
        case 1:
          if ((int)A(j, 0) == t+1) x(j) = 1;
          else x(j) = 0;
          break;
        case 2:
          if ((int)A(j, 0) >= t+1) x(j) = x1((int)A(j, 0) - t - 1);
          break;
        case 3:
          if ((int)A(j, 0) >= t+1) x(j) = x2((int)A(j, 0) - t - 1);
          break;
        default:
          if ((int)A(j, 0) >= t+1) x(j) = x1((int)A(j, 0) - t - 1) 
            + A(j, 4)*x(j);
        }
      }
      XX += x*x.t(); 
      Xy += x*a1(t);
    }
    XX = inv_sympd(XX);
    b = XX*Xy;
    d1 = dot(a, a) - dot(Xy, XX*Xy);
    vb = (d1/(N-k))*XX;
    for (j = k - 1; j > 0; j--) {
      if (fabs(b(j)/sqrt(vb(j, j))) < c) A.shed_row(j);
    }
  }

    
  return A;
    
}
