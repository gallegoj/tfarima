#include "RcppArmadillo.h"
#include "stsm.h"
using namespace arma;

bool bandchol(const arma::mat &A, arma::mat &L, double &detL)
{
  int m = A.n_rows;
  int n = A.n_cols;
  
  int i, j, k;
  double dot = 0.0;
  
  detL = 1.0;
  for (j = 0; j < m; j++){
    dot = 0.0;
    for(i = 1; i < n; i++)
      dot += L(j,i)*L(j,i);
    
    L(j,0) = A(j,0) - dot;
    if (L(j,0) < 0.0)
      return false;
    
    for (i = 1; i < n; i++) {
      if(j+i < m){
        dot = 0;
        for (k = 1; k < n; k++)
          if(i+k < n)
            dot += L(j+i,i+k)*L(j,k);
        L(j+i,i) = A(j+i,i) - dot;
      }
    }
    L(j,0) = sqrt(L(j,0));
    
    for(i = 1; i < n; i++)
      if(j+i < m)
        L(j+i, i) /= L(j,0);
    detL *= L(j, 0);
  }
  
  if(detL <= 0.0)
    return false;
  return true;
}

void bandcholsol(const arma::mat &L, arma::colvec &x)
{
  int i, k, m, n;
  double sum;
  m = L.n_rows;
  n = L.n_cols;
  
  for(i = 0; i < m; i++){
    sum = x(i);
    for (k=i-1; k>-1; k--)
      if(i-k < n)
        sum -= L(i,i-k)*x(k);
    x(i) = sum/L(i,0);
  }
}

// [[Rcpp::export]]
double llrfC(const arma::colvec &w, const arma::colvec &d, const arma::mat &A,
             const arma::mat &Sv, double s2u, bool s2star) {
  int i, j, k, n;
  double sum, detL, ll, g0, s2;
  n = w.n_elem;
  k = d.n_elem;
  colvec g(k);
  mat G(n, k);
  mat L(n, k);
  colvec e = w;
  
  for (i = 0; i < k; i++) {
    sum = 0;
    for (j = i; j < k; j++) 
      sum += d(j)*d(j-i);
    g(i) = sum;
  }
  g *= s2u;
  for (i = 0; i < k-1; i++) {
    sum = 0;
    for (j = i; j < k-1; j++)
      sum += dot(A.row(j)*Sv, A.row(j-i));
    g(i) += sum;
  }
  g0 = g(0);
  g = g/g0;
  for (j = 0; j < k; j++)
    for (i = 0; i < n; i++)
      G(i, j) = g(j);
  if (!bandchol(G, L, detL)) 
    return datum::nan;
  bandcholsol(L, e);
  s2 = dot(e, e)/(g0*n);
  if (s2star)
    return s2;
  ll = -0.5*n*(1.0 + log(2*datum::pi) + log(s2) + log(g0)) - log(detL);
  return ll;
}
