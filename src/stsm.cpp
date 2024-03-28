#include "RcppArmadillo.h"
#include "pol.h"
#include "res.h"
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
double llrfC(const arma::colvec &w, const arma::colvec &nabla,
             const arma::rowvec &b, const arma::mat &C,
             const arma::mat &S, arma::colvec &s2, bool cform) {

  int h, i, j, k, p, n;
  double sum, detL, ll, g0;
  n = w.n_elem;
  k = b.n_elem;
  colvec phi, psi, g1;
  rowvec varphi(k+1);
  colvec g2(k+1);
  mat I = eye(k, k);
  mat C0 = I;
  mat A(k, k);
  mat Phi0, Theta0, Psi0, G0, Psi;
  mat G(n, k+1, fill::zeros);
  mat L(n, k+1);
  colvec e = w;
  // Reduced form
  varphi(0) = 1;
  A.row(0) = b;
  for (i = 1; i < k; i++) {
    varphi(i) = -trace(C*C0)/i;
    C0 = C*C0 + varphi(i)*I;
    A.row(i) = b*C0;
  }
  varphi(i) = -trace(C*C0)/i;
  if (varphi.n_elem > nabla.n_elem) {
    phi = polydivC(varphi, nabla, false);
    p = phi.n_elem - 1;
    e = condresC(w, phi, colvec(1, fill::ones));
  } else {
    p = 0;
  }
  if (cform)
    A = join_rows(A.t(), colvec(k, fill::zeros));
  else
    A = join_rows(colvec(k, fill::zeros), A.t());
  A = join_cols(varphi, A);
  k += 1;
  // Autocovariances MA part  
  for (i = 0; i < k; i++) {
    sum = 0;
    for (j = i; j < k; j++)
      sum += dot(A.col(j).t()*S, A.col(j-i));
    g2(i) += sum;
  }
  
  if (p > 0) {
    // Psi weights
    Psi = A;
    for (i = 1; i < k; i++ ) {
      for (j = 1; j <= p; j++) {
        if (i-j < 0) break;
        Psi.col(i) -= phi(j)*Psi.col(i-j);
      }
    }
    
    Phi0 = zeros(p, p);
    for (i = 0; i < p; i++) 
      for (j = i; j < p; j++)
        Phi0(i, j-i) = -phi(j+1);
    
    G0 = zeros(p+1, p+1);
    for (i = 0; i <= p; i++) {
      for (j = 0; j <= i; j++)
        G0(i,j) = phi(i-j);
      for (j = i+1; j <= p; j++)
        G0(i,j-i) += phi(j);
    }

    g1 = zeros(p+1);
    for (i = 0; i <= p; i++) {
      sum = 0;
      for (j = i; j < k; j++)
        sum += dot(A.col(j).t()*S, Psi.col(j-i));
      g1(i) += sum;
    }
    g1 = solve(G0, g1);
    g1.shed_row(p);
    G0 = toeplitz(g1);
    G0 = Phi0*G0*Phi0.t();

    C0 = zeros(k, k);
    for (i = 0; i < p; i++)
      for (j = 0; j < p; j++)
        C0(i, j) += G(i,j);

    Theta0.zeros(k-1,k*(k-1));
    for (i = 0; i < k-1; i++) {
      for (j = i; j < k-1; j++)
        for (h = 0; h < k; h++)
          Theta0(i,(j-i)*k+h) = A(j+1, h);
    }
    Psi0.zeros(p,k*(k-1));
    for (i = 0; i < p; i++) {
      for (j = i; j < k-1; j++)
        for (h = 0; h < k; h++)
          Psi0(i, j*k + h) = Psi(j-i, h);
    }
    Phi0 = Phi0*Psi0*kron(eye(k-1,k-1), S)*Theta0.t();

    for (i = 0; i < p; i++) {
      for (j = 0; j < k-1; j++) {
        C0(i,j) += Phi0(i,j);
        C0(j,i) += Phi0(i,j);
      }
    }
  }

  if (p > 0) {
    for (i = 0; i < k; i++)
      for (j = 0; j <= i; j++)
        G(i, i-j) = C0(i,j);
  }

  for (j = 0; j < k; j++)
    for (i = j; i < n; i++)
      G(i, j) += g2(j);
  g0 = G.col(0).max();
  G /= g0;
  if (!bandchol(G, L, detL)) 
    return datum::nan;
  bandcholsol(L, e);
  s2(0) = dot(e, e)/(g0*n);
  ll = -0.5*n*(1 + log(2*datum::pi) + log(s2(0)) + 2*log(detL)/n + log(g0) );
  return ll;
 
}

// [[Rcpp::export]]
arma::colvec resrfC(const arma::colvec &w, const arma::colvec &nabla,
             const arma::rowvec &b, const arma::mat &C,
             const arma::mat &S, arma::colvec& s2, bool cform) {
  
  int h, i, j, k, p, n;
  double sum, detL, ll, g0;
  n = w.n_elem;
  k = b.n_elem;
  colvec phi, psi, g1;
  rowvec varphi(k+1);
  colvec g2(k+1);
  mat I = eye(k, k);
  mat C0 = I;
  mat A(k, k);
  mat Phi0, Theta0, Psi0, G0, Psi;
  mat G(n, k+1, fill::zeros);
  mat L(n, k+1);
  colvec e = w;
  
  // Reduced form
  varphi(0) = 1;
  A.row(0) = b;
  for (i = 1; i < k; i++) {
    varphi(i) = -trace(C*C0)/i;
    C0 = C*C0 + varphi(i)*I;
    A.row(i) = b*C0;
  }
  varphi(i) = -trace(C*C0)/i;
  if (varphi.n_elem > nabla.n_elem) {
    phi = polydivC(varphi, nabla, false);
    p = phi.n_elem - 1;
    e = condresC(w, phi, colvec(1, fill::ones));
  } else {
    p = 0;
  }
  if (cform)
    A = join_rows(A.t(), colvec(k, fill::zeros));
  else
    A = join_rows(colvec(k, fill::zeros), A.t());
  A = join_cols(varphi, A);
  k += 1;
  // Autocovariances MA part  
  for (i = 0; i < k; i++) {
    sum = 0;
    for (j = i; j < k; j++)
      sum += dot(A.col(j).t()*S, A.col(j-i));
    g2(i) += sum;
  }
  
  if (p > 0) {
    // Psi weights
    Psi = A;
    for (i = 1; i < k; i++ ) {
      for (j = 1; j <= p; j++) {
        if (i-j < 0) break;
        Psi.col(i) -= phi(j)*Psi.col(i-j);
      }
    }
    
    Phi0 = zeros(p, p);
    for (i = 0; i < p; i++) 
      for (j = i; j < p; j++)
        Phi0(i, j-i) = -phi(j+1);
    
    G0 = zeros(p+1, p+1);
    for (i = 0; i <= p; i++) {
      for (j = 0; j <= i; j++)
        G0(i,j) = phi(i-j);
      for (j = i+1; j <= p; j++)
        G0(i,j-i) += phi(j);
    }
    
    g1 = zeros(p+1);
    for (i = 0; i <= p; i++) {
      sum = 0;
      for (j = i; j < k; j++)
        sum += dot(A.col(j).t()*S, Psi.col(j-i));
      g1(i) += sum;
    }
    g1 = solve(G0, g1);
    g1.shed_row(p);
    G0 = toeplitz(g1);
    G0 = Phi0*G0*Phi0.t();
    
    C0 = zeros(k, k);
    for (i = 0; i < p; i++)
      for (j = 0; j < p; j++)
        C0(i, j) += G(i,j);
    
    Theta0.zeros(k-1,k*(k-1));
    for (i = 0; i < k-1; i++) {
      for (j = i; j < k-1; j++)
        for (h = 0; h < k; h++)
          Theta0(i,(j-i)*k+h) = A(j+1, h);
    }
    Psi0.zeros(p,k*(k-1));
    for (i = 0; i < p; i++) {
      for (j = i; j < k-1; j++)
        for (h = 0; h < k; h++)
          Psi0(i, j*k + h) = Psi(j-i, h);
    }
    Phi0 = Phi0*Psi0*kron(eye(k-1,k-1), S)*Theta0.t();
    
    for (i = 0; i < p; i++) {
      for (j = 0; j < k-1; j++) {
        C0(i,j) += Phi0(i,j);
        C0(j,i) += Phi0(i,j);
      }
    }
  }
  
  if (p > 0) {
    for (i = 0; i < k; i++)
      for (j = 0; j <= i; j++)
        G(i, i-j) = C0(i,j);
  }
  
  for (j = 0; j < k; j++)
    for (i = j; i < n; i++)
      G(i, j) += g2(j);
  g0 = G.col(0).max();
  G /= g0;
  bandchol(G, L, detL);
  bandcholsol(L, e);
  s2(0) = dot(e, e)/(g0*n);
  return e*pow(detL, 1/n)/sqrt(g0);
}

// [[Rcpp::export]]
bool kf0C(const arma::colvec &y, const arma::colvec &b, const arma::mat &C,
         const arma::mat &S, const arma::colvec &x1, const arma::mat &P1,
         arma::colvec &v, arma::colvec &s2) {
  // STS model must be in future form
  int i, j, k, n, t;
  double sum, s11;
  k = b.n_elem;
  n = y.n_elem;
  if (v.n_elem < n) return false;
  colvec K(k), S12(k), a(k);
  mat S22(k, k), Pa(k, k);
  s11 = S(0,0);
  for (i = 0; i < k; i++) {
    S12(i) = S(i+1, 0);
    for (j = 0; j <= i; j++)
      S22(i, j) = S(i+1, j+1);
  }
  a = x1; Pa = P1;
  for (t = 0; t < n; t++) {
    v(t) = y(t) - dot(b, a);
    s2(t) = dot(b.t()*Pa, b) + s11;
    if (s2(t) < 0)
      return false;
    K = (C*Pa*b + S12)/s2(t);
    a = C*a + K*v(t);
    Pa = C*Pa*C.t() + S22 - (C*Pa*b + S12)*K.t();
    Pa = (Pa + Pa.t())/2;
  }
  return true;
}

// [[Rcpp::export]]
bool kfC(const arma::colvec &y, const arma::colvec &b, const arma::mat &C,
         const arma::mat &S, arma::colvec &x1, arma::mat &P1,
         arma::colvec &v, arma::colvec &s2, arma::mat &X, arma::mat &PX, 
         bool cform, bool filtered, bool xn) {
  // STS model must be in future form
  int i, j, k, n, t;
  double sum, s11;
  k = b.n_elem;
  n = y.n_elem;
  if (v.n_elem < n) return false;
  colvec K(k), S12(k), a(k), a1(k);
  mat S22(k, k), Pa(k, k), Pa1(k,k);
  s11 = S(0,0);
  for (i = 0; i < k; i++) {
    S12(i) = S(i+1, 0);
    for (j = 0; j <= i; j++)
      S22(i, j) = S(i+1, j+1);
  }
  a = x1; Pa = P1;
  for (t = 0; t < n; t++) {
    if (!cform && !filtered) {
      for (i = 0; i < k; i++) {
        X(t, i) = a(i);
        for (j = 0; j < k; j++)
          PX(t*k+i, j) = Pa(i, j);
      }
    }
    v(t) = y(t) - dot(b, a);
    s2(t) = dot(b.t()*Pa, b) + s11; 
    if (!cform && filtered) {
      a1 = Pa*b;
      Pa1 = Pa - a1*a1.t()/s2(t);
      a1 = a + a1*v(t)/s2(t);
      for (i = 0; i < k; i++) {
        X(t, i) = a1(i);
        for (j = 0; j < k; j++)
          PX(t*k+i, j) = Pa1(i, j);
      }
    }
    K = (C*Pa*b + S12)/s2(t);
    a = C*a;
    Pa1 = Pa;
    Pa = C*Pa*C.t() + S22; 
    if (cform && !filtered) {
      for (i = 0; i < k; i++) {
        X(t, i) = a(i);
        for (j = 0; j < k; j++)
          PX(t*k+i, j) = Pa(i, j);
      }
    }
    a += K*v(t);
    Pa -= (C*Pa1*b + S12)*K.t(); 
    Pa = (Pa + Pa.t())/2;
    if (cform && filtered) {
      for (i = 0; i < k; i++) {
        X(t, i) = a(i);
        for (j = 0; j < k; j++)
          PX(t*k+i, j) = Pa(i, j);
      }
    }
  }
  
  if (xn) {
    x1 = a;
    P1 = Pa;
  }
  
  return true;
}

// [[Rcpp::export]]
bool ksC(const arma::colvec &y, const arma::colvec &b, const arma::mat &C,
         const arma::mat &S, const arma::colvec &x1, const arma::mat &P1,
         arma::mat &X, arma::mat &PX, bool cform) {
  
  int i, j, k, n, t;
  k = b.n_elem;
  n = y.n_elem;
  
  colvec v(n), s2v(n);
  colvec x1c = x1;
  mat P1c = P1;
  if (!kfC(y, b, C, S, x1c, P1c, v, s2v, X, PX, false, false, cform))
    return false;
  colvec r(k, fill::zeros), K(k), a(k), S12(k);
  mat N(k, k), L(k, k), Pa(k, k);
  for (i = 0; i < k; i++)
    S12(i) = S(i+1, 0);
  
  for (t = n-1; t > -1; t--) {
    for (i = 0; i < k; i++) {
      a(i) = X(t, i);
      for (j = 0; j < k; j++)
        Pa(i, j) = PX(t*k+i, j);
    }
    K = (C*Pa*b + S12)/s2v(t);
    L = C - K*b.t();
    r = b*v(t)/s2v(t) + L.t()*r;
    N = b*b.t()/s2v(t) + L.t()*N*L;
    a += Pa*r;
    Pa -= Pa*N*Pa.t();
    for (i = 0; i < k; i++) {
      X(t, i) = a(i);
      for (j = 0; j < k; j++)
        PX(t*k+i, j) = Pa(i, j);
    }
  }
  
  if (cform) {
    for (t = 0; t < n - 1; t++) {
      for (j = 0; j < k; j++) 
        X(t, j) = X(t+1, j);
    }
    for (j = 0; j < k; j++) 
      X(n-1, j) = x1c(j);
    for (t = 0; t < (n - 1)*k; t++) {
      for (j = 0; j < k; j++) 
        PX(t, j) = PX(t+k, j);
    }
    for (i = 0; i < k; i++) 
      for (j = 0; j < k; j++) 
        PX((n-1)*k+i, j) = P1c(i, j);
  }
  
  return true;
}
