#include "RcppArmadillo.h"
#include "pol.h"
#include "res.h"
#include "stsm.h"
using namespace arma;

// [[Rcpp::export]]
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
      dot += L(j, i)*L(j, i);
    
    L(j,0) = A(j,0) - dot;
    if (L(j, 0) <= 0.0)
      return false;
    L(j, 0) = std::sqrt(L(j, 0));
    
    for (i = 1; i < n; i++) {
      if(j+i < m){
        dot = 0;
        for (k = 1; k < n - i; k++)
          //if(i+k < n)
            dot += L(j+i, i+k)*L(j, k);
        L(j+i, i) = (A(j+i, i) - dot)/L(j, 0);
      }
    }
    detL *= L(j, 0);
  }
  
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
    //for (k=i-1; k>-1; k--)
    //  if(i-k < n)
    for (k = std::max(0, i - (n - 1)); k < i; k++)     
        sum -= L(i, i-k)*x(k);
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
  colvec g2(k+1, fill::zeros);
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
  if (!bandchol(G, L, detL)) {
    return datum::nan;
  }
  bandcholsol(L, e);
  s2(0) = dot(e, e)/(g0*n);
  ll = -0.5*n*(1 + log(2*datum::pi) + log(s2(0)) + 2*log(detL)/n + log(g0) );
  return ll;
 
}

// [[Rcpp::export]]
double llucaC(arma::colvec &w, const arma::colvec &phi, const arma::mat &A,
              const arma::mat &S, arma::colvec &s2, int res) {
  
  int h, i, j, k, p, q, r, n;
  double sum, detL, ll, g0;
  n = w.n_elem;
  q = A.n_cols - 1;
  p = phi.n_elem - 1;
  r = std::max(p, q + 1);  
  k = A.n_rows;
  colvec psi, g1;
  colvec g2(q+1, fill::zeros);
  mat Phi0, Theta0, Psi0, G0, Psi;
  mat G(n, r, fill::zeros);
  mat L(n, r, fill::zeros);
  colvec e = w;

  // Initial conditions
  if (p > 0) {
    e = condresC(w, phi, colvec(1, fill::ones));
    // Psi weights
    Psi = A;
    for (i = 1; i < q; i++ ) {
      for (j = 1; j <= p && i-j>=0; j++) {
        Psi.col(i) -= phi(j)*Psi.col(i-j);
      }
    }
    
    // Phi0*G0*Phi0'
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
      for (j = i; j <= q; j++)
        sum += dot(A.col(j).t()*S, Psi.col(j-i));
      g1(i) += sum;
    }
    g1 = solve(G0, g1);
    g1.shed_row(p);
    G0 = Phi0*toeplitz(g1)*Phi0.t();
    for (i = 0; i < p; i++)
      for (j = 0; j <= i; j++)
        G(i, i - j) = G0(i, j);
    
    // cov(w0, a0)
    Theta0.zeros(q, q*k);
    for (i = 0; i < q; i++) {
      for (j = i; j < q; j++)
        for (h = 0; h < k; h++)
          Theta0(i, (j-i)*k+h) = A(h, j+1);
    }
    
    Psi0.zeros(p, q*k);
    for (i = 0; i < p; i++) {
      for (j = i; j < q; j++)
        for (h = 0; h < k; h++)
          Psi0(i, j*k+h) = Psi(h, j-i);
    }
    Phi0 = Phi0*Psi0*kron(eye(q, q), S)*Theta0.t();
    Psi0.zeros(r, r);
    for (i = 0; i < p; i++) {
      for (j = 0; j < q; j++) {
        Psi0(i, j) += Phi0(i, j);
        Psi0(j, i) += Phi0(i, j);
      }
    }
    for (i = 0; i < r; i++)
      for (j = 0; j <= i; j++)
        G(i, i - j) += Psi0(i, j);
  }

  // Autocovariances MA part  
  for (i = 0; i <= q; i++) {
    sum = 0;
    for (j = i; j <= q; j++)
      sum += dot(A.col(j).t()*S, A.col(j-i));
    g2(i) += sum;
  }
  for (j = 0; j <= q; j++)
    for (i = j; i < n; i++)
      G(i, j) += g2(j);
  g0 = G.col(0).max();
  if (g0 > 1) G /= g0;
  else g0 = 1;
  if (!bandchol(G, L, detL))  {
    return datum::nan;
  }
  bandcholsol(L, e);
  s2(0) = dot(e, e)/n;
  if (detL > 0) {
    ll = -0.5*n*(1.0 + log(2.0*datum::pi) + log(s2(0)) + 2.0*log(detL)/n);
    if (res == 1) w = e*pow(detL, 1.0/n)/sqrt(g0);
    else if (res == 2) w = e/sqrt(g0);
  } else { 
    ll = -0.5*n*(1.0 + log(2.0*datum::pi) + log(s2(0)) + 2.0*mean(log(L.col(0))));  
    if (res == 1) w = e*prod(pow(L.col(0), 1.0/n))/sqrt(g0);
    else if (res == 2) w = e/sqrt(g0);
  }
  s2(0) /= g0;
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
         const arma::mat &S, const arma::colvec &x0, const arma::mat &P0,
         arma::colvec &v, arma::colvec &s2) {
  // STS model must be in lagged form
  int i, j, k, n, t;
  double sum, s11;
  k = b.n_elem;
  n = y.n_elem;
  if (v.n_elem < n) return false;
  colvec K(k), S12 = S.submat(1, 0, k, 0), a(k);
  mat S22 = S.submat(1, 1, k, k), Pa(k, k);
  s11 = S(0,0);
  a = x0; Pa = P0;
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
         const arma::mat &S, const arma::colvec &x0, const arma::mat &P0,
         arma::colvec &v, arma::colvec &s2, arma::mat &X, arma::mat &PX, 
         bool filtered) {
  // SS model must be in lagged form: 
  // y(t) = b'x(t-1) + u(t); x(t) = Cx(t-1) + v(t)
  
  int i, j, k, n, t;
  double sum, s11;
  k = b.n_elem;
  n = y.n_elem;
  if (v.n_elem < n || b.n_elem != k || C.n_rows != C.n_cols) return false;
  if (S.n_rows != S.n_cols || S.n_rows != k + 1) return false;
  if (x0.n_elem != k || P0.n_rows != k || P0.n_rows != P0.n_cols) return false;
  
  colvec K(k), S12 = S.submat(1, 0, k, 0), x = x0, x1(k);
  mat S22 = S.submat(1, 1, k, k), P = P0, P1(k,k), L(k,k);
  s11 = S(0,0);
  for (t = 0; t < n; t++) {
    v(t) = y(t) - dot(b, x);
    s2(t) = dot(b.t()*P, b) + s11; 
    if (std::abs(s2(t)) < 1e-10) return false;    
    K = (C*P*b + S12)/s2(t);
    x1 = C*x;
    P1 = C*P*C.t() + S22;
    if (!filtered) { // forecasting
      X.row(t) = x1.t();
      PX.rows(t*k, (t+1)*k-1) = P1;
    }
    x = x1 + K*v(t);
    P = P1 - K*K.t()*s2(t);
    P = (P + P.t())/2;
    if (filtered) {
      X.row(t) = x.t();
      PX.rows(t*k, (t+1)*k-1) = P;
    }
  }
  
  return true;
}

// [[Rcpp::export]]
bool ksC(const arma::colvec &y, const arma::colvec &b, const arma::mat &C,
         const arma::mat &S, arma::colvec &x0, arma::mat &P0,
         arma::mat &X, arma::mat &PX) {
  // SS model must be in lagged form: 
  // y(t) = b'x(t-1) + u(t); x(t) = Cx(t-1) + v(t)
  int i, j, k, n, t;
  k = b.n_elem;
  n = y.n_elem;
  
  colvec v(n), s2v(n);
  if (!kfC(y, b, C, S, x0, P0, v, s2v, X, PX, true))
    return false;
  colvec r(k, fill::zeros), K(k), a(k), S12 = S.submat(1, 0, k, 0);
  mat N(k, k, fill::zeros), L(k, k), Pa(k, k);

  for (t = n - 1; t >= 0; t--) {
    Pa = PX.rows(t * k, (t + 1) * k - 1);
    X.row(t) += r.t()*Pa;
    PX.rows(t * k, (t + 1) * k - 1) -= Pa*N*Pa.t();
    K = (C*Pa*b + S12)/s2v(t);
    L = C - K*b.t();
    r = b*v(t)/s2v(t) + L.t()*r;
    N = b*b.t()/s2v(t) + L.t()*N*L;
  }
  x0 += Pa*r;
  P0 -= P0*N*P0.t();
  return true;
}
