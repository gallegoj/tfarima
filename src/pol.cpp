#include "RcppArmadillo.h"
#include "pol.h"
using namespace arma;

// Polynomial roots
//
// \code{polyrootsC} C function to compute the roots of a lag polynomial.
//
// @param pol Numeric vector, c(1, coef_1, ..., coef_p).
// 
// @return \code{polyrootsC} returns a matrix with five columns showing 
// the real and imaginary parts and the modulus, the frequency and the period
// of each root.  
// 
// @section Warning:
// This C function is mainly used to restrict the parameters into the 
// admissible region during the estimation of ARIMA models. 
// 
// [[Rcpp::export]]
arma::mat polyrootsC(const arma::colvec &pol) {
  cx_vec r = roots(pol);
  return sortrootsC(r);
}

bool simeqC(double x, double y, double tol) {
  if (fabs(x - y) < tol) return true;
  else return false;
}

bool ltC(double x, double y, double tol) {
  if (std::fabs(x - y) < tol) return false; 
  return (x < y);
}

// [[Rcpp::export]]
arma::mat sortrootsC(const arma::cx_colvec &r) {
  int h, i, j, p;
  double d;
  std::complex <double> cx;
  
  p = r.n_elem;
  Col<int> indx(p);
  indx.fill(1);
  mat T(p, 6);
  
  for (j = 0; j < p; j++) {
    cx = r(j);
    T(j, 0) = cx.real();  T(j, 1) = cx.imag();
    T(j, 2) = abs(cx);
    if (simeqC(T(j, 2), 0.0)) T(j, 3) = 0.0;
    else T(j, 3) = acos(T(j, 0)/T(j, 2)) / (2.0 * datum::pi);    
    if ( simeqC(pow(T(j, 1), 2), 0) ) {
      T(j, 1) = 0;
      if (T(j, 0) > 0) T(j, 3) = 0;
      else T(j, 3) = 0.5;
    }
    T(j, 4) = (simeqC(T(j, 3), 1e-6)) ? datum::inf : 1.0/T(j, 3);    
    T(j, 5) = 1;
  }
  
  // Sort by frequency and modulus
  for (i = 0; i <  p; i++) {
    for (j = i + 1; j <  p; j++) {
      if ( ltC(T(j, 2), T(i, 2)) || 
           ( T(j, 3) < T(i, 3) && simeqC(T(j, 2), T(i, 2)) ) ) {
        for (h = 0; h < 6; h++) {
          d = T(i, h);
          T(i, h) = T(j, h);
          T(j, h) = d;
        }
      }
    }
  }
  
  // Multiplicity
  for (i = 0; i <  p; i++) {
    if (indx(i)) {
      for (j = i + 1; j <  p; j++) {
        if ( simeqC(T(j, 3),  T(i, 3)), 1e-6)  {
          if ( simeqC(T(j, 0), T(i, 0), 1e-6)  &&
               simeqC(T(j, 1), T(i, 1), 1e-6) ) {
              T(i, 5) += 1;
              indx(j) = 0;
          }
        } else {
          break;
        }
      }
    }
  }
  
  T = T.rows(find(indx > 0));  
  return T;
  
}

// [[Rcpp::export]]
arma::mat roots2polC(arma::mat A, bool check) {
  int m, i, j, r;
  arma::vec pol(1, arma::fill::ones);
  arma::vec pol1(2, arma::fill::zeros);
  arma::vec pol2(3, arma::fill::zeros);
  pol(0) = pol1(0) = pol2(0) = 1;
  r = A.n_rows - 1;
  for (i = 0; i <  r; i++) {
    if (check) {
      if(A(i, 2) < 1) {
        A(i, 2) = 1/A(i, 2);
        A(i, 0) = 1/A(i, 0);
      }
     }
    if ( simeqC(A(i, 2), A(i+1, 2)) && simeqC(A(i, 1), -A(i+1, 1)) ) {
      pol2(1) = -2*cos(2*datum::pi*A(i, 3))/A(i, 2);
      pol2(2) = pow(1/A(i, 2), 2);
      m = int(A(i, 5)); 
      for (j = 0; j < m; j++)
        pol = polymultC(pol, pol2);
      ++i;
    } else {
      pol1(1) = -1/A(i, 0);
      m = int(A(i, 5)); 
      for (j = 0; j < m; j++)
        pol = polymultC(pol, pol1);
    } 
  }
  
  if (i == r) {
    if (check) {
      if(A(i, 2) < 1) {
        A(i, 2) = 1/A(i, 2);
        A(i, 0) = 1/A(i, 0);
      }
    }
    pol1(1) = -1/A(i, 0);
    for (j = 0; j < int(A(i, 5)); j++)
      pol = polymultC(pol, pol1);
  }
  return pol;
  
}


// Check parametric admissibility
//
// \code{admregC} C function to check if the roots of a lag 
// polynomial lie inside the admissibe region.
//
// @param pol Numeric vector, c(1, coef_1, ..., coef_p).
// @param ar logical. If TRUE, roots must lie outside the unit circle;
//       if FALSE, roots can also lie on the unit circle. 
// 
// @return \code{admregC} returns TRUE or FALSE.  
// 
// [[Rcpp::export]]
bool admregC(const arma::colvec &pol, bool ar) {
  
  int j;
  int p = pol.n_elem-1;
  
  mat A(p, p, fill::zeros);
  cx_vec eigval;
  
  for (j = 0; j < p; j++)
    A(0, j) = -pol(j+1);
  for (j = 1; j < p; j++)
    A(j, j-1) = 1.0;
  
  eig_gen(eigval, A);
  if (ar) {
    for (j = 0; j < p; j++) {
      if (abs(eigval(j))>= 1.0)
        return false;
    }
  } else {
    for (j = 0; j < p; j++) {
      if (abs(eigval(j))>1.0)
        return false;
    }
  }
  
  return true;
  
}

// Polynomial multiplication
//
// \code{polymultC} Computes the product of two lag polynomials 
//  c(B) = a(B)b(B).
//
// @param pol1,pol2 are numeric vectors with the coefficients of 
// two polynomials.
// 
// @return \code{polymultC} returns a column vector with the coefficients
// of the product polynomial. 
// 
// @section Warning:
// This C function is mainly used to unscramble the AR, I amd MA operators. 
// 
// [[Rcpp::export]]
arma::colvec polymultC(const arma::colvec &pol1, const arma::colvec &pol2) {
  int r = pol1.n_elem;
  int s = pol2.n_elem;
  arma::colvec pol(r+s-1, fill::zeros);
  
  for (int i = 0; i < r; i++) {
    for (int j = 0; j < s; j++) {
      pol(i+j)  += pol1(i)*pol2(j);
    }
  }
  return pol;
}

// Polynomial division
//
// \code{polydivC} computes the quotient of two lag polynomials c(B) =
// a(B)/b(B).
//
// @param pol1,pol2 are numeric vectors with the coefficients of two
//   polynomials.
// @param rem logical. If true, the function returns the remainder of the
//   division return; if false, the quotient is returned.
// @ param tol tolerance to check if a value is null.
//
// @return \code{polydivC} returns a column vector with the coefficients of the
//   quotient or remainder polynomial.
//
//   [[Rcpp::export]]
arma::colvec polydivC(const arma::colvec &pol1, const arma::colvec &pol2, 
                      bool rem, double tol) {
  
  int i, j, l1, l2, l3;
  double d;
  l1 = pol1.n_elem - 1;
  l2 = pol2.n_elem - 1;
  if (l1 < l2) {
    if (rem) return pol1;
    else return colvec(1, fill::zeros);
  }
  l3 = l1 - l2;
  
  colvec p1 = pol1;
  colvec p2 = pol2;
  colvec q(l3+1, fill::zeros);
  
  for (i = 0; i <= l3; i++) {
    d = p1(l1-i)/p2(l2);
    q(l3-i) = d;
    for (j = 0; j <= l2; j++)
      p1(l1-i-j) -= d * p2(l2 - j);
  }
  
  if (!rem) {  
    for (j = l3; j > 0; j--) { 
      if(simeqC(q(j), 0, tol)) q.shed_row(j);
      else break;
    }
    return q;  
  } else {
      for (j = l1; j > 0; j--) { 
        if(simeqC(p1(j), 0, tol)) p1.shed_row(j);
        else break;
      }
      return p1;  
  }
}

// [[Rcpp::export]]
arma::colvec polygcdC(const arma::colvec &pol1, const arma::colvec &pol2,
                      double tol) {
  
  if (pol2.n_elem > pol1.n_elem)
    return polygcdC(pol2, pol1, tol);
  colvec p1 = pol1;
  colvec p2 = pol2;
  colvec p3;
  while(!(simeqC(p2(0), 0, tol) && p2.n_elem == 1) ) { //    
    p3 = polydivC(p1, p2, true, tol);
    p1 = p2;
    p2 = p3;
  }
  return p1/p1(0);  
}

// [[Rcpp::export]]
arma::colvec polyprsC(const arma::colvec &pol1, const arma::colvec &pol2) {
  
  if (pol2.n_elem > pol1.n_elem) return polyprsC(pol2, pol1);
  int h, l1, l2;
  double d, b;
  colvec p1 = pol1;
  colvec p2 = pol2;
  colvec p3;
  d = 0;
  while(!simeqC(p2(0), 0, 1.490116e-08)) {
    l1 = p1.n_elem - 1; l2 = p2.n_elem - 1; 
    h = l1 - l2;
    b = pow(-1, d+1)*p1(l1)*pow(h, d);
    h = h*pow(p2(l2)/h, d);
    p3 = polydivC(p1, p2, true);
    p1 = p2;
    p2 = p3/b;
  }
  if(p1.n_elem == 1) 
    return colvec(1, fill::ones);
  else {
    return p1;  
  }
}


// Raise a polynomial to some power
//
// \code{polyraiseC} expands a polynomial to the d-th power.
//
// @param pol is a numeric vector and d is a positive integer.
// @param d polynomial is raised to the power d.
// 
// @return \code{polyraiseC} returns a numeric vector with the coefficients
// of the expanded polynomial. 
// 
// [[Rcpp::export]]
arma::colvec polyraiseC(const arma::colvec &pol, int d) {

  if (d == 1) {
    return pol;
  } else if (d == 2) {
    return polymultC(pol, pol);
  } else {
    vec pol1 = polymultC(pol, pol);
    for (int i = 2; i < d; i++) {
      pol1 = polymultC(pol1, pol);
    }
    return pol1;
  }
}

// Rational polynomial.
//
// \code{polyratioC} computes a rational polynomial of degree d from 
// the ratio of two lag polynomials c(B) = a(B)/b(B).
//
// @param num Numerator polynomial, c(1, a_1, ..., a_p).
// @param den Denominator polynomial, c(1, b_1, ..., b_q).
// @param d Degree of the rational polynomial, integer.
// 
// @return \code{polyratioC} returns a numeric vector of dimension "d+1" with the coefficients
// of the rational polynomial. 
// 
// [[Rcpp::export]]
arma::colvec polyratioC(const arma::colvec &num, const arma::colvec &den, int d) {
  int i, j;
  int p = num.n_elem-1;
  int q = den.n_elem-1;
  double x;
  
  if (d < 0) d = p + q;
  
  vec pol(d+1, fill::zeros);
  
  for (j = 0; j <= d; j++) {
    x = 0;
    for (i = 1; i <= q; i++){
      if (j-i>-1) {
        x  += den(i)*pol(j-i);
      }				
    }
    if (j <= p) {
      pol(j) = num(j)-x;				
    }
    else {
        pol(j) = -x;			
    }
  }
  
  return pol;
}


