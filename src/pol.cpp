#include "RcppArmadillo.h"
#include "pol.h"
using namespace arma;

// [[Rcpp::export]]
double polyevalC(const arma::colvec &pol, double z) {
  int i, p;
  double d = 0;
  p = pol.n_elem-1;
  for (i=0; i<p; i++) d = d*z + pol(i);
  return d;
}

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
  if (fabs(x - y) < tol) return false;
  else if (x < y) return true;
  else return false;
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
    T(j, 3) = acos( T(j, 0)/T(j, 2) )/(2.0*datum::pi);
    if ( simeqC(T(j, 3), 0) ) T(j, 3) = 0; 
    T(j, 4) = 1.0/T(j, 3);
    T(j, 5) = 1;
  }
  
  // Sort by frequency and modulus
  for (i = 0; i <  p; i++) {
    for (j = i + 1; j <  p; j++) {
      if ( ltC(T(j, 2), T(i, 2)) || ( T(j, 3) < T(i, 3) && simeqC(T(j, 2), T(i, 2)) ) ) {
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
        if ( simeqC(T(j, 3),  T(i, 3))  )  {
          if ( simeqC(T(j, 0), T(i, 0))  && simeqC(T(j, 1), T(i, 1)) ) {
              T(i, 5) += 1;
              indx(j) = 0;
          }
        } else {
          break;
        }
      }
    }
  }
  
  for (i = p-1; i >  -1; i--) {
    if (indx(i) == 0) T.shed_row(i);
  }
  
  return T;
  
}

// [[Rcpp::export]]
arma::mat combinerootsC(arma::mat T) {
  int h, i, j, p;
  double d;

  p = T.n_rows; 
  Col<int> indx(p);
  indx.fill(1);

  // Sort by frequency and modulus
  for (i = 0; i <  p; i++) {
    for (j = i + 1; j <  p; j++) {
      if ( ltC(T(j, 2), T(i, 2)) || (ltC(T(j, 3), T(i, 3)) && simeqC(T(j, 2), T(i, 2)) ) ) {
        for (h = 0; h < 6; h++) {
          d = T(i, h);
          T(i, h) = T(j, h);
          T(j, h) = d;
        }
      }
    }
  }

  // Greatest Multiplicity
  for (i = 0; i <  p; i++) {
    if (indx(i)) {
      for (j = i + 1; j <  p; j++) {
        if ( simeqC(T(j, 3),  T(i, 3))  )  {
          if ( simeqC(T(j, 0), T(i, 0))  && simeqC(T(j, 1), T(i, 1)) ) {
            if (T(i, 5) < T(j, 5) ) indx(i) = 0;
            else indx(j) = 0;
          }
        } else {
          break;
        }
      }
    }
  }
  
  for (i = p-1; i >  -1; i--) {
    if (indx(i) == 0) T.shed_row(i);
  }
  
  return T;
  
}

// [[Rcpp::export]]
arma::mat roots2polC(arma::mat T) {
  int h, i, j, p, r, s;
  double d;
  
  p = T.n_rows; 
  mat A(p, 5, fill::zeros);

  r = 0;    
  for (i = 0; i <  p; i++) {
    if (T(i, 5) > 0) {
      s = 1;
      h = 1;
      d = T(i, 5);
      for (j = i + 1; j <  p; j++) {
        if (T(j, 5) > 0 && simeqC(T(j, 2), T(i, 2)) ) {
          h++;
          if (T(j, 5) < d)  d = T(j, 5);
          if (T(j, 4) > s && std::isfinite(T(j, 4)))  s = T(j, 4);
        }
      }
        
      if (h == 1) {
        A(r, 0) = 1; A(r, 1) = d; A(r, 2) = 1; A(r, 3) = T(r, 2);
      } else if (h==2) {
        
      } else {
        if (simeqC(T(i, 3), 0) ) {
          A(r, 0) = d; A(r, 1) = d; A(r, 2) = 1; A(r, 3) = T(i, 2);
        } else {
          
        }
      }
      
      for (j = i; j <  p; j++)
        if (T(j, 5) > 0 && simeqC(T(j, 2), T(i, 2)) )
          T(j, 5) -= d;
        
      i--;
      r++;
    }
  }
  
  return A;
  
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
  
  int j, p;
  p = pol.n_elem-1;
  
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
  
  int i, j, r, s;

  r = pol1.n_elem - 1;
  s = pol2.n_elem - 1;
  
  vec pol(r+s+1, fill::zeros);

  for (i = 0; i <= r; i++) {
    for (j = 0; j <= s; j++) {
      pol(i+j)  += pol1(i)*pol2(j);
    }
  }
  
  return pol;
}


// [[Rcpp::export]]
arma::colvec polydivC(const arma::colvec &pol1, const arma::colvec &pol2, bool rem) {
  
  int i, j, l1, l2, l3;
  double d;
  l1 = pol1.n_elem - 1;
  l2 = pol2.n_elem - 1;
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
      if(q(j) == 0) q.shed_row(j);
      else break;
    }
    return q;  
  } else {
      for (j = l1; j > 0; j--) { 
        if(p1(j) == 0) p1.shed_row(j);
        else break;
      }
      return p1;  
  }
}

// [[Rcpp::export]]
arma::colvec polygcdC(const arma::colvec &pol1, const arma::colvec &pol2) {
  
  if (pol2.n_elem > pol1.n_elem) return polygcdC(pol2, pol1);
  colvec p1 = pol1;
  colvec p2 = pol2;
  colvec p3;
  while(!simeqC(p2(0), 0, 1.490116e-08) || p2.n_elem > 1) {
    p3 = polydivC(p1, p2, true);
    p1 = p2;
    p2 = p3;
  }
  if(p1.n_elem == 1) 
    return colvec(1, fill::ones);
  else {
    return p1; // /p1(0);  
  }
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

// Polynomial factors
//
// \code{polyfactorsC} C function called by the R function polyroots 
// to compute the roots of a polynomial in the lag operator.
//
// @param pol Numeric vector, c(1, coef_1, ..., coef_p).
// 
// @return \code{polyfactorsC} returns a matrix with five columns showing 
// the real and imaginary parts and the modulus, frequency and period
// of each root.  
// 
// @section Warning:
// Use the R function polyroots insted of polyrootsC.
// 
// [[Rcpp::export]]
arma::mat polyfactorsC(const arma::colvec &pol) {

  int h, i, r;
  mat A = polyrootsC(pol);
  r = A.n_rows;
  mat B(r, 4, fill::zeros);
  h = 0;
  r--;
  for (i = 0; i <  r; i++) {
    if ( (A(i, 2) == A(i+1, 2)) && (A(i, 1) == -A(i+1, 1)) ) {
      B(h, 0) = 2;
      B(h, 1) = 1;
      B(h, 2) = -2*A(i, 0);
      B(h, 3) = pow(A(i, 2), 2);
      i++;
    } else {
      B(h, 0) = 1;
      B(h, 1) = 1;
      B(h, 2) = -A(i, 0);
    } 
    h++;
  }
  
  if (i == r) {
    B(h, 0) = 1;
    B(h, 1) = 1;
    B(h, 2) = -A(i, 0);
    h++;
  }
  
  B.resize(h, 4);
  
  return B;
  
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


