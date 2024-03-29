% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{wold.pol}
\alias{wold.pol}
\title{Wold polynomial}
\usage{
wold.pol(a, type = c("wold", "palindromic", "cramer-wold"), tol = 1e-05)
}
\arguments{
\item{type}{character indicating the type of polynomial: (1) Wold polynomial, 
(2) Palindromic polynomial and (3) Cramer-Wold factor.}

\item{tol}{tolerance to check if an autocovariance is zero.}

\item{x}{numeric vector, coefficients of a palindromic or a Wold polynomial.}
}
\value{
Numeric vector.
}
\description{
Transforming a palindromic polymonial into a Wold polynomial/ 
Computing the Cramer-Wold factorization
}
\details{
\code{wold.pol} can be used with three purposes:

(1) to transform a self-reciprocal or palindromic polynomial
 a_0 + a_1(B+F) + a_2(B^2+F^2) + ... + a_p(B^p+F^p) 
 into a Wold polynomial
 b_0 + b_1(B+F) + b_2(B+F)^2 + ... + b_p(B+F)^p;
 
 (2) to revert the previous transformation to obtain the palindromic
 polynominal from a Wold polynomial and
 
 (3) to compute the Cramer-Wold factorization: b(B+F) = c(B)c(F).
}
\examples{
wold.pol(c(6, -4, 1))
}
