#' MA process compatible with a vector of autocovariances
#'
#' \code{autocov2MA} computes the error variance and the MA polynomial from a
#' vector of autocovriances.
#'
#' @param x a numeric vector with q autocovariances.
#' @param method character. Two methods are provided to compute the MA
#'   parameters: roots, based on the roots of the autocovariance generating
#'   function; acov, based on the non-linear system of equations relating
#'   autocovariances and MA parameters.
#' @param tol tolerance to check if an autocovariance is null.
#'
#' @return A vector of the form c(s2, 1, -theta1, ..., -thetaq).
#'
#' @examples
#' ma1 <- um(ma = "1 - 0.8B", sig2 = 0.5)
#' autocov2MA(autocov(ma1, 1))
#' autocov2MA(autocov(ma1, 1), "acov")
#' @export
autocov2MA <- function(x, method = c("roots", "acov"), tol = 1e-5) {
  method <- match.arg(method, c("roots", "acov"))
  if (method == "acov") {
    code <- 0
    ma <- as.numeric(acovtomaC(x, code))
    ma <- c(s2 = ma[1]^2, ma= c(1, -ma[-1]))
    names(ma)[-1] <- paste0("ma", 0:(length(ma)-2))
    return(ma)
  }
  
  x <- as.numeric(x)
  n <- length(x)
  while(abs(x[n]/x[1]) < tol && n > 0) {
    x <- x[-n]
    n <- n - 1
  }
  if (n == 0) return(0)
  stopifnot(x[1] > 0)
  if (n == 1)
    return(x)
  x <- wold.pol(x, type = "c")
  names(x)[-1] <- paste0("ma", 0:(length(x)-2))
  return(x)
}

#' Wold polynomial
#' 
#' Transforming a palindromic polymonial into a Wold polynomial/ 
#' Computing the Cramer-Wold factorization
#'
#' \code{wold.pol} can be used with three purposes:
#' 
#' (1) to transform a self-reciprocal or palindromic polynomial
#'  a_0 + a_1(B+F) + a_2(B^2+F^2) + ... + a_p(B^p+F^p) 
#'  into a Wold polynomial
#'  b_0 + b_1(B+F) + b_2(B+F)^2 + ... + b_p(B+F)^p;
#'  
#'  (2) to revert the previous transformation to obtain the palindromic
#'  polynominal from a Wold polynomial and
#'  
#'  (3) to compute the Cramer-Wold factorization: b(B+F) = c(B)c(F). 
#'
#' @param x numeric vector, coefficients of a palindromic or a Wold polynomial.
#' @param type character indicating the type of polynomial: (1) Wold polynomial, 
#' (2) Palindromic polynomial and (3) Cramer-Wold factor.
#' @param tol tolerance to check if an autocovariance is zero.
#' @return Numeric vector.
#'
#' @examples
#' wold.pol(c(6, -4, 1))

#' @export
wold.pol <- function(x, type = c("wold", "palindromic", "cramer-wold"), 
                     tol = 1e-5) {
  type <- tolower(type)
  type <- match.arg(type)
  if (startsWith("cramer-wold", type)) {
    if (x[1] == 0) 
      return(0)
    l <- length(x)
    while(abs(x[l]/x[1]) < tol) {
      x <- x[-l]
      l <- l - 1
    }
    if (l == 1)
      return(x)
    x <- wold.pol(x)
    r <- polyroot(x)
    n <- length(r)
    r <- sapply(r, function(x) {
      d <- sqrt(x^2 - 4)
      r1 <- (x + d)/2
      r2 <- (x - d)/2
      if (Mod(r1) >= Mod(r2)) return(c(r1, r2))
      else return(c(r2, r1))
    })
    r <- t(r)
    print(cbind(r, Mod(r), Arg(r) ))
    if (n > 1) {
      i <- (abs(Mod(r[, 1]) - 1) < tol) & ( (abs(Re(r[, 1])) < tol) | (abs(Im(r[, 1])) > tol) )
      print(i)
      n1 <- sum(i)
      print(n1)
      if (n1 > 1 && (n1 %% 2) == 0) {
        r1 <- r[i, ]
        r1 <- r1[order(abs(Arg(r1[, 1]))), ]
        print(r1)
        for (j in seq(1, n1, 2)) {
          if (sign(Arg(r1[j, 1])) == sign(Arg(r1[j+1, 1]))) 
            r2 <- (r1[j, ] + r1[j+1, ])/2
          else
            r2 <- (r1[j, ] + r1[j+1,2:1])/2
          r1[j, ] <- r2
          r1[j+1, ] <- c(r2[2], r2[1])
        }
        print(r1)
        r[i, ] <- r1
      }
    }
    r <- r[, 1]
    ma <- 1
    for (r1 in r)
      ma <- c(ma, 0) - c(0, ma/r1)
    ma <- Re(ma)
    ma <- c(x[l]/ma[length(ma)], ma)
    if (!is.finite(ma[1]))
      wold.pol(x[-l], "c")
    else return(ma)
  }
  x <- as.numeric(x)
  n <- length(x)
  if (n < 3) return(x)
  if (startsWith("palindromic", type)) {
    b <- rep(0, n)
    b[1:2] <- x[1:2]
    for (i in 3:n) {
      for (j in 0:(i-1)) {
        if (i-2*j<1) break
        b[i-2*j] <- b[i-2*j] + x[i]*choose(i-1,j)
      }
    }
    return(b)
  }
  p0 <- rep(0, n)
  p1 <- rep(0, n)
  p2 <- rep(0, n)
  b <- rep(0, n)
  b[1:2] <- x[1:2]
  p0[1] <- 2
  p1[2] <- 1
  for (i in 3:n) {
    p2[1:i] <- c(0, p1[1:(i-1)]) - p0[1:i]
    b[1:i] <- b[1:i] + x[i]*p2[1:i]
    p0 <- p1
    p1 <- p2
  }
  return(b)
}

#' Standardize time series
#'
#' \code{std} standardizes a time series.
#' 
#' @param x a \code{ts} object.
#'
#' @return The standardized time series.
#' 
#' @export
std <- function(x) {
  (x - mean(x))/sd(x)
}

#' Outlier dates
#' 
#' \code{outlierDates} shows the indeces and dates of outliers.
#' 
#' @param x an \code{ts} object.
#' @param c critical value to determine whether or not an observation is an outlier.     
#'
#' @return 
#' 
#' A table with the indices, dates and z-scores of the outliers.
#'
#' @export
outlierDates <- function(x, c = 3) {
  if (!is.ts(x)) x <- as.ts(x)
  freq <- frequency(x)
  indx <- seq(x)
  x <- (x - mean(x))/sd(x)
  if (freq == 12) 
    d <- cbind(indx = indx, year = floor(time(x)), month = cycle(x))
  else if (freq == 4) 
    d <- cbind(indx = indx, year = floor(time(x)), quarter = cycle(x))
  else if (freq > 1) 
    d <- cbind(indx = indx, year = floor(time(x)), season = cycle(x))
  else if (start(x)[1] > 1)
    d <- cbind(indx = indx, year = floor(time(x)))
  else
    d <- cbind(indx = indx)
  d[abs(x) > c, ]
}

#' Value of a time series at a date
#' 
#' \code{tsvalue} select a value from a time series by date.
#' 
#' @param x an \code{ts} object.
#' @param date the time of the specific observation, c(year, month/quarter).     
#'
#' @return 
#' 
#' The value of the observation, double.
#'
#' @export
tsvalue <- function(x, date) {
  stopifnot(is.ts(x))
  start <- start(x)
  s <- frequency(x)
  n <- length(x)
  if (s > 1) t <- (date[1] - start[1] + 1)*s - (start[2] - 1) - (s - date[2])
  else t <- (date[1] - start[1] + 1)*s
  stopifnot(t > 0 || t <= n)
  x[n]
}

#' Annual sum
#'
#' \code{S} generates the annual sum of a monthly or quarterly time series.
#'
#' @param x an \code{ts} object.
#' @param extend logical. If TRUE, the transformed series is extended with NA's
#'   to have the same length as the origianl series.
#' @return
#'
#' The transformed time series, a \code{ts} object.
#'
#' @export
S <- function(x, extend = TRUE) {
  if (is.mts(x)) {
    y <- sapply(1:ncol(x), function(k) {
      S(x[, k], extend)
    })
    return(y)
  }
  stopifnot(is.ts(x))
  start <- start(x)
  s <- frequency(x)
  n <- length(x)
  y <- sapply(s:n, function(k) sum(x[(k-s+1):k])) 
  if (extend) y <- c(rep(NA, s), y)
  ts(y, end = end(x), frequency = s)
}

YearSeason <- function(x, year.digits = 4) {
  if (!is.ts(x)) x <- as.ts(x)
  freq <- frequency(x)
  n <- length(x)
  y <- as.character(floor(time(x)))
  s <- as.character(cycle(x))
  if (year.digits == 2 && nchar(y[1]) == 4 && nchar(y[n]) == 4)
    y <- substr(y, 3, 4)
  s[nchar(s) == 1] <- paste(0, s[nchar(s)==1], sep = "")
  if ( freq == 4 ||freq == 12) 
    return(paste(y, s, sep = "."))
  else
    return(y)
}

yearSsum <- function(x, bc = FALSE) {
  if (!is.ts(x)) x <- as.ts(x)
  s <- frequency(x)
  ts(diffC(x, rep(1, s), bc), end = end(x), frequency = s)  
} 
