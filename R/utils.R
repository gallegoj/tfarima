#' MA process compatible with a vector of autocovariances
#'
#' \code{autocov2MA} computes the error variance and the MA polynomial from a
#' vector of autocovriances.
#'
#' @param x a vector of q autocovariances.
#'
#' @return A vector of the form c(sig2, 1, -theta1, ..., -thetaq).
#' 
#' @examples 
#' ma1 <- um(ma = "1 - 0.8B", sig2 = 0.5)
#' autocov2MA(autocov(ma1, 1))
#' @export
autocov2MA <- function(x) {
  x <- as.numeric(x)
  n <- length(x)
  while(x[n] == 0 && n > 0) {
    x <- x[-n]
    n <- n - 1
  }
  if (n == 0) return(0)
  stopifnot(x[1] > 0)
  if (n == 1)
    return(x)
  b <- rep(0, n)
  b[1:2] <- x[1:2]
  if (n > 2) {
    p0 <- rep(0, n)
    p1 <- rep(0, n)
    p2 <- rep(0, n)
    p0[1] <- 2
    p1[2] <- 1
    for (i in 3:n) {
      p2[1:i] <- c(0, p1[1:(i-1)]) - p0[1:i]
      b[1:i] <- b[1:i] + x[i]*p2[1:i]
      p0 <- p1
      p1 <- p2
    }
  }
  
  r <- polyroot(b)
  l <- length(r)
  r <- sapply(r, function(x) {
    c((x + sqrt(x^2 - 4))/2, (x - sqrt(x^2 - 4))/2)
  })
  r <- as.vector(r)
  r <- r[order(Mod(r))]
  r <- r[1:l]
  ma <- 1
  for (r1 in r)
    ma <- c(ma, 0) - c(0, ma*r1)
  ma <- Re(ma)
  ma <- c(b[n]/ma[length(ma)], ma)
  if (!is.finite(ma[1]))
    autocov2MA(x[-n])
  else return(ma)
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
