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
