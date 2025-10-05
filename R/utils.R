## tfarima/R/utils.R
## Jose L Gallego (UC)

#' Convert autocovariances to MA parameters
#'
#' Computes MA polynomial coefficients and error variance from autocovariances.
#'
#' @param x Numeric vector of autocovariances (length q).
#' @param method Estimation method: "roots" (Godolphin 1976) or "acov" (Wilson
#'   1969). Default is "roots".
#' @param tol Tolerance for zero autocovariance. Default 1e-5.
#'
#' @return Named vector: c(s2, ma0=1, ma1=-theta1, ..., maq=-thetaq).
#'
#' @references Godolphin, E. J. (1976). On the Cramér-Wold factorization.
#' \emph{Biometrika}, \strong{63}(2), 367-372. \doi{10.2307/2335982}
#'
#' Tunnicliffe Wilson, G. (1969). Factorization of the covariance generating
#' function of a pure moving average process. \emph{SIAM Journal on Numerical
#' Analysis}, \strong{6}(1), 1-7. \doi{10.1137/0706001}
#'
#' @examples
#' ma1 <- um(ma = "1 - 0.8B", sig2 = 0.5)
#' autocov2MA(autocov(ma1, 1))
#' autocov2MA(autocov(ma1, 1), method = "acov")
#' @export
autocov2MA <- function(x, method = c("roots", "acov"), tol = 1e-5) {
  method <- match.arg(method, c("roots", "acov"))
  x <- as.numeric(x)
  
  if (method == "acov") {
    code <- 0
    ma <- as.numeric(acovtomaC(x, code, tol, 100))
    ma <- c(s2 = ma[1]^2, ma= c(1, -ma[-1]))
    names(ma)[-1] <- paste0("ma", 0:(length(ma)-2))
    return(ma)
  }
  
  n <- length(x)
  while(n > 0 && abs(x[n]/x[1]) < tol) {
    x <- x[-n]
    n <- n - 1
  }
  if (n == 0) return(0)
  if (x[1] <= 0) stop("First autocovariance must be positive")  
  if (n == 1) return(x)
  res <- wold.pol(x, type = "c")
  names(res)[-1] <- paste0("ma", 0:(length(res)-2))
  return(res)
}

#' Standardize time series
#'
#' Centers and scales a time series to zero mean and unit variance.
#'
#' @param x A \code{ts} object or numeric vector.
#'
#' @return Standardized time series with same class as input.
#'
#' @export
std <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
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
  startx <- start(x)
  freq <- frequency(x)
  indx <- seq(x)
  xstar <- std(x)
  indx <- which(xstar > c)
  if (length(indx) == 0) return(NULL)
  if (startx[1] == 1 && freq == 1)
    res <- cbind(indx = indx)
  else
    res <- cbind(indx = indx, year = floor(time(x)[indx]))
    
  if (freq == 12) 
    res <- cbind(res, month = cycle(x)[indx])
  else if (freq == 4) 
    res <- cbind(res, quarter = cycle(x)[indx])
  else if (freq > 1) 
    res <- cbind(res, season = cycle(x)[indx])
  cbind(res, std_res = xstar[indx])
}

#' Extract time series value by date
#'
#' Retrieves the value of a time series at a specific date.
#'
#' @param x A \code{ts} object.
#' @param date Numeric vector: c(year) for annual data or
#'   c(year, period) for seasonal data.
#'
#' @return Numeric value at the specified date.
#'
#' @export
tsvalue <- function(x, date) {
  if (!is.ts(x)) stop("x must be a ts object")
  
  start <- start(x)
  freq <- frequency(x)
  n <- length(x)
  
  # Calculate position
  if (freq > 1) {
    if (length(date) < 2) stop("Date must include year and period for seasonal data")
    t <- (date[1] - start[1]) * freq + (date[2] - start[2]) + 1
  } else {
    t <- (date[1] - start[1]) + 1
  }
  
  if (t < 1 || t > n) {
    stop("Date is outside the time series range")
  }
  
  x[t]
}

#' Annual (rolling) sum
#'
#' Computes rolling sum over one year for monthly/quarterly data.
#'
#' @param x A \code{ts} object.
#' @param extend If TRUE, pads result with NAs to match original length.
#'   Default is TRUE.
#'
#' @return Time series of annual sums with same frequency as input.
#'
#' @export
S <- function(x, extend = TRUE) {
  if (is.mts(x)) {
    result <- sapply(seq_len(ncol(x)), function(k) S(x[, k], extend))
    return(ts(result, start = start(x), frequency = frequency(x)))
  }
  
  if (!is.ts(x)) stop("x must be a ts object")
  
  freq <- frequency(x)
  n <- length(x)
  
  if (n < freq) {
    warning("Series shorter than one year")
    return(ts(numeric(0), frequency = freq))
  }
  
  # Rolling sum
  y <- sapply(freq:n, function(k) sum(x[(k - freq + 1):k], na.rm = FALSE))
  
  if (extend) {
    y <- c(rep(NA, freq - 1), y)
  }
  
  ts(y, end = end(x), frequency = freq)
}

#' Format year-season labels
#'
#' Creates formatted labels combining year and seasonal period.
#'
#' @param x A \code{ts} object.
#' @param year.digits Number of digits for year (2 or 4). Default is 4.
#'
#' @return Character vector of formatted date labels.
#'
#' @noRd
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
