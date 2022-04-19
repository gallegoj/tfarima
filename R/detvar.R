## tfarima/R/detvar.R
## Jose L Gallego (UC)

#' Calendar variables
#'
#' \code{CalendarVar} creates a set of deterministic variables to capture
#' calendar effects.
#'
#' @param x an object of class \code{ts} used to determine the sample period and
#'   frequency.
#' @param form a character indicated the set of calendar variables: td, td7,
#'   td6, wd.
#' @param ref a non-negative integer indicating the reference day.
#' @param lom logical. If TRUE length of the month effect is also estimated.
#' @param lpyear logical. If TRUE a leap year effect is also estimated.
#' @param easter logical. If TRUE an additional deterministic variable is
#'   generated to capture Easter effects.
#' @param len duration of the Easter, integer.
#' @param easter.mon logical. It is TRUE if Holy Monday is a public holiday.
#' @param n.ahead number of additional observations to extend the sample period.
#'
#' @return An object of class \code{mts} or \code{ts}.
#' @references
#'
#' Bell, W.R. and Hillmer, S.C. (1983) “Modeling time series with calendar
#' variation”, Journal of the American Statistical Society, Vol. 78, pp.
#' 526–534.
#'
#' @examples
#'
#' Y <- rsales
#' X <- CalendarVar(Y, easter = TRUE)
#'
#' @export
CalendarVar <- 
function(x, form = c("dif", "td", "td7", "td6", "wd",  "wd2", "null"), ref = 0,
         lom = TRUE, lpyear = TRUE, easter = FALSE, len = 4, easter.mon = FALSE,
         n.ahead = 0) 
{

  form <- tolower(form)[1]
  names <- c()
  if (form != "null") {  
    xreg <- extra.days(x, -1, n.ahead)
    names <- colnames(xreg)
    if (ref < 0 || ref > 6) ref <- 0
    ref <- ref + 1
    if (form == "dif"||form == "td") {
      xreg <- sapply((1:7)[-ref], function(k) xreg[, k] - xreg[, ref])
      names <- sapply((1:7)[-ref], function(k) {
        paste0(names[k], "_", names[ref])
      })
      if (lom && lpyear) {
        if (form == "td") lom <- FALSE
        else lpyear <- FALSE
      }
    } else if (form == "td7") {
      lom <- FALSE
      lpyear <- FALSE
    } else if (form == "td6") {
      xreg <- xreg[, -ref]
      names <- names[-ref]
    } else if (form == "wd") {
      xreg <- matrix(rowSums(xreg[,2:6])-5*rowSums(xreg[,c(1, 7)])/2, ncol = 1)
      names <- "wkdays_wknd"
      if (lpyear) lom <- FALSE
    } else {
      stop("invalid option for tranding days")
    }
  } else xreg <- NULL

  if (lom) {
    Lom <- month.length(x, n.ahead) - 28
    xreg <- cbind(xreg, Lom)
    names <- c(names, "lom")    
  }

  if (lpyear) {
    ly <- lpyear(x, TRUE, n.ahead = n.ahead)
    xreg <- cbind(xreg, ly)
    names <- c(names, "lpyear")    
  }
  
  if (easter) {
    Easter <- EasterVar(x, len, easter.mon, n.ahead = n.ahead)
    xreg <- cbind(xreg, Easter)
    names <- c(names, paste0("Easter", len, ifelse(easter.mon, "M", "")))    
  }

  if (ncol(xreg) > 1)
    xreg <- ts(xreg, start = start(x), frequency = frequency(x))
  colnames(xreg) <- names
  xreg
}


#' Intervention variables
#'
#' \code{InterventionVar} creates an intervention variable to capture the effect
#' of an external event.
#'
#' @param Y an object of class \code{ts} used to determine the sample period and
#'   frequency.
#' @param type a character indicating the type of intervention variables: (P)
#'   pulse, (S) step, (R).
#' @param date the date of the event, c(year, month).   
#' @param n.ahead number of additional observations to extend the sample period.
#'
#' @return An intervention variable, a 'ts' object.
#'
#' @references
#'
#' G. E. P. Box, G. C. Tiao, “Intervention Analysis with Applications to
#' Economic and Environmental Problems”, Journal of the American Statistical
#' Association, Vol. 70, No. 349. (Mar., 1975), pp. 70-79.
#'
#' @examples
#' 
#' Y <- seriesJ$Y
#' P58 <- InterventionVar(Y, date = 58, type = "P")
#'
#' @export
InterventionVar <- function(Y, date, type = c("P", "S", "R"), n.ahead = 0) {
  type <- match.arg(type)
  start <- start(Y)
  s <- frequency(Y)
  nY <- length(Y) + n.ahead
  x <- ts(integer(nY), start = start, frequency = s)
  if (s > 1) n <- (date[1] - start[1] + 1)*s - (start[2] - 1) - (s - date[2])
  else n <- (date[1] - start[1] + 1)*s
  stopifnot(n > 0 || n <= nY)
  x[n] <- 1
  if (type != "P") {
    x <- cumsum(x)
    if (type == "R") x <- cumsum(x)
    x <- ts(x, start = start, frequency = s)
  }
  x
}

#' Trigonometric variables
#'
#' \code{sincos} creates an full set of trigonometric variables.
#'
#' @param Y an object of class \code{ts} used to determine the sample period and
#'   frequency.
#' @param constant logical indicator to include a column of ones.
#' @param n.ahead number of additional observations to extend the sample period.
#'
#' @return A matrix of trigonometric variables.
#'
#' @examples
#'
#' Y <- AirPassengers
#' P58 <- sincos(Y)
#'
#' @export
sincos <- function(Y, n.ahead = 0, constant = FALSE) {
  start <- start(Y)
  s <- frequency(Y)
  stopifnot(s>1)
  n <- length(Y) + n.ahead
  t <- 2*pi*(start[2]:(start[2] + n - 1))/s
  k = floor((s-1)/2)
  if (constant) {
    X <- rep(1, n)
    names <- c("constant")
  } else {
    X <- NULL
    names <- c()
  }
  if (s > 2) {
    for (j in 1:k) { 
        X <- cbind(X, cos(j*t), sin(j*t))
        names <- c(names, paste0("c", j), paste0("s", j))
    }
  }
  if (s %% 2 == 0) {
    X <- cbind(X, cos(s*t/2))
    names <- c(names, paste0("c", s/2))
  }
  colnames(X) <- names
  X <- ts(X, start = start, frequency = s)
}

#' Seasonal dummies
#'
#' \code{sdummies} creates an full set of seasonal dummies.
#'
#' @param Y an object of class \code{ts} used to determine the sample period and
#'   frequency.
#' @param ref the reference season, positive integer
#' @param constant logical indicator to include a column of ones.
#' @param n.ahead number of additional observations to extend the sample period.
#'
#' @return A matrix of trigonometric variables.
#'
#' @examples
#'
#' Y <- AirPassengers
#' P58 <- sincos(Y)
#'
#' @export
sdummies <- function(Y, ref = 1, constant = FALSE, n.ahead = 0) {
  start <- start(Y)
  s <- frequency(Y)
  stopifnot(s>1)
  stopifnot(ref >=1 && ref <= s)
  n <- length(Y) + n.ahead
  m <- cycle(Y)
  d <- (m == ref)*1L
  i <- (1:s)[-ref]
  D <- sapply(i, function(k) (m==k) - d)
  names <- paste0("D", i, paste0("_D", ref))
  if (constant) {
    D <- cbind(rep(1, n), D)
    names <- c("constant", names)
  }
  colnames(D) <- names
  D <- ts(D, start = start, frequency = s)
}


easter.date <- function(year, easter.mon = FALSE) {
  a <- year %% 19
  b <- year %% 4
  c <- year %% 7
  d <- (19 * a + 24) %% 30
  e <- (2 * b + 4 * c + 6 * d + 5) %%  7
  day <- 22 + d + e
  if (easter.mon) day <- day + 1
  if (day <= 31) month = 3
  else {
    day <- day - 31
    month <- 4
  }
  c(month, day)
}

EasterVar <- function(x, len = 4, easter.mon = FALSE, n.ahead = 0) {
  if (frequency(x) != 12) stop("function only implemented for monthly ts")
  start <- start(x)
  len <- abs(len)
  if (len < 1||len > 22) len <- 7 + (easter.mon*1L)
  n <- length(x) + n.ahead
  x <- ts(double(n), start = start, frequency = 12)
  end <- end(x)  
  # first year
  e <- easter.date(start[1], easter.mon)
  if (start[2] <= e[1]) {
    n <- e[1] - start[2] + 1
    if (e[1] == 3) x[n] <- 1
    else if (e[2] > len) x[n] <- 1
    else {
      x[n] <- e[2]/len
      if (n > 1) x[n-1] <- 1 - x[n]
    }
  }

  for (y in (start[1]+1):(end[1]-1)) {
    e <- easter.date(y, easter.mon)
    n <- (y - start[1] + 1)*12 - (start[2] - 1) - (12 - e[1])
    if (e[1] == 3) x[n] <- 1
    else if (e[2] > len) x[n] <- 1
    else {
      x[n] <- e[2]/len
      x[n-1] <- 1 - x[n]
    }
  }
  
  # last year
  if (end[2] > 2) {
    y <- end[1]
    e <- easter.date(y, easter.mon)
    n <- (y - start[1] + 1)*12 - (start[2] - 1) - (12 - e[1])
    if (e[1] <= end[2]) {
      if (e[1] == 3) x[n] <- 1
      else if (e[2] > len) x[n] <- 1
      else {
        x[n] <- e[2]/len
        x[n-1] <- 1 - x[n]
      }
    } else if (end[2] == 3 && e[2] < len) {
      x[n] <- 1- e[2]/len
    }
  }

  x  

}

month.length <- function(Y, n.ahead = 0) {
  if (frequency(Y) != 12) stop("function only implemented for monthly ts")
  start <- start(Y)
  n <- length(Y) + n.ahead
  x <- ts(double(n), start = start, frequency = 12)
  
  y <- start[1]
  m <- start[2]
  for (t in 1:n) {
    if (m == 2) {
      if( (y%%400 == 0) || (y%%100 == 0) || (y%%4 == 0)) x[t] <- 29
      else x[t] <- 28
    } else if (m == 4 || m == 6 || m == 9 || m == 11) x[t] <- 30
    else x[t] <- 31
    m <- m + 1
    if (m > 12) {
      m <- 1
      y <- y +1
    }
  }
  
  x

}

lpyear <- function(Y, wd = FALSE, n.ahead = 0) {
  if (frequency(Y) != 12) stop("function only implemented for monthly ts")
  start <- start(Y)
  n <- length(Y) + n.ahead
  x <- ts(double(n), start = start, frequency = 12)
  y <- start[1]
  m <- start[2]
  for (t in 1:n) {
    if (m == 2) {
      if (wd) {
        if( (y%%400 == 0) || (y%%100 == 0) || (y%%4 == 0)) x[t] <- 0.75
        else x[t] <- -0.25
      } else {
        if( (y%%400 == 0) || (y%%100 == 0) || (y%%4 == 0)) x[t] <- 1
      }
    } 
    m <- m + 1
    if (m > 12) {
      m <- 1
      y <- y +1
    }
  }
  x
}

extra.days <- function(x, day, n.ahead = 0) {
  if (day < 0 || day > 6) {
    X <- sapply(0:6, function(k) {
      extra.days(x, k, n.ahead)
    })
    colnames(X) <- c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")
    return(ts(X, start = start(x), frequency = 12))
  }
  if (frequency(x) != 12) stop("function only implemented for monthly ts")
  start <- start(x)
  n <- length(x) + n.ahead
  x <- ts(rep(4, n), start = start, frequency = 12)
  
  y <- start[1]
  m <- start[2]
  for (t in 1:n) {
    a <- (14 - m) %/% 12
    b <- y - a
    c <- m + 12*a - 2
    d <- (1 + b + (b%/%4) - (b%/%100) + (b%/%400) + ((31*c)%/%12)) %% 7
    if (m == 2L) {
      if (day == d) {
        if( (y%%400 == 0) || (y%%100 == 0) || (y%%4 == 0)) x[t] <- 5
      }
    } else if (m == 4L ||m == 6L ||m == 9L ||m == 11L) {
      if (d == 6L) {
        if(day == 6L || day == 0L) x[t] <- 5
      } else {
        if (day == d || day == d + 1) x[t] <- 5
      }
    } else {
      if (d == 5L) {
        if(day == 5L || day == 6L || day == 0L) x[t] <- 5
      } else if (d == 6L) {
        if(day == 6L || day == 0L || day == 1L) x[t] <- 5
      } else{
        if(day == d || day == d + 1 || day == d + 2) x[t] <- 5
      }         
    }
    m <- m + 1
    if (m > 12) {
      m <- 1
      y <- y + 1
    }
  }
  x - 4
}

