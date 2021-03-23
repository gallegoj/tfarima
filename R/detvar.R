## tfarima/R/detvar.R
## Jose L Gallego (UC)

#' Calendar variables
#'
#' \code{CalendarVar} creates a set of deterministic variables to capture
#' calendar effects.
#'
#' @param x an object of class \code{ts} used to determine the sample period and
#'   frequency.
#' @param form a character indicated the set of calendar variables.
#' @param easter logical. If TRUE an additional deterministic variable is
#'   generated to capture Easter effects.
#' @param n.ahead number of additional observations to extend the sample period.
#'
#' @return A matrix of explanatory variables.
#'
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
function(x, form = c("dif", "lom", "td"), easter = FALSE, n.ahead = 0) 
{
  form <- tolower(form)[1]
  Sun <- num.days(x, 0, n.ahead = n.ahead) - 4
  Mon <- num.days(x, 1, n.ahead = n.ahead) - 4
  Tue <- num.days(x, 2, n.ahead = n.ahead) - 4
  Wed <- num.days(x, 3, n.ahead = n.ahead) - 4
  Thu <- num.days(x, 4, n.ahead = n.ahead) - 4
  Fri <- num.days(x, 5, n.ahead = n.ahead) - 4
  Sat <- num.days(x, 6, n.ahead = n.ahead) - 4
  Lom <- month.length(x, n.ahead) - 28
  
  if (form == "td") {
    xreg <- data.frame(Sun, Mon, Tue, Wed, Thu, Fri, Sat)
  } else if (form == "lom") {
    xreg <- data.frame(Lom, Mon, Tue, Wed, Thu, Fri, Sat)
  } else {
    xreg <- data.frame(Lom, Mon_Sun = Mon - Sun, Tue_Sun = Tue - Sun,
                       Wed_Sun = Wed - Sun, Thu_Sun = Thu - Sun,
                       Fri_Sun = Fri - Sun, Sat_Sun = Sat - Sun)
  }
  if (easter) {
    Easter <- EasterVar(x, n.ahead = n.ahead)
    xreg <- data.frame(Easter, xreg)
  }
  
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
  frequency <- frequency(Y)
  nY <- length(Y) + n.ahead
  x <- ts(integer(nY), start = start, frequency = frequency)
  s <- frequency
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

easter.date <- function(year) {
  a <- year %% 19
  b <- year %% 4
  c <- year %% 7
  d <- (19 * a + 24) %% 30
  e <- (2 * b + 4 * c + 6 * d + 5) %%  7
  day <- 22 + d + e
  if (day <= 31) month = 3
  else {
    day <- day -31
    month <- 4
  }
  
  c(month, day)
  
}

EasterVar <- function(Y, len = 4, n.ahead = 0) {
  if (frequency(Y) != 12) stop("function only implemented for monthly ts")
  start <- start(Y)
  n <- length(Y) + n.ahead
  x <- ts(double(n), start = start, frequency = 12)
  end <- end(x)  
  # first year
  e <- easter.date(start[1])
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
    e <- easter.date(y)
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
    e <- easter.date(y)
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

num.days <- function(Y, day, n.ahead = 0) {
  if (frequency(Y) != 12) stop("function only implemented for monthly ts")
  start <- start(Y)
  n <- length(Y) + n.ahead
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
        if(day == 5L || day == 6L || d==0L) x[t] <- 5
      } else if (d == 6L) {
        if(day == 6L || day == 0L || day == 1L) x[t] <- 5
      } else{
        if(day == d || day == d + 1 || day == d + 2) x[t] <- 5
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

