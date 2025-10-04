## tfarima/R/detvar.R
## Jose L Gallego (UC)

#' Calendar variables
#'
#' `CalendarVar()` creates a set of deterministic regressors to capture calendar
#' effects (trading/working days, length-of-month, leap-year and Easter).
#'
#' @param x A `ts` object used to determine start, length and frequency.
#' @param form Character selecting the set of calendar variables: `"dif"`
#'   (differences wrt reference day), `"td"` (6 dummies + lom; omits reference),
#'   `"td7"` (7 dummies), `"td6"` (6 dummies; omits reference), `"wd"` (weekdays
#'   vs weekend), or `"null"` (no trading-day regressors).
#' @param ref Non-negative integer (0–6) indicating the reference day (0 =
#'   Sunday, 1 = Monday, …, 6 = Saturday). Ignored unless `form` needs it.
#' @param lom Logical. If `TRUE` include a length-of-month regressor.
#' @param lpyear Logical. If `TRUE` include a leap-year regressor.
#' @param easter Logical. If `TRUE` include an Easter regressor.
#' @param len Integer duration for Easter effect (days). Typical values: 4–8.
#' @param easter.mon Logical. `TRUE` if Holy Monday is a public holiday.
#' @param n.ahead Integer. Extra observations to extend the sample (forecast
#'   horizon).
#'
#' @return An object of class `ts` or `mts` with the requested regressors.
#'
#' @references Bell, W.R. and Hillmer, S.C. (1983) “Modeling time series with
#'   calendar variation”, *Journal of the American Statistical Association*, 78,
#'   526–534.
#'
#' @examples
#' X <- CalendarVar(AirPassengers, form = "wd", easter = TRUE, len = 5)
#'
#' @export
CalendarVar <- function(
    x,
    form = c("dif", "td", "td7", "td6", "wd", "null"),
    ref = 0,
    lom = TRUE,
    lpyear = TRUE,
    easter = FALSE,
    len = 4,
    easter.mon = FALSE,
    n.ahead = 0
) {
  stopifnot(is.ts(x))
  stopifnot(frequency(x) == 12)
  form <- match.arg(form[1], c("dif", "td", "td7", "td6", "wd", "null"))
  nms <- c()
  xreg <- NULL
  if (form != "null") {  
    xreg <- extra.days(x, -1, n.ahead)
    nms <- colnames(xreg)
    if (ref < 0 || ref > 6) ref <- 0
    ref <- ref + 1
    if (form == "dif"||form == "td") {
      xreg <- sapply((1:7)[-ref], function(k) xreg[, k] - xreg[, ref])
      nms <- sapply((1:7)[-ref], function(k) {
        paste0(nms[k], "_", nms[ref])
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
      nms <- nms[-ref]
    } else if (form == "wd") {
      xreg <- matrix(rowSums(xreg[,2:6])-5*rowSums(xreg[,c(1, 7)])/2, ncol = 1)
      nms <- "wkdays_wknd"
      if (lpyear) lom <- FALSE
    } 
  }

  if (lom) {
    Lom <- month.length(x, n.ahead) - 28
    xreg <- cbind(xreg, Lom)
    nms <- c(nms, "lom")    
  }

  if (lpyear) {
    ly <- leapyear(x, wd = TRUE, n.ahead = n.ahead)
    xreg <- cbind(xreg, ly)
    nms <- c(nms, "leapyear")    
  }
  
  if (easter) {
    Easter <- EasterVar(x, len, easter.mon, n.ahead = n.ahead)
    xreg <- cbind(xreg, Easter)
    nms <- c(nms, paste0("Easter", len, ifelse(easter.mon, "M", "")))    
  }
  if (is.null(xreg))
    return(NULL)

  xreg <- as.matrix(xreg)
  colnames(xreg) <- nms
  ts(xreg, start = start(x), frequency = frequency(x))
}

#' Intervention variables
#'
#' `InterventionVar()` creates pulse, step, or ramp variables at a given date.
#'
#' @param Y A `ts` object used to determine start, length and frequency.
#' @param date Either a single positive index within the sample, or a vector
#'   `c(year, month)` for monthly series. For non-seasonal series, `c(year)` is
#'   also accepted.
#' @param type One of `"P"` (pulse), `"S"` (step), or `"R"` (ramp = cumulative
#'   step).
#' @param n.ahead Integer. Extra observations to extend the sample.
#'
#' @return A `ts` intervention variable.
#'
#' @references Box, G.E.P. and Tiao, G.C. (1975) “Intervention analysis with
#'   applications to economic and environmental problems”, *JASA*, 70(349),
#'   70–79.
#'
#' @examples
#' # Pulse at March 1958:
#' P <- InterventionVar(AirPassengers, date = c(1958, 3), type = "P")
#' # Or by index within the extended sample (here no extension):
#' P2 <- InterventionVar(AirPassengers, date = 123, type = "P")
#'
#' @export
InterventionVar <- function(Y, date, type = c("P", "S", "R"), n.ahead = 0) {
  stopifnot(is.ts(Y))
  type <- match.arg(type)
  start <- start(Y)
  s <- frequency(Y)
  nY <- length(Y) + n.ahead
  x <- ts(integer(nY), start = start, frequency = s)
  if (length(date) == 1 && s > 1) {
    n <- as.integer(date)
  } else if (s > 1L) {
    n <- (date[1] - start[1] + 1)*s - (start[2] - 1) - (s - date[2])
  } else {
    if (date[1] < start[1]) n <- as.integer(date[1])
    else n <- (date[1] - start[1] + 1)
  }
  stopifnot(n > 0 && n <= nY)
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
#' `sincos()` creates a full set of seasonal trigonometric regressors
#' (cos/sin harmonics) for a seasonal frequency.
#'
#' @param Y A seasonal `ts` object.
#' @param n.ahead Integer. Extra observations to extend the sample.
#' @param constant Logical. If `TRUE`, include an intercept column.
#'
#' @return A `ts` matrix with trigonometric variables (and intercept if requested).
#'
#' @examples
#' X <- sincos(AirPassengers, constant = TRUE)
#'
#' @export
sincos <- function(Y, n.ahead = 0, constant = FALSE) {
  stopifnot(is.ts(Y))
  start <- start(Y)
  s <- frequency(Y)
  stopifnot(s > 1)
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
  X <- as.matrix(X)
  colnames(X) <- names
  X <- ts(X, start = start, frequency = s)
}

#' Seasonal dummies
#'
#' `sdummies()` creates a full set of seasonal dummies (reference-coded).
#'
#' @param Y A seasonal `ts` object.
#' @param ref Reference season (1..frequency).
#' @param constant Logical. If `TRUE`, include an intercept column.
#' @param n.ahead Integer. Extra observations to extend the sample.
#'
#' @return A `ts` matrix with seasonal dummies (and intercept if requested).
#'
#' @examples
#' D <- sdummies(AirPassengers, ref = 1, constant = TRUE, n.ahead = 24)
#'
#' @export
sdummies <- function(Y, ref = 1, constant = FALSE, n.ahead = 0) {
  stopifnot(is.ts(Y))
  start <- start(Y)
  s <- frequency(Y)
  stopifnot(s > 1)
  stopifnot(ref >= 1 && ref <= s)
  n <- length(Y) + n.ahead
  m0 <- start[2]
  m <- ((m0 - 1):(m0 - 1 + n - 1)) %% s + 1
  d <- (m == ref)*1L
  i <- (1:s)[-ref]
  D <- sapply(i, function(k) (m==k) - d)
  names <- paste0("D", i, paste0("_D", ref))
  if (constant) {
    D <- cbind(rep(1L, n), D)
    names <- c("constant", names)
  }
  colnames(D) <- names
  D <- ts(D, start = start, frequency = s)
}

# --- helpers ------------------
#' @noRd
easterDate <- function(year, easter.mon = FALSE) {
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

#' @noRd
EasterVar <- function(x, len = 4, easter.mon = FALSE, n.ahead = 0) {
  stopifnot(is.ts(x))
  if (frequency(x) != 12) stop("function only implemented for monthly ts")
  start <- start(x)
  len <- abs(len)
  if (len < 1||len > 22) len <- 7 + (easter.mon*1L)
  n <- length(x) + n.ahead
  x <- ts(double(n), start = start, frequency = 12)
  end <- end(x)  
  # first year
  e <- easterDate(start[1], easter.mon)
  if (start[2] <= e[1]) {
    n <- e[1] - start[2] + 1
    if (e[1] == 3) x[n] <- 1
    else if (e[2] > len) x[n] <- 1
    else {
      x[n] <- e[2]/len
      if (n > 1) x[n-1] <- 1 - x[n]
    }
  }
  
  # middle years
  for (y in (start[1]+1):(end[1]-1)) {
    e <- easterDate(y, easter.mon)
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
    e <- easterDate(y, easter.mon)
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

#' @noRd
month.length <- function(Y, n.ahead = 0) {
  stopifnot(is.ts(Y))
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

#' @noRd
leapyear <- function(Y, wd = FALSE, n.ahead = 0) {
  stopifnot(is.ts(Y))
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

#' @noRd
extra.days <- function(x, day, n.ahead = 0) {
  stopifnot(is.ts(x))
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

