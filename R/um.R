## tfarima/R/um.R
## Jose L Gallego (UC)

#' Univariate (ARIMA) model
#'
#' \code{um} creates an S3 object representing a univariate ARIMA model, which
#' can contain multiple AR, I and MA polynomials, as well as parameter
#' restrictions.
#'
#' @param z an object of class \code{ts}.
#' @param ar list of stationary AR lag polynomials.
#' @param i list of nonstationary AR (I) polynomials.
#' @param ma list of MA polynomials.
#' @param mu mean of the stationary time series.
#' @param sig2 variance of the error.
#' @param bc logical. If TRUE logs are taken.
#' @param fit logical. If TRUE, model is fitted.
#' @param envir the environment in which to look for the time series z when it
#'   is passed as a character string.
#' @param warn logical. If TRUE, a warning is displayed for non-admissible
#'   models.
#' @param ... additional arguments.
#'
#' @return An object of class \code{um}.
#'
#' @references
#'
#' Box, G.E.P., Jenkins, G.M., Reinsel, G.C. and Ljung, G.M. (2015) Time Series
#' Analysis: Forecasting and Control. John Wiley & Sons, Hoboken.
#'
#' @examples
#'
#' ar1 <- um(ar = "(1 - 0.8B)")
#' ar2 <- um(ar = "(1 - 1.4B + 0.8B^2)")
#' ma1 <- um(ma = "(1 - 0.8B)")
#' ma2 <- um(ma = "(1 - 1.4B + 0.8B^2)")
#' arma11 <- um(ar = "(1 - 1.4B + 0.8B^2)", ma = "(1 - 0.8B)")
#'
#' @export
um <- function(z = NULL, ar = NULL, i = NULL, ma = NULL, mu = NULL, sig2 = 1.0, 
               bc = FALSE, fit = TRUE, envir = parent.frame(), warn = TRUE,
               ...) {

  call <- match.call()
  if (is.numeric(z)){
    z <- deparse(substitute(z))
  }

  if (!is.null(ar)) {
    ar <- lagpol0(ar, "ar", envir = envir)
    phi <- polyexpand(ar)
  } else {
    phi <- 1.0
  }
  names(phi) <- paste0("[phi", 0:(length(phi)-1), "]")

  if (!is.null(i)) {
    i <- lagpol0(i, "i", envir)
    nabla <- polyexpand(i)
  } else {
    nabla <- 1.0
  }

  if (!is.null(ma)) {
    ma <- lagpol0(ma, "ma", envir = envir)
    theta <- polyexpand(ma)
  } else {
    theta <- 1.0
  }
  names(theta) <- paste("[theta", 0:(length(theta)-1), "]", sep = "")
  param <- lapply(c(ar, ma), function(x) unlist(x$param))
  param <- unname(param)
  param <- unlist(param)
  param <- param[!duplicated(names(param))]
  if (!is.null(mu)) {
    if (!all(names(param) != "mu")) 
      stop("'mu' can not be used as an ARMA parameter")
    param <- c(mu = mu, param)
  }
  param <- as.list(param)

  if (is.null(names(sig2))) names(sig2) <- "sig2"
  
  mod <- list(z = z, call = call, phi = phi, nabla = nabla, theta = theta,
              mu = mu, sig2 = sig2, bc = bc, ar = ar, i = i, ma = ma,  
              param = param, k = length(param), kar = length(ar), 
              kma = length(ma),
              p = length(phi)-1, d = length(nabla)-1, q = length(theta)-1,
              optim = NULL, method = NULL, is.adm = TRUE)
  class(mod) <- "um"
  mod <- .update_um(mod, unlist(param, use.names = TRUE))
  if (mod$is.adm) {
    if (!is.null(z) && fit) mod <- fit.um(mod, envir=envir, ...)
    else mod$b <- param.um(mod)
  } else {
    if (warn) warning("non-admisible model")
  }
  
  return(mod)
  
}

#' Airline Model (SARIMA(0,1,1)x(0,1,1)s)
#'
#' Creates a seasonal ARIMA model with the structure popularized by Box and
#' Jenkins using airline passenger data: (0,1,1)x(0,1,1)s.
#'
#' @param z A \code{ts} object (must have frequency > 1).
#' @param bc Logical. If TRUE, applies Box-Cox (log) transformation.
#' @param sma Character. Specification for seasonal MA operator. Options are:
#'   standard, generalized or factorized. See manual for more details.
#' @param ... Additional arguments passed to \code{\link{um}}.
#'
#' @return A \code{um} object with airline model specification.
#'
#' @details This is a convenience function equivalent to: \code{um(z, bc = bc, i
#' = list(1, c(1, s)), ma = list(1, c(1, s)), ...)} where s = frequency(z).
#'
#' @seealso \code{\link{um}}
#' @export
airline <- function(z, bc = FALSE, 
                    sma = c("standard", "generalized", "factorized"), ...) {
  z.name <- deparse(substitute(z))
  sma <- match.arg(sma)
  if (!is.ts(z)) {
    stop("'z' must be a ts object")
  }
  s <- frequency(z)
  if (s < 2) {
    stop("Airline model requires seasonal data (frequency > 1)")
  }
  if (startsWith("standard", sma)) {
    um1 <- um(z, bc = bc, i = list(1, c(1, s)), 
              ma = list(theta = 1, THETA =c(1, s)), ...)
  } else if (startsWith("generalized", sma)) {
    um1 <- um(z, bc = bc, i = list(1, c(1, s)), 
              ma = list(theta = 2, THETA = paste0(s)), ...)
  } else {
    THETA = paste0("(1:", floor(s/2), ")/", s)
    um1 <- um(z, bc = bc, i = list(1, c(1, s)), 
              ma = list(theta = 2, THETA = THETA), ...)
  }
  um1$z <- z.name
  um1
}


#' Convert \code{arima} into \code{um}.
#'
#' \code{as.um} converts an object of class \code{arima} into an object 
#' of class \code{um}.
#'
#' @param arima an object of class \code{arima}.
#' @param ... additional arguments.
#' 
#' @return 
#' An object of class \code{um}.
#' 
#' @examples
#' z <- AirPassengers
#' a <- arima(log(z), order = c(0,1,1), 
#' seasonal = list(order = c(0,1,1), frequency = 12))
#' um1 <- as.um(a)
#' 
#' @export
as.um <- function(arima, ...) {
  if (is.ssm(arima))
    return(.ssm2um(arima, ...))
  else if(is.uc(arima))
    return(.uc2um(arima, ...))
  else if(inherits(arima, "ucarima"))
    return(.ucarima2um(arima, ...))
  
  b <- arima$coef
  arma <- arima$arma
  
  p <- arma[1]
  if ( p > 0 ) {
    param <- as.list(b[1:p])
    names(param) <- names(arima$coef)[1:p]
    ar <- lagpol(param = param)
  } else ar <- NULL
  
  q <- arma[2]
  if ( q > 0 ) {
    param <- as.list(-b[(p+1):(p+q)])
    names(param) <- names(arima$coef)[(p+1):(p+q)]
    ma <- lagpol(param = param)
  } else ma <- NULL
  
  s <- arma[5]
  if ( s > 0) {
    P <- arma[3]
    if ( P > 0 ) {
      param <- as.list(b[(p+q+1):(p+q+P)])
      names(param) <- names(arima$coef)[(p+q+1):(p+q+P)]
      ar <- list(ar, lagpol(param = param, s = s))
    }
    Q <- arma[4]
    if ( Q > 0 ) { 
      param <- as.list(-b[(p+q+P+1):(p+q+P+Q)])
      names(param) <- names(arima$coef)[(p+q+P+1):(p+q+P+Q)]
      ma <- list(ma, lagpol(param = param, s = s))
    }
  }
  
  d <- arma[6]
  D <- arma[7]
  if (d > 0 && D > 0) i <- list(d, c(D, s))
  else if (d > 0) i <- d
  else if (D > 0) i <- c(D, s)
  else i <- NULL
  
  if (any(names(arima$coef) == "intercept")) {
    mu <- unname(arima$coef["intercept"])
  } else mu <- NULL

  um <- um(ar = ar, i = i, ma = ma, mu = mu, sig2 = arima$sigma2)
  um$z <- arima$series
  um
}

#' Theoretical autocovariances of an ARMA model
#'
#' \code{autocov} computes the autocovariances of an ARMA model.
#'
#' @param mdl an object of class \code{um} or \code{ucm}.
#' @param lag.max maximum lag for autocovariances.
#' @param ... additional arguments.
#'  
#' @return 
#' A numeric vector.
#' 
#' @note 
#' The I polynomial is ignored.
#' 
#' @examples
#' ar1 <- um(ar = "1-0.8B")
#' autocov(ar1, lag.max = 13)
#' 
#' @export
autocov <- function (mdl, ...) { UseMethod("autocov") }

#' @rdname autocov
#' @export
autocov.um <- function(mdl, lag.max = 10, ...) {
  stopifnot(inherits(mdl, "um"))
  g <- as.numeric( tacovC(mdl$phi, mdl$theta, mdl$sig2, lag.max) )
  names(g) <- paste("gamma", 0:lag.max, sep = "")
  g
}

#' Theoretical simple/partial autocorrelations of an ARMA model
#'
#' \code{autocorr} computes the simple/partial autocorrelations of an ARMA model.
#'
#' @param x an object of class \code{um}.
#' @param lag.max maximum lag for autocovariances.
#' @param par logical. If TRUE partial autocorrelations are computed.
#' @param ... additional arguments.
#' 
#' @return 
#' A numeric vector.
#' 
#' @note 
#' The I polynomial is ignored.
#' 
#' @examples
#' ar1 <- um(ar = "1-0.8B")
#' autocorr(ar1, lag.max = 13)
#' autocorr(ar1, lag.max = 13, par = TRUE)
#' 
#' @export
autocorr <- function (x, ...) { UseMethod("autocorr") }

#' @rdname autocorr
#' @export
autocorr.um <- function(x, lag.max = 10, par = FALSE, ...) {
  stopifnot(inherits(x, "um"))
  if (!par) {
    g <- as.numeric( tacovC(x$phi, x$theta, x$sig2, lag.max) )
    names(g) <- paste("rho", 0:lag.max, sep = "")
    g <- g/g[1]
    g
  } else {
    g <- as.numeric( pacorrC(x$phi, x$theta, lag.max) )
    names(g) <- paste("phi", 1:lag.max, ".", 1:lag.max, sep = "")
    g
  }
}

#' Coefficients of a univariate model
#'
#' \code{coef} extracts the "coefficients" from a um object.
#'
#' @param object a \code{um} object.
#' @param ... other arguments.
#'
#' @return A numeric vector.
#'
#' @export
coef.um <- function(object, ...) {
  param.um(object)
}

##' Add or Replace Inputs in Models
#'
#' Adds new inputs to transfer function or univariate models.
#'
#' @param mdl A \code{um} or \code{tfm} object.
#' @param xreg Optional matrix of exogenous regressors.
#' @param inputs Optional list of \code{tf} objects (only for \code{tfm}).
#' @param y Optional \code{ts} object for output series.
#' @param envir Environment for evaluation. Default is calling environment.
#' @param ... Additional arguments passed to model constructor.
#' 
#' @return A \code{tfm} object.
#' 
#' @seealso \code{\link{um}}, \code{\link{tfm}}
#' 
#' @export
setinputs <- function (mdl, ...) { UseMethod("setinputs") }

#' @rdname setinputs
#' @export
setinputs.um <- function(mdl, xreg = NULL, inputs = NULL, y = NULL, 
                          envir = NULL, ...) {
  if (is.null (envir)) envir <- parent.frame ()
  if (!is.null(y)) mdl$z <- deparse(substitute(y))
  y <- z.um(mdl, z = y, envir = envir)
  tfm1 <- tfm(output = y, noise = mdl, fit = FALSE, new.name = FALSE)
  setinputs.tfm(tfm1, xreg, inputs)  
}

#'Calendar effects
#'
#'\code{calendar} extends the ARIMA model \code{um} by including a set of
#'deterministic variables to capture the calendar variation in a monthly time
#'series. Two equivalent representations are available: (i) D0, D1, ..., D6,
#'(ii) L, D1-D0, ..., D6-D0 where D0, D2, ..., D6 are deterministic variables
#'representing the number of Sundays, Mondays, ..., Saturdays, L = D0 + D1 + ...
#'+ D6 is the of the month. Alternatively, the Leap Year indicator (LPY) can be
#'included instead of L. The seven trading days can also be compacted into two
#'variables: week days and weekends. Optionally, a deterministic variable to
#'estimate the Easter effect can also be included, see "\code{\link{easter}}".
#'
#'@param mdl an object of class \code{\link{um}} or \code{\link{tfm}}.
#'@param y a time series.
#'@param form representation for calendar effects: (1) \code{form = dif}, L,
#'  D1-D0, ..., D6-D0; (2) \code{form = td}, LPY, D1-D0, ..., D6-D0; (3)
#'  \code{form = td7}, D0, D2, ..., D6; (4) \code{form = td6}, D1, D2, ..., D6;
#'  (5) \code{form = wd}, (D1+...+D5) - 2(D6+D0)/5.
#'@param ref a integer indicating the the reference day. By default, ref = 0.
#'@param lom,lpyear a logical value indicating whether or not to include the
#'  lom/lead year indicator.
#'@param easter logical. If \code{TRUE} an Easter effect is also estimated.
#'@param len the length of the Easter, integer.
#'@param easter.mon logical. TRUE indicates that Easter Monday is a public
#'  holiday.
#'@param n.ahead a positive integer to extend the sample period of the
#'  deterministic variables with \code{n.ahead} observations, which could be
#'  necessary to forecast the output.
#'@param p.value estimates with a p-value greater than p.value are omitted.
#'@param envir environment in which the function arguments are evaluated. If
#'  NULL the calling environment of this function will be used.
#'@param ... other arguments.
#'
#'@return An object of class "\code{\link{tfm}}".
#'
#'@references W. R. Bell & S. C. Hillmer (1983) Modeling Time Series with
#'  Calendar Variation, Journal of the American Statistical Association, 78:383,
#'  526-534, DOI: 10.1080/01621459.1983.10478005
#'
#' @examples
#' Y <- tfarima::rsales
#' um1 <- um(Y, i = list(1, c(1, 12)), ma = list(1, c(1, 12)), bc = TRUE)
#' tfm1 <- calendar(um1)
#'
#'@export
calendar <- function (mdl, ...) { UseMethod("calendar") }

#' @rdname calendar
#' @export
calendar.um <-
function(mdl, y = NULL, form = c("dif", "td", "td7", "td6", "wd"),
         ref = 0, lom = TRUE, lpyear = TRUE, easter = FALSE, len = 4, 
         easter.mon = FALSE, n.ahead = 0, p.value = 1, 
         envir = parent.frame (), ...)
{
  if (is.null(y)) y <- z.um(mdl, envir = envir)
  else mdl$z <- deparse(substitute(y))
  if (frequency(y) != 12) stop("function only implemented for monthly ts")

  n.ahead <- abs(n.ahead)
  xreg <- CalendarVar(y, form, ref, lom, lpyear, easter, len, easter.mon, n.ahead)
  tfm1 <- tfm(y, xreg = xreg, noise = mdl, new.name = FALSE, envir = envir, ...)
  if (p.value < 0.999) {
    p <- summary.tfm(tfm1, p.values = TRUE)
    p <- (p[1:ncol(xreg)] <= p.value)
    if (all(p)) return(tfm1)
    if (any(p)) {
      xreg <- xreg[, p]
      tfm1 <- tfm(y, xreg = xreg, noise = tfm1$noise, envir=envir, new.name = FALSE, ...)
      return(tfm1)
    }
    return(mdl)
  }
  return(tfm1)
}

#' Diagnostic checking
#'
#' \code{diagchk} displays tools for diagnostic checking.
#'
#' @param mdl an object of class \code{um} or \code{tfm}.
#' @param ... additional arguments.
#'
#' @export
diagchk <- function (mdl, ...) { UseMethod("diagchk") }


#' @rdname diagchk
#' @param mdl an object of class \code{um}.
#' @param z optional, an object of class \code{ts}.
#' @param method character; "exact" or "conditional" residuals.
#' @param lag.max integer; maximum number of lags for ACF/PACF.
#' @param lags.at numeric vector; specific lags in ACF/PACF plots.
#' @param freq.at numeric vector; specific frequencies in (cum) periodogram plot.
#' @param std logical; if TRUE standardized residuals are used.
#' @param envir environment in which the function arguments are evaluated.
#'    If NULL the calling environment of this function will be used.
#'    
#' @examples
#' z <- AirPassengers
#' airl <- um(z, i = list(1, c(1,12)), ma = list(1, c(1,12)), bc = TRUE)
#' diagchk(airl)
#' 
#' @export
diagchk.um <- function(mdl, z = NULL, method = c("exact", "cond"),
                       lag.max = NULL, lags.at = NULL, freq.at = NULL,
                       std = TRUE, envir=NULL, ...) {
  if (is.null (envir)) envir <- parent.frame ()
  u <- residuals.um(mdl, z, envir=envir)
  ide(u, graphs = c("plot", "hist", "acf", "pacf", "cpgram"), ylab = "u",
      lag.max = lag.max, lags.at = lags.at, freq.at = freq.at,
      std = std, envir = envir, ...)
}


#' Graphs for ARMA models
#'
#' \code{display} shows graphs characterizing one or a list of ARMA models.
#'
#' @param um an object of class \code{um} or a list of these objects.
#' @param lag.max number of lags for ACF/PACF.
#' @param lags.at the lags of the ACF/PACF at which tick-marks are to be drawn.
#' @param n.freq number of frequencies for the spectrum.
#' @param log.spec logical. If TRUE log spectrum is computed.
#' @param graphs vector of graphs.
#' @param byrow orientation of the graphs.
#' @param eq logical. If TRUE the model equation is used as title.
#' @param cex double. Font size for equation text.
#' @param ... additional arguments.
#' 
#' @examples
#' um1 <- um(ar = "(1 - 0.8B)(1 - 0.8B^12)")
#' um2 <- um(ma = "(1 - 0.8B)(1 - 0.8B^12)")
#' display(list(um1, um2))
#' 
#' @export 
display <- function (um, ...)
UseMethod("display")

#' @rdname display
#' @export
display.um <- 
function(um, lag.max = 25, n.freq = 501, log.spec = FALSE, lags.at = NULL,  
         graphs = c("acf", "pacf", "spec"), byrow = FALSE, eq = TRUE,
         cex = 1.25, ...) 
{
  if (inherits(um, "um")) {
    n.um <- 1
    um <- list(um)
  } else if( all(sapply(um, is.um)) ) {
    n.um <- length(um)
  } else {
    stop("Invalid um object")    
  }
  
  graphs <- tolower(graphs)
  graphs <- unique(graphs)
  graphs <- match.arg(graphs, c("acf", "pacf", "spec"), several.ok = TRUE)
  n.graphs <- length(graphs)
  if (n.graphs < 1) {
    graphs <- c("acf", "pacf", "spec")
    n.graphs <- 3
  }

  if (lag.max < 1) lag.max <- 25
  if (n.freq < 10) n.freq <- 501
  
  if (!is.null(lags.at)) {
    if (length(lags.at) == 1 && lags.at[1] > 1) {
      lags.at = seq(lags.at, lag.max, lags.at)
    } 
  }  
  
  Lag <- 1:lag.max
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  m <- matrix(1:(n.um*n.graphs), nrow = n.graphs, ncol = n.um, byrow = FALSE)
  if (byrow) m <- t(m)
  
  if (eq) {
    if (byrow){
      m <- cbind(1:n.um, m + n.um)
      layout(m, widths = c(lcm(1), rep(1, n.graphs+1)))
      rot = 90
    } else{
      m <- rbind(1:n.um, m + n.um)
      layout(m, heights = c(lcm(1), rep(1, n.graphs+1)))
      rot = 0
    }
    par(mar = c(0, 0, 0, 0))
    for (i in 1:n.um) {
      plot.new()
      text(.5, .5, parse(text = eq.um(um[[i]])), cex = cex,
           font = 2, srt = rot)
    }
  }
  else {
    layout(m)
  }
  par(mar = c(3, 3, 1, 0.8), mgp = c(1.5, 0.6, 0))
  
  for (i in 1:n.um) {
    mod <- um[[i]]
    for (j in 1:n.graphs) {
      if (graphs[j] == "acf") {
        ACF <- as.numeric( tacovC(mod$phi, mod$theta, mod$sig2, lag.max) )
        ACF <- ACF[-1]/ACF[1]
        if (!is.null(lags.at)) {
          plot(Lag, ACF, type = "h", ylim = c(-1, 1), xaxt="n", ...)
          axis(1, at = lags.at)
        } else 
          plot(Lag, ACF, type = "h", ylim = c(-1, 1),  ...)
        abline(h = 0)
        
      } else if (graphs[j] == "pacf") {
        PACF <- as.numeric( pacorrC(mod$phi, mod$theta, lag.max) )
        if (!is.null(lags.at)) {
          plot(Lag, PACF, type = "h", ylim = c(-1, 1), xaxt="n", ...)
          axis(1, at = lags.at)
        } else 
          plot(Lag, PACF, type = "h", ylim = c(-1, 1),  ...)
        abline(h = 0)
      } else {
        A = spectrumC(mod$phi, mod$theta, mod$sig2, n.freq)
        if (log.spec) {
          plot(A[, 1], log(A[, 2]), type = "l", ylab = "log-SPEC",
               xlab = "Freq.", ...)
          abline(h = 0)
        }
        else
          plot(A[, 1], A[, 2], type = "l", ylab = "Spec", xlab = "Freq", ...)
      }
    }
  }
  
  invisible(NULL)
}

#' @rdname display
#' @export
display.default <- function(um, ...) {
  if( all(sapply(um, is.um)) ) {
    display.um(um, ...)
  } else {
    stop("Invalid um object")
  }
}

#' Equation of a univariate ARIMA model
#'
#' \code{equation} prints the equation of an object of class um.
#'
#' @param x an object of class \code{um}, a list of these objects or an object
#'   of class \code{ucarima}.
#' @param ... additional arguments.
#'
#' @examples
#' equation(um(ar = "(1 - 0.8B)"))
#'
#' @export
equation <- function (x, ...) UseMethod("equation")

#' @rdname equation
#' @param unscramble logical. If TRUE, AR, I and MA polynomials are unscrambled.
#' @param digits integer. Number of significant digits.
#' @param z character. Symbol for time series.
#' @param a character. Symbol for error.
#' @param width integer. Maximum width for line wrapping. If NULL, uses console width.
#' @export
equation.um <- function(x, unscramble = FALSE, digits = 4, 
                        z = "z", a = "a", width = NULL, ...) {
  
  txt <- ""
  if (unscramble) {
    if (x$p > 0) 
      txt <- paste0(txt, "(", as.character.lagpol(x$phi, ...), ")")
    if (x$d > 0) 
      txt <- paste0(txt, "(", as.character.lagpol(x$nabla, ...), ")")
    if (x$bc) txt <- paste0(txt, "log(", z, "_t) = ")
    else txt <- paste0(txt, z, "_t = ")
    if (x$q > 0) 
      txt <- paste0(txt, "(", as.character.lagpol(x$theta, ...), ")")
  } else {
    if (x$p > 0) {
      for (i in 1:x$kar) {
        txt <- paste0(txt, "(", as.character.lagpol(x$ar[[i]]$pol, ...), ")")
        if (x$ar[[i]]$p > 1)
          txt <- paste0(txt, "^", x$ar[[i]]$p)
      }
    }
    if (x$d > 0) {
      for (i in 1:length(x$i)) {
        txt <- paste0(txt, "(", as.character.lagpol(x$i[[i]]$pol, ...), ")")
        if (x$i[[i]]$p > 1)
          txt <- paste0(txt, "^", x$i[[i]]$p)
      }
    }
    if (x$bc) txt <- paste0(txt, "log(", z, "_t) = ")
    else txt <- paste0(txt, z, "_t = ")
    if (x$q > 0) {
      for (i in 1:x$kma) {
        txt <- paste0(txt, "(", as.character.lagpol(x$ma[[i]]$pol, ...), ")")
        if (x$ma[[i]]$p > 1)
          txt <- paste0(txt, "^", x$ma[[i]]$p)
      }
    }
  }
  
  equation_txt <- paste0(txt, a, "_t, s2", a, " = ", signif(x$sig2, digits = digits))
  
  if (is.null(width)) {
    width <- getOption("width", 80)
  }
  
  lines <- strwrap(equation_txt, width = width)
  
  cat(paste(lines, collapse = "\n"), "\n")
  
  invisible(equation_txt)
}

#' Easter effect
#'
#' \code{easter} extends the ARIMA model \code{um} by including a regression
#' variable to capture the Easter effect.
#'
#' @param um an object of class \code{\link{um}}.
#' @param z a time series.
#' @param len a positive integer specifying the duration of the Easter.
#' @param easter.mon logical. If TRUE Easter Monday is also taken into account. 
#' @param n.ahead a positive integer to extend the sample period of the
#'   Easter regression variable with \code{n.ahead} observations, which could be
#'   necessary to forecast the output.
#' @param envir environment in which the function arguments are evaluated.
#'    If NULL the calling environment of this function will be used.
#' @param ... other arguments.
#'
#' @return An object of class "\code{\link{tfm}}".
#' 
#' @examples
#' Y <- rsales
#' um1 <- um(Y, i = list(1, c(1, 12)), ma = list(1, c(1, 12)), bc = TRUE)
#' tfm1 <- easter(um1)
#' @export
easter <- function (um, ...) { UseMethod("easter") }

#' @rdname easter
#' @export
easter.um <- function(um, z = NULL, len = 4, easter.mon = FALSE, n.ahead = 0, 
                      envir = NULL, ...) {
  call <- match.call()
  if (is.null (envir)) envir <- parent.frame ()
  if (is.null(z)) 
    z <- z.um(um, envir = envir)
  else
    um$z <- deparse(substitute(z))
  
  if (frequency(z) != 12) stop("function only implemented for monthly ts")
  n.ahead <- abs(n.ahead)
  xreg <- matrix(EasterVar(z, len, easter.mon, n.ahead = n.ahead), ncol = 1)
  colnames(xreg) <- paste0("Easter", len, ifelse(easter.mon, "M", ""))
  tfm1 <- tfm(z, xreg = xreg, noise = um, envir = envir, new.name = FALSE)
  tfm1$call <- call
  tfm1
}

#' Factorized form of a univariate ARIMA model
#'
#' \code{factorize} .
#'
#' @param um an object of class \code{um} or a list of these objects.
#' @param ... additional arguments.
#' 
#' @examples
#' factorize(um(ar = "(1 - 0.8B)"))
#' 
#' @export 
factorize <- function (um, ...) UseMethod("factorize")

#' @rdname factorize
#' @param full logical value. If TRUE, lag polynomials are completely
#'   factorized. Otherwise, they are factored isolating positive real
#'   roots and grouping the remaining roots.
#' @export
factorize.um <- function(um, full = TRUE, ...) {
  if (!is.null(um$ar)) ar <- factors(as.lagpol(um$phi), full = full)
  else ar <- NULL
  if (!is.null(um$i)) i <- factors(as.lagpol(um$nabla), full = full)
  else i <- NULL
  if (!is.null(um$ma)) ma <- factors(as.lagpol(um$theta), full = full)
  else ma <- NULL
  um1 <- um(bc = um$bc, ar = ar, i = i, ma = ma, sig2 = um$sig2, ...)
  um1$z <- um$z
  return(um1)
} 

#'Estimation of the ARIMA model
#'
#'\code{fit} fits the univariate model to the time series z.
#'
#'@param mdl an object of class \code{\link{um}} or \code{\link{tfm}}.
#'@param method Exact/conditional maximum likelihood.
#'@param optim.method  the \code{method} argument of the \code{optim}
#' function.
#'@param show.iter logical value to show or hide the estimates at the
#' different iterations.
#'
#'@return An object of class "um" with the estimated parameters.
#'
#'@note
#' The \code{um} function estimates the corresponding ARIMA model when a time
#' series is provided. The \code{fit} function is useful to fit a model to
#' several time series, for example, in a Monte Carlo study.
#'
#' @export
fit <- function (mdl, ...) { UseMethod("fit") }

#' @rdname fit
#' @param z a time series.
#' @param envir environment in which the function arguments are evaluated.
#'    If NULL the calling environment of this function will be used.
#' @examples
#' z <- AirPassengers
#' airl <- um(i = list(1, c(1, 12)), ma = list(1, c(1, 12)), bc = TRUE)
#' airl <- fit(airl, z)
#' @export
fit.um <- function(mdl, z = NULL, method = c("exact", "cond"),
                   optim.method = "BFGS", show.iter = FALSE, envir=NULL, ... ) {

  stopifnot(inherits(mdl, "um"))
  if (is.null (envir)) envir <- parent.frame ()
  if (!is.stationary.um(mdl)) stop("Non-stationary AR preestimates")
  if (!is.invertible.um(mdl)) warning("Non-invertible MA preestimates")
  if ((mdl$p > 50 || mdl$q > 50) && length(method) > 1) method <- "cond"
  method <- match.arg(method)
  exact <- method == "exact"

  if (is.null(z)) {
    if (is.null(mdl$z)) stop("argment z required")
    else z <- eval(parse(text = mdl$z), envir)
  } else {
    mdl$z <- deparse(substitute(z))
  }
  
  logLik.arma <- function(b) {
    mdl <<- .update_um(mdl, b)
    if (mdl$is.adm) {
      if (is.null(mdl$mu)) {
        if (exact) ll <- ellarmaC(w, mdl$phi,mdl$theta)
        else ll <- cllarmaC(w, mdl$phi,mdl$theta)
      } else {
        if (exact) ll <- ellarmaC(w-mdl$mu, mdl$phi,mdl$theta)
        else ll <- cllarmaC(w-mdl$mu, mdl$phi,mdl$theta)
      }
    } else ll <- ll0
    
    if (show.iter) print(c(ll, b, mdl$is.adm ))
    
    mdl$is.adm <<- TRUE
    return(-ll)
  }
  if (mdl$bc & min(z) <= 0) mdl$bc <- FALSE
  w <- diffC(z, mdl$nabla, mdl$bc)
  b <- param.um(mdl)
  ll0 <- -logLik.arma(b)
  opt <- stats::optim(b, logLik.arma, method = optim.method, hessian = F)  
  if(opt$convergence > 0)
    warning(gettextf("possible convergence problem: optim gave code = %d",
                     opt$convergence), domain = NA)
  
  b <- opt$par
  mdl <- .update_um(mdl, b)

  n <- length(w)
  if (is.null(mdl$mu)) {
    if (exact) mdl$sig2 <- ssrC(w, mdl$phi, mdl$theta)/n
    else mdl$sig2 <- cssrC(w, mdl$phi, mdl$theta)/n
  } else {
    if(exact) mdl$sig2 <- ssrC(w-mdl$mu, mdl$phi, mdl$theta)/n
    else mdl$sig2 <- cssrC(w-mdl$mu, mdl$phi, mdl$theta)/n
  }
  
  mdl$optim <- opt
  mdl$method <- method
  mdl$b <- param.um(mdl)

  return(mdl)
  
}

#' Log-likelihood of an ARIMA model
#' 
#' \code{logLik} computes the exact or conditional log-likelihood of object of 
#' the class \code{um}.   
#' 
#' @param object an object of class \code{um}.
#' @param z an object of class \code{ts}.
#' @param method exact or conditional.
#' @param ... additional arguments.
#' 
#' @return 
#' The exact or conditional log-likelihood.
#' @export
logLik.um <-function(object, z = NULL, method = c("exact", "cond"), ...) {
  method <- match.arg(method)
  w <- w.um(object, z, TRUE)
  if (method == "exact") ll <- ellarmaC(w, object$phi, object$theta)
  else ll <- cllarmaC(w, object$phi, object$theta)
  
  return(ll)
}

#' @export
AIC.um <- function(object, z = NULL, method = c("exact", "cond"), ..., k = 2) {
  method <- match.arg(method)
  w <- w.um(object, z, TRUE)
  if (method == "exact") ll <- ellarmaC(w, object$phi, object$theta)
  else ll <- cllarmaC(w, object$phi, object$theta)
  b <- param.um(object)
  aic <- (-2.0*ll+k*length(b))/length(w)
}

#' Modifying a TF or an ARIMA model
#'
#' \code{modify} modifies an object of class \code{um} or \code{tfm} 
#' by adding and/or removing lag polynomials.
#' 
#' @param mdl an object of class \code{um} or \code{tfm}.
#' @inheritParams um
#' 
#' @return 
#' An object of class \code{um} or \code{um}.
#' 
#' @export
modify <- function (mdl, ...) { UseMethod("modify") }

#' @rdname modify
#' @examples
#' um1 <- um(ar = "(1 - 0.8B)")
#' um2 <- modify(um1, ar = list(0, "(1 - 0.9B)"), ma = "(1 - 0.5B)")
#' @export
modify.um <- 
function(mdl, ar = NULL, i = NULL, ma = NULL, mu = NULL, sig2 = NULL, 
         bc = NULL, ...) 
{
  stopifnot(is.um(mdl))
  
  if (is.null(mu)) mu <- mdl$mu
  else if (mu == 0) mu <- NULL
  if (is.null(sig2)) sig2 = mdl$sig2
  if (is.null(bc)) bc <- mdl$bc
  
  newop <- function(oldop, op) {
    if (is.null(op)) return(oldop)
    if (is.null(oldop)) return(op)
    if (is.lagpol(op)) return(c(oldop, list(op)))
    if (is.lagpol.list(op)) return(c(oldop, op))
    if (is.numeric(op)) op <- list(op)
    if (is.character(op)) op <- list(op)
    if (is.list(op)) {
      op1 <- lapply(op, function(x) {
        if (is.numeric(x)) {
          if (any(x <= 0)) x[x<=0]
          else NULL
        } else NULL
      })
      op2 <- lapply(op, function(x) {
        if (is.numeric(x)) {
          if (all(x <= 0)) NULL
          else x[!(x <= 0)]
        } else x
      })
      op1[sapply(op1, is.null)] <- NULL
      if (length(op1) != 0) {
        op1 <- unlist(op1)
        if (max(op1) == 0)
          oldop <- NULL
        else
          oldop <- oldop[op1]
      }
      oldop <- c(oldop, op2)
      oldop[sapply(oldop, is.null)] <- NULL
      return(oldop)
    }
    stop("invalid lag operator")
  }
  
  ar <- newop(mdl$ar, ar)
  i <- newop(mdl$i, i)
  ma <- newop(mdl$ma, ma)
  um(mdl$z, ar = ar, i = i, ma = ma, mu = mu, bc = bc, sig2 = mdl$sig2, ...)
  
}

#' Unscramble I polynomial
#' 
#' \code{nabla} multiplies the I polynomials of an object of 
#' the \code{um} class.
#'
#' @param x an object of class \code{um}.
#' @param ... additional arguments.
#' 
#' @return 
#' A numeric vector \code{c(1, a1, ..., ad)}
#' 
#' @note 
#' This function returns the member variable \code{um$nabla}.
#' 
#' @examples
#' um1 <- um(i = "(1 - B)(1 - B^12)")
#' nabla(um1)
#' @export
nabla <- function (x, ...) { UseMethod("nabla") }

#' @rdname nabla
#' @export
nabla.um <- function(x, ...) {
  stopifnot(inherits(x, "um"))
  as.lagpol(x$nabla)
}

#' Intervention analysis/Outlier treatment
#'
#' \code{intervention} estimates the effect of a intervention at a known time.
#'
#' @param mdl an object of class \code{\link{um}} or \code{\link{tfm}}.
#' @param y a "ts" object, optional.
#' @param type the type intervention (pulse, step, ramp) or the type of outlier
#'   (AO, LS, TC, IO).
#' @param time the date of the intervention, in format c(year, season).
#' @param n.ahead a positive integer to extend the sample period of the
#'   intervention variable with \code{n.ahead} observations, which could be
#'   necessary to forecast the output.
#' @param envir the environment in which to look for the time series z when it
#'   is passed as a character string.
#' @param ... additional arguments.
#' @return an object of class "\code{\link{tfm}}" or a table.
#'
#' @export
intervention <- function (mdl, ...) { UseMethod("intervention") }

#' @rdname intervention
#' @export
intervention.um <- function(mdl, y = NULL, type, time, n.ahead = 0, 
                             envir = parent.frame (), ...) {
  if (!is.null(y)) mdl$z <- deparse(substitute(y))
  y <- z.um(mdl, y, envir)
  tfm1 <- tfm(y, noise = mdl, fit = FALSE, new.name = FALSE, envir = envir)
  intervention.tfm(tfm1, NULL, type, time, n.ahead, envir, ...)
}
  
#' Outliers detection at known/unknown dates
#'
#' \code{outliers} performs a detection of four types of anomalies (AO, TC, LS
#' and IO) in a time series described by an ARIMA model. If the dates of the
#' outliers are unknown, an iterative detection process like that proposed by
#' Chen and Liu (1993) is conducted.
#'
#' @param mdl an object of class \code{\link{um}} or \code{\link{tfm}}.
#' @param y an object of class \code{ts}, optional.
#' @param types a vector with the initials of the outliers to be detected,
#'   c("AO", "LS", "TC", "IO").
#' @param dates a list of dates c(year, season). If \code{dates = NULL}, an
#'   iterative detection process is conducted.
#' @param c a positive constant to compare the z-ratio of the effect of an
#'   observation and decide whether or not it is an outlier. This argument is
#'   only used when \code{dates = NULL}.
#' @param calendar logical; if true, calendar effects are also estimated.
#' @param easter logical; if true, Easter effect is also estimated.
#' @param resid type of residuals (exact or conditional) used to identify
#'   outliers.
#' @param n.ahead a positive integer to extend the sample period of the
#'   intervention variables with \code{n.ahead} observations, which could be
#'   necessary to forecast the output.
#' @param p.value estimates with a p-value greater than p.value are omitted.
#' @param tc.fix a logical value indicating if the AR coefficient in the
#'   transfer function of the TC is estimated or fix.
#' @param envir environment in which the function arguments are evaluated.
#'    If NULL the calling environment of this function will be used.
#' @param ... other arguments.
#' @return an object of class "\code{\link{tfm}}" or a table.
#' @export
outliers <- function (mdl, ...) { UseMethod("outliers") }

#' @rdname outliers
#' @examples
#' Y <- rsales
#' um1 <- um(Y, i = list(1, c(1, 12)), ma = list(1, c(1, 12)), bc = TRUE)
#' outliers(um1)
#' @export
outliers.um <- function(mdl, y = NULL, types = c("AO", "LS", "TC", "IO"), 
                        dates = NULL, c = 3, calendar = FALSE, easter = FALSE, 
                        resid = c("exact", "cond"), n.ahead = 0, 
                        p.value = 1, tc.fix = TRUE, envir = NULL, ...) {
  if (is.null (envir)) envir <- parent.frame ()
  if (!is.null(y)) mdl$z <- deparse(substitute(y))
  y <- z.um(mdl, y, envir)
  tfm1 <- tfm(y, noise = mdl, fit = FALSE, new.name = FALSE, envir = envir)
  outliers.tfm(tfm1, NULL, types, dates, c, calendar, easter, resid, n.ahead, 
               p.value, tc.fix, envir, ...)
}


#' Unscramble AR polynomial
#' 
#' \code{phi} multiplies the AR polynomials of an object of 
#' the \code{um} class.
#'
#' @param x an object of class \code{um}.
#' @param ... additional arguments.
#' 
#' @return 
#' A numeric vector \code{c(1, a1, ..., ad)}
#' 
#' @note 
#' This function returns the member variable \code{um$phi}.
#' 
#' @examples
#' um1 <- um(ar = "(1 - 0.8B)(1 - 0.5B)")
#' phi(um1)
#' 
#' @export
phi <- function (x, ...) { UseMethod("phi") }

#' @rdname phi
#' @export
phi.um <- function(x, ...) {
  stopifnot(inherits(x, "um"))
  as.lagpol(x$phi)
}

#' Pi weights of an AR(I)MA model
#'
#' \code{pi.weights} computes the pi-weights of an AR(I)MA model.
#'
#' @param um an object of class \code{um}.
#' @param lag.max largest AR(Inf) coefficient required.
#' @param var.pi logical. If TRUE (FALSE), the I polynomials is 
#' considered (ignored).
#' @param ... additional arguments.
#' 
#' @return 
#' A numeric vector.
#' 
#' @examples
#' um1 <- um(i = "(1 - B)(1 - B^12)", ma = "(1 - 0.8B)(1 - 0.8B^12)")
#' pi.weights(um1, var.pi = TRUE)
#' 
#' @export
pi.weights <- function (um, ...) { UseMethod("pi.weights") }

#' @rdname pi.weights
#' @export
pi.weights.um <- function(um, lag.max = 10, var.pi = FALSE, ...) {
  stopifnot(inherits(um, "um"))
  if (var.pi) 
    piB <- as.numeric( polyratioC(polymultC(um$phi, um$nabla), 
                                  um$theta, lag.max) )
  else
    piB <- as.numeric( polyratioC(um$phi, um$theta, lag.max) )
  names(piB) <- paste("pi", 0:lag.max, sep = "")
  piB
}


#' Forecasts from an ARIMA model
#'
#' \code{predict} computes point and interval predictions for a time series 
#' from models of class  \code{\link{um}}.
#'
#' @param object an object of class \code{\link{um}}.
#' @param z an object of class \code{\link{ts}}.
#' @param ori the origin of prediction. By default, it is the last observation.
#' @param n.ahead number of steps ahead.
#' @param level confidence level.
#' @param i transformation of the series \code{z} to be forecasted. It is a
#' lagpol as those of a \code{\link{um}} object.  
#' @param envir environment in which the function arguments are evaluated.
#'    If NULL the calling environment of this function will be used.
#' @param ... additional arguments.
#' @return An object of class "\code{\link{predict.um}}".
#' @examples
#' Z <- AirPassengers
#' um1 <- um(Z, i = list(1, c(1, 12)), ma = list(1, c(1, 12)), bc = TRUE)
#' p <- predict(um1, n.ahead = 12)
#' p
#' plot(p, n.back = 60)
#'
#' @export predict.um
#' @export
predict.um <- function(object, z = NULL, ori = NULL, n.ahead = 1, level = 0.95,
                       i = NULL, envir=NULL, ...) {

  if (is.null (envir)) envir <- parent.frame ()
  if (is.null(z)) {
    if (is.null(object$z)) stop("argument z required")
    else {
      z <- eval(parse(text = object$z), envir)
    }
  }

  if (is.null(object$mu)) mu <- 0
  else mu <- object$mu
  
  if (!is.null(i)) {
    i <- lagpol0(i, "i", envir=envir)
    nabla <- polyexpand(i)
    nabla <- polydivC(object$nabla, nabla, FALSE, 1e-5)
    object$nabla <- nabla
    object$i <- list(as.lagpol(nabla))
    mu <- mu*sum(nabla)
    z <- ts(diffC(z, nabla, object$bc), end = end(z), frequency = frequency(z))
    if (object$bc){
      z <- z*100
      object$sig2 <- object$sig2*100^2
    }
    object$bc <- FALSE
  }
  
  if (is.null(ori)) ori <- length(z)
  if (n.ahead < 1) n.ahead <- 1
  
  X <-  forecastC(z, object$bc, mu, object$phi, object$nabla,
                  object$theta, object$sig2, ori, n.ahead)
  t <- (ori+1):(ori+n.ahead)
  z <- ts(X[, 1], start = start(z), frequency = frequency(z))
  
  dates <- time(zoo::as.zoo(z))
  f <- X[t, 1]
  if (any(level <= 0 || level >= 1)) level[level <= 0 || level >= 1] <- 0.95
  level <- unique(level)
  cv <- qnorm((1-level)/2, lower.tail = F)
  se <- sqrt(X[t, 4])
  if (object$bc) {
    low <- sapply(cv, function(x) f*exp(-x*se)) 
    upp <- sapply(cv, function(x) f*exp(x*se))
  } else {
    low <- sapply(cv, function(x) f - x*se) 
    upp <- sapply(cv, function(x) f + x*se)
  }
  
  out <- list(z = z, rmse = se, low = low, upp = upp, level = level,
              dates = dates, ori = ori, n.ahead = n.ahead, 
              ori.date = dates[ori])
  class(out) <- "predict.um"
  out  
}

#' Print univariate models
#' @rdname print
#' @param x An object of class um.
#' @param rows integer. Number of rows printed.
#' @param ... Additional arguments.
#' @export
print.predict.um <- function(x, rows = NULL, ...) {
  stopifnot(inherits(x, "predict.um"))
  if (is.null(rows)) rows <- 1:x$n.ahead
  t <- x$ori + rows
  df <- data.frame(x$z[t], x$rmse[rows])
  nms <- c("Forecast", "RMSE")
  k <- length(x$level)
  if (!is.matrix(x$low)) {
    x$low <- matrix(x$low, nrow = 1)
    x$upp <- matrix(x$upp, nrow = 1)
  }
  for (i in 1:k) {
    df1 <- data.frame(x$low[rows, i], x$upp[rows, i])
    nms <- c(nms, paste(x$level[i]*100, "% LB", sep = ""), 
             paste(x$level[i]*100, "% UB", sep = ""))
    df <- data.frame(df, df1)
  }
  colnames(df) <- nms
  rownames(df) <- x$dates[t]
  print(df)
}

#' @export
plot.predict.um <- function(x, n.back = 0, xlab = "Time", ylab = "z", main = "",
                            symbol= TRUE, pch = 16, col = c("black", "blue", "red"), ...) {
  n.back <- as.integer(n.back)
  n2 <- length(x$z)
  if (n.back > 0)  {
    if (n.back > x$ori) n.back <- x$ori
    n1 <- x$ori-n.back+1
  } else {
    n1 <- x$ori+1
  }
  z <- ts(x$z[n1:n2], end = end(x$z), frequency = frequency(x$z))
  zf <- ts(x$z[(x$ori+1):n2], end = end(x$z), frequency = frequency(x$z))
  zf.low <- ts(x$low, end = end(x$z), frequency = frequency(x$z)) 
  zf.upp <- ts(x$upp, end = end(x$z), frequency = frequency(x$z)) 
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  if (main != "")
    par(mar = c(3,3,3,0.8), mgp = c(1.5,0.6,0))
  else
    par(mar = c(3,3,1.5,0.8), mgp = c(1.5,0.6,0))
  
  if (length(col) != 3) col = c("black", "blue", "red")
  plot.ts(z, type = "n", ylim = c(min(z, x$low), max(z, x$upp)), 
          xlab = xlab, ylab = ylab, main = main)
  if (n.back > 0) { 
    z <- ts(x$z[n1:x$ori], end = x$ori.date, frequency = frequency(x$z))
    abline(h = mean(z), col = "gray")
    if (symbol) lines(z, type = "o", pch = 16, col = col[1])
    z <- ts(x$z[n1:(x$ori+1)], start = start(z), frequency = frequency(z))
    lines(z, col = col[1])
  }
  
  if (symbol) lines(zf, type = "o", pch = pch, cex = 0.6, col = col[2])
  else lines(zf, col = col[2])
  
  k <- ncol(x$low)
  for (i in 1:k) {
    lines(zf.low[, k-i+1], col = col[3], lty = i)
    lines(zf.upp[, k-i+1], col = col[3], lty = i)
  }
  
  invisible(NULL)
  
}

#' Print univariate models
#' @rdname print
#' @param arima logical. If TRUE, lag ARIMA polynomials are printed.
#' @export
print.um <- function(x, arima = FALSE, ...) {
  stopifnot(inherits(x, "um"))
  if (arima) {
    if (x$p > 0) {
      cat("AR polynomial:\n")
      print(phi(x))
    }
    if (x$d > 0) {
      cat("I polynomial:\n")
      print(nabla(x))
    }
    if (x$q > 0) {
      cat("MA polynomial:\n")
      print(theta(x))
    }
    print(x$sig2)
  } else {
    if (is.null(names(x$sig2))) names(x$sig2) <- "sig2"
    if (is.null(x$optim) || length(x$param) == 0)
      equation(x)
    else
      print(summary(x), short = TRUE, ...)
  } 
}


#' Psi weights of an AR(I)MA model
#'
#' \code{psi} computes the psi-weights of an AR(I)MA model.
#'
#' @param um an object of class \code{um}.
#' @param lag.max Largest MA(Inf) coefficient required.
#' @param var.psi logical. If TRUE the I polynomials is also inverted. If FALSE
#'   it is ignored.
#' @param ... additional arguments.   
#'
#' @return A numeric vector.
#'
#' @examples
#' um1 <- um(i = "(1 - B)(1 - B^12)", ma = "(1 - 0.8B)(1 - 0.8B^12)")
#' psi.weights(um1)
#' psi.weights(um1, var.psi = TRUE)
#' @export
psi.weights <- function (um, ...) { UseMethod("psi.weights") }

#' @rdname psi.weights
#' @export
psi.weights.um <- function(um, lag.max = 10, var.psi = FALSE, ...) {
  stopifnot(inherits(um, "um"))
  if (var.psi) 
    psi <- as.numeric( polyratioC(um$theta, polymultC(um$phi, um$nabla),
                                  lag.max) )
  else
    psi <- as.numeric( polyratioC(um$theta, um$phi, lag.max) )
  names(psi) <- paste("psi", 0:lag.max, sep = "")
  psi
}


#' Residuals of the ARIMA model
#'
#' \code{residuals} computes the exact or conditional residuals.
#'
#' @param object an object of class \code{um}.
#' @param z an object of class \code{ts}.
#' @param method exact/conditional residuals.
#' @param envir environment in which the function arguments are evaluated.
#'    If NULL the calling environment of this function will be used.
#' @param ... additional arguments.
#'
#' @return An object of class \code{um}.
#'
#' @examples
#' z <- AirPassengers
#' airl <- um(z, i = list(1, c(1, 12)), ma = list(1, c(1, 12)), bc = TRUE)
#' r <- residuals(airl)
#' summary(r)
#' @export
residuals.um <- function(object, z = NULL, method = c("exact", "cond"), envir=NULL, ...) {

  stopifnot(inherits(object, "um"))
  if (is.null (envir)) envir <- parent.frame ()
  if (is.null(z)) {
    if (is.null(object$z)) stop("argument z required")
    else {
      z <- eval(parse(text = object$z), envir)
    }
  }

  method <- match.arg(method)
  if (!is.null(object$method)) method <- object$method
  
  
  w <- diffC(z, object$nabla, object$bc)
  if (is.null(object$mu)) {
    if (method == "cond")
      res <- condresC(w, object$phi, object$theta, TRUE)
    else 
      res <- exactresC(w, object$phi, object$theta)
  } else {
    if (method == "cond")
      res <- condresC(w-object$mu, object$phi,object$theta, TRUE)
    else 
      res <- exactresC(w-object$mu, object$phi,object$theta)
  }
  
  if (is.ts(z))
    res <- ts(res[, 1], end = end(z), frequency = frequency(z))
  else
    res <- res[, 1]
  
  return(res)
  
}

#' @rdname roots
#' @param opr character. Operators for which roots are computed. Options: "arma",
#'   "arma", "ar", "ma", "i" or "arima".
#' @examples
#' um1 <- um(ar = "(1 - 0.8B)(1 - 0.8B^12)")
#' roots(um1)
#' @export
roots.um <- function(x, opr = c("arma", "ar", "ma", "i", "arima"), ...) {
  stopifnot(inherits(x, "um"))
  opr <- match.arg(opr)
  
  t <- list()
  if (startsWith(opr, "ar") & x$kar) {
    t <- lapply(x$ar, roots.lagpol)
  }
  
  if ( regexpr("i", opr)[1] > 0 & !is.null(x$i) ) {
    t <- c(t, lapply(x$i, roots.lagpol))
  }
  
  if (endsWith(opr, "ma") & x$kma) {
    t <- c(t, lapply(x$ma, roots.lagpol))
  }
  t
}


#' Simulate Time Series from ARIMA or Transfer Function Models
#' 
#' Generates random time series from ARIMA (\code{um}) or transfer function (\code{tfm}) models.
#' 
#' @param mdl An object of class \code{\link{um}} or \code{\link{tfm}}.
#' @param ... Additional arguments.
#' 
#' @return A \code{ts} object with the simulated time series.
#' 
#' @seealso \code{\link{sim.um}}, \code{\link{sim.tfm}}
#' @export
sim <- function (mdl, ...) { UseMethod("sim") }

#' @rdname sim
#' @param n Number of observations to simulate.
#' @param z0 Initial conditions for nonstationary series. Default is \code{NULL} (zero initial conditions).
#' @param n0 Number of initial observations to discard as burn-in. Default is \code{0}.
#' @param a Optional vector of innovations with length \code{n + n0}. If \code{NULL}, 
#'   innovations are drawn from \eqn{N(0, \sigma^2)}.
#' @param seed Random seed for reproducibility.
#' @param envir Environment for argument evaluation. Default is \code{parent.frame()}.
#'
#' @examples
#' # AR(1) model
#' mdl1 <- um(ar = "1 - 0.8B", sig2 = 1)
#' z1 <- sim(mdl1, n = 100, seed = 123)
#' 
#' # ARIMA(0,1,1) with burn-in
#' mdl2 <- um(i = 1, ma = "1 - 0.5B", sig2 = 1)
#' z2 <- sim(mdl2, n = 100, n0 = 50, seed = 456)
#'
#' @export
sim.um <- function(mdl, n = 100, z0 = NULL, n0 = 0, a = NULL, seed = NULL, 
                   envir = parent.frame(), ...) {
  stopifnot(inherits(mdl, "um"), n > length(mdl$nabla))
  if (is.null(mdl$mu)) mu <- 0
  else mu <- mdl$mu
  if (is.null(z0)) {
    if (is.null(mdl$z)) z0 <- 0
    else z0 <- eval(parse(text = mdl$z), envir)
  }
  if (!is.null(seed)) set.seed(seed)
  if (n0 < 0) n0 <- abs(n0)
  if (!is.null(a)) {
    stopifnot(length(a) == (n+n0))
  } else a <- rnorm(n + abs(n0), 0, sqrt(mdl$sig2))
  if (mdl$p > 0 || mdl$q > 0)
    a0 <- rnorm(max(c(mdl$p, mdl$q)), 0, sqrt(mdl$sig2))
  else
    a0 <- 0
  z <- simC(a, a0, mdl$bc, mu, mdl$phi, mdl$nabla, mdl$theta, z0)
  z <- as.numeric(z)
  if (n0 > 0) z <- z[-(1:n0)]
  if (is.ts(z0)) z <- ts(z, start = start(z0), frequency = frequency(z0)) 
  else z <- ts(z)
  return(z)
}  


#' Spectrum of an ARMA model
#'
#' \code{spec} computes the spectrum of an ARMA model.
#'
#' @param um an object of class \code{um}.
#' @param nabla logical. If TRUE, the pseudospectrum for a non stationary ARIMA
#'   model is calculated. By default, the spectrum is computed for the
#'   stationary ARMA model.
#' @param n.freq number of frequencies.
#' @param ... additional parameters.
#'
#' @return A matrix with the frequencies and the power spectral densities.
#'
#' @note The I polynomial is ignored.
#'
#' @examples
#' um1 <- um(i = "(1 - B)(1 - B^12)", ma = "(1 - 0.8B)(1 - 0.8B^12)")
#' s <- spec(um1, lag.max = 13)
#' @export
spec <- function (um, ...) { UseMethod("spec") }

#' @rdname spec
#' @export
spec.um <- function(um, nabla = FALSE, n.freq = 501, ...) {
  stopifnot(inherits(um, "um"))
  if (nabla)
    A <- spectrumC(polymultC(um$phi, um$nabla), um$theta, um$sig2, n.freq) 
  else
    A <- spectrumC(um$phi, um$theta, um$sig2, n.freq) 
  colnames(A) <- c("w", "fw")
  A
}

#' Summary of um model
#'
#' \code{summary} prints a summary of the estimation and diagnosis.
#'
#' @param object an object of class \code{um}.
#' @param z an object of class \code{ts}.
#' @param method exact/conditional maximum likelihood.
#' @param digits number of significant digits to use when printing.
#' @param envir environment in which the function arguments are evaluated.
#'    If NULL the calling environment of this function will be used.
#' @param ... additional arguments.
#'
#' @return
#' A list with the summary of the estimation and diagonosis.
#' 
#' @examples
#' z <- AirPassengers
#' airl <- um(z, i = list(1, c(1,12)), ma = list(1, c(1,12)), bc = TRUE)
#' summary(airl)
#'
#' @export
summary.um <- function(object, z = NULL, method = c("exact", "cond"),
                       digits = max(3L, getOption("digits") - 3L), envir=NULL, ...) {

  stopifnot(inherits(object, "um"))
  if (is.null (envir)) envir <- parent.frame ()
  model.name <- deparse(substitute(object))
  args <- list(...)
  if(is.null(args[["table"]])) table <- FALSE
  else table <- args[["table"]]
  
  if (is.null(z)) {
    if (is.null(object$z)) stop("argument z required")
    else {
      z <- eval(parse(text = object$z), envir)
    }
  }

  method <- match.arg(method)
  if (!is.null(object$method)) method <- object$method
  
  logLik.arma <- function(b) {
    object <<- .update_um(object, b)
    if (is.null(object$mu)) {
      if (method == "cond") ll <- cllarmaC(w, object$phi, object$theta)
      else ll <- ellarmaC(w, object$phi, object$theta)
    } else {
      if (method == "cond") ll <- ellarmaC(w-object$mu, object$phi, object$theta)
      else ll <- cllarmaC(w-object$mu, object$phi, object$theta)
    }
    ll
  }
  
  res.arma <- function(b) {
    object <<- .update_um(object, b)
    if (is.null(object$mu)) {
      if (method == "cond") res <- condresC(w, object$phi, object$theta, TRUE)
      else res <- gresC(w, object$phi, object$theta)
    } else {
      if (method == "cond") res <- condresC(w - object$mu, object$phi, object$theta, TRUE)
      else res <- gresC(w - object$mu, object$phi, object$theta)
    }
    as.vector(res)
  }
  
  w <- diffC(z, object$nabla, object$bc)
  N <- length(z)
  n <- length(w)
  b <- param.um(object)
  ll <- logLik.arma(b)
  aic <- (-2.0*ll+2*length(b))/n
  bic <- (-2*ll+log(n)*length(b))/n
  res <- res.arma(b)
  J <- numDeriv::jacobian(res.arma, b)
  g <- t(J) %*% res
  varb <- solve(t(J) %*% J)*(sum(res^2)/length(res))

  if (is.null(object$mu)) {
    if (method == "cond") {
      res <- as.numeric(condresC(w, object$phi, object$theta, TRUE))
      ssr <- sum(res^2)
    } else{
      res <- as.numeric(exactresC(w, object$phi, object$theta))
      ssr <- ssrC(w, object$phi, object$theta)
    }
  } else {
    if (method == "cond") {
      res <- as.numeric(condresC(w-object$mu, object$phi, object$theta, TRUE))
      ssr <- sum(res^2)
    } else{
      res <- exactresC(w-object$mu, object$phi, object$theta)
      ssr <- ssrC(w-object$mu, object$phi, object$theta)
    }
  }
  
  if (object$k > 0) b.names <- names(object$param)
  
  object$sig2 <- ssr/length(w)
  se <- sqrt(diag(varb))
  z.ratio <- b/se
  p <- pnorm(abs(z.ratio), lower.tail = F)*2
  
  X <- cbind(b, g, se, z.ratio, p)
  colnames(X) <- c("Estimate", "Gradient", "Std. Error", "z Value", "Pr(>|z|)")
  rownames(X) <- b.names
  if (table) return(X)
  
  tss <- var(z)*(N-1)
  mean.resid <- mean(res)
  rss <- var(res)*(length(res)-1)
  sd.resid <- sd(res)
  z.mean.resid <- mean.resid*sqrt(length(res))/sd.resid
  p.mean.resid <- pnorm(abs(z.mean.resid), lower.tail = F)*2
  if (is.ts(z)) {
    start <- start(z)
    end <- end(z)
    s <- frequency(z)
  } else {
    start <- NULL
    end <- NULL
    s <- NULL
  }
  
  Q1 <- Box.test(res, lag = object$p+object$q+1, type = "Ljung-Box", fitdf = object$p+object$q)
  Q2 <- Box.test(res, lag = as.integer(n/4)+object$p+object$q, type = "Ljung-Box",
                 fitdf = object$p+object$q)
  Q <- c(Q1 = Q1, Q2 = Q2)
  
  g <- rep(3, length(res))
  g[1:floor(length(res)/3)] <- 1
  g[(floor(length(res)/3)+1):floor(2*length(res)/3)] <- 2
  h.stat <- bartlett.test(res, g = g)  
  
  if (is.ts(z))
    res <- ts(res, end = end(z), frequency = frequency(z))
  out <- list(call = object$call, model.name = model.name, ml.method = object$method,
              z = object$z, aic = aic, bic = bic, gradient = g, varb = varb,
              logLik = ll, N = N, n = n, n.par = object$p+object$q,
              mean.resid = mean.resid, sd.resid = sd.resid, 
              z.mean.resid = z.mean.resid, p.mean.resid = p.mean.resid,
              resid = res, rss = ssr, sig2 = object$sig2, Q = Q, h.stat = h.stat,
              tss = tss, table = X, start = start, end = end, s = s)
  class(out) <- "summary.um"
  out
}

#' Print Summary of Univariate Model
#'
#' Print method for objects of class \code{summary.um}.
#'
#' @param x A \code{summary.um} object.
#' @param stats Logical. If TRUE, prints diagnostic statistics.
#' @param short Logical. If TRUE, prints abbreviated output.
#' @param digits Number of significant digits.
#' @param ... Additional arguments.
#'
#' @seealso \code{\link{summary.tfm}}
#'
#' @export
print.summary.um <- function(x, stats = TRUE, short = FALSE,
                             digits = max(3L, getOption("digits") - 3L), ...) {
  stopifnot(inherits(x, "summary.um")||inherits(x, "summary.tfm"), ...) 
  
  if (short) {
    print(x$table[, c(1, 3)])
    cat("\nlog likelihood: ", x$logLik)
    cat("\nResidual standard error: ", x$sd.resid)
    cat("\naic:", x$aic)
    return(invisible(NULL))
  }
  
  cat("\nModel:\n", x$model.name, " <- ", deparse(x$call, width.cutoff = 75L), "\n")
  cat("\nTime series:\n")
  if (!is.null(x$start)) {
    cat(" ", x$z, " (", paste(x$start, collapse = "."), " - ", 
        paste(x$end, collapse = "."), ")\n")
  } else {
    cat(x$z, "\n")
  }
  
  cat("\nMaximum likelihood method:\n", x$ml.method, "\n")
  
  cat("\nCoefficients:\n")
  printCoefmat(x$table, digits = digits, na.print = "NA")
  cat("\n")
  if (stats) {
    w.left <- 20
    w.right <- 10
    cat(format("Total nobs", width = w.left, justify = "left"),
        format(x$N, digits = NULL, width = w.right, justify = "right"), 
        format("Effective nobs", width = w.left, justify = "left"),
        format(x$n, digits = NULL, width = w.right, justify = "right"), "\n")
    
    cat(format("log likelihood", width = w.left, justify = "left"),
        format(x$logLik, digits = digits, width = w.right, justify = "right"),
        format("Error variance", width = w.left, justify = "left"),
        format(x$sig2, digits = digits, width = w.right, justify = "right"), "\n")
    
    cat(format("Mean of residuals", width = w.left, justify = "left"),
        format(x$mean.resid, digits = digits, width = w.right, justify = "right"),
        format("SD of the residuals", width = w.left, justify = "left"),
        format(x$sd.resid, digits = digits, width = w.right, 
               justify = "right"), "\n")
    
    label <- "z-test for residuals"
    cat(format(label, width = w.left, justify = "left"),
        format(x$z.mean.resid, digits = digits, width = w.right,
               justify = "right"),
        format("p-value", width = w.left, justify = "left"),
        format(x$p.mean.resid, digits = digits, width = w.right,
               justify = "right"), "\n")
    
    label <- paste("Ljung-Box Q(", x$Q$Q1.parameter, ") st.", sep = "")
    cat(format(label, width = w.left, justify = "left"),
        format(x$Q$Q1.statistic, digits = digits, width = w.right,
               justify = "right"),
        format("p-value", width = w.left, justify = "left"),
        format(x$Q$Q1.p.value, digits = digits, width = w.right,
               justify = "right"), "\n")
    
    label <- paste("Ljung-Box Q(", x$Q$Q2.parameter, ") st.", sep = "")
    cat(format(label, width = w.left, justify = "left"),
        format(x$Q$Q2.statistic, digits = digits, width = w.right,
               justify = "right"),
        format("p-value", width = w.left, justify = "left"),
        format(x$Q$Q2.p.value, digits = digits, width = w.right,
               justify = "right"), "\n")
    
    cat(format("Barlett H(3) stat.", width = w.left, justify = "left"),
        format(x$h.stat$statistic, digits = digits, width = w.right,
               justify = "right"),
        format("p-value", width = w.left, justify = "left"),
        format(x$h.stat$p.value, digits = digits, width = w.right, 
               justify = "right"), "\n")  
    
    cat(format("AIC", width = w.left, justify = "left"),
        format(x$aic, digits = digits, width = w.right, justify = "right"),
        format("BIC", width = w.left, justify = "left"),
        format(x$bic, digits = digits, width = w.right, justify = "right"), "\n")
  }
  
  invisible(NULL)
  
}  


#' Seasonal adjustment
#' 
#' \code{seasadj} removes the seasonal component of time series.
#' 
#' @inheritParams decomp
#' 
#' @return \code{seasadj} returns a seasonal adjusted time series.
#' 
#' @export
seasadj <- function (mdl, ...) { UseMethod("seasadj") }

#' @param envir environment in which the function arguments are evaluated.
#'    If NULL the calling environment of this function will be used.
#' @rdname seasadj
#' @examples
#' Y <- AirPassengers
#' um1 <- um(Y, bc = TRUE, i = list(1, c(1,12)), ma = list(1, c(1,12)))
#' Y <- seasadj(um1)
#' ide(Y)
#' @export
seasadj.um <-
function(mdl, z = NULL, method = c("mixed", "forecast", "backcast"), envir=NULL, ...)
{
  stopifnot(mdl$p+mdl$d > 1)
  if (is.null (envir)) envir <- parent.frame ()

  if (is.null(z)) {
    if (is.null(mdl$z)) stop("argument z required")
    else {
      z <- eval(parse(text = mdl$z), envir)
    }
  }
  method <- match.arg(method)
  if (is.null(mdl$mu)) mu <- 0
  else mu <- mdl$mu
  
  ariroots <- unlist( lapply( c(mdl$ar, mdl$i), function(x) roots(x, FALSE)) )
  if (method == "forecast") type = 1
  else if(method == "backcast") type = 2
  else type = 3
  end <- end(z)
  s <- frequency(z)
  z <- seasadjC(z, mdl$bc, mu, mdl$phi, mdl$nabla, mdl$theta, mdl$sig2,
                ariroots, type)
  z <- ts(z, end = end, frequency = s )
  
  return(z)
  
}

#' Wiener-Kolmogorov filter
#'
#' \code{wkfilter} extracts a signal for a time series described by an ARIMA
#' model given the ARIMA model for the signal.
#'
#' @param object an object of class \code{\link{um}}.
#' @param ... additional arguments.
#'
#' @return An object of class \code{ts} containing the estimated signal.
#'
#' @export
wkfilter <- function (object, ...) { UseMethod("wkfilter") }

#' @rdname wkfilter
#' @param um.uc ARIMA models for the observed time series and the
#'   unobserved component (signal).
#' @param z an optional \code{ts} object. If \code{NULL}, the time series to be
#'   filtered is contained in the \code{um.z} object.
#' @param output character, output of the function: `"series"` (default) returns
#'   the filtered time series, or `"filter"` returns the filter coefficients.
#' @param tol numeric tolerance used in polynomial divisions. Default is
#'   \code{1e-5}.
#' @param envir environment to get \code{z} when not provided.
#'
#' @examples
#' um1 <- um(AirPassengers, bc = TRUE, i = list(1, c(1,12)), ma = list(1, c(1,12)))
#' msx1 <- msx(um1, i = "(1-B)2")
#' trend <- wkfilter(um1, msx1$signal1)
#' seas <- wkfilter(um1, msx1$signal2)
#' @export
wkfilter.um <- function(object, um.uc, z = NULL, output = c("series", "filter"), 
                        tol = 1e-5, envir = parent.frame(), ...) {
  output <- match.arg(output)
  stopifnot(inherits(object, "um") && inherits(um.uc, "um"))
  stopifnot(object$sig2 > 0 && um.uc$sig2 > 0)
  if (!is.ts(z)) z <- z.um(object, z, envir)
  stopifnot(length(z) > object$p - object$d)
  
  .filter <- function(um, z, num, A, r) {
    n <- length(z)
    p <- predict(um, z = z, n.ahead = um$q + r)
    w <- if (um$bc) condresC(log(p$z), num, 1, FALSE) 
    else condresC(p$z, num, 1, FALSE)
    b <- c(w[(n+um$q-um$p-um$d+1):(n+um$q)], rep(0, um$q))
    x0 <- solve(A, b)
    if (um$q > 0) {
      x <- condres0C(w[1:(n+um$q-um$p-um$d)], 1, um$theta, 0, x0, FALSE)
      x <- c(x, x0)
    } else x <- w
    return(x[1:n])
  }  
  
  start <- start(z)
  s <- frequency(z)
  phi <- polydivC(object$phi, um.uc$phi, FALSE, tol)
  nabla <- polydivC(object$nabla, um.uc$nabla, FALSE, tol)
  num <- polymultC(um.uc$theta, polymultC(phi, nabla))
  r <- length(num) - 1
  if (object$q < r)
    th <- c(object$theta, rep(0, r - object$q))
  else {
    th <- object$theta
    r <- object$q
  }
  g <- tacovC(1, num, 1, r)
  A1 <- sapply(0:r, function(x) c(rep(0, x), th[1:(r+1-x)]))
  A2 <- sapply(r:0, function(x) c(rep(0, x), th[(r+1):(x+1)]))
  A <- A1 + A2
  num <- rev(as.numeric(solve(A, rev(g))))
  r <- object$p + object$d
  if (r < object$q) 
    r <- object$q
  
  A1 <- sapply(1:(object$p+object$d), function(i) {
    c(rep(0, i-1), object$theta, rep(0, object$p+object$d-i))
  })  
  A1 <- t(A1)
  if (object$q > 0) {
    phi <- polymultC(object$phi, object$nabla)  
    phi <- rev(as.numeric(phi))
    A2 <- sapply(1:object$q, function(i) {
      c(rep(0, i-1), phi, rep(0, object$q-i))
    })  
    A2 <- t(A2)
    A <- rbind(A1, A2)
  } else A <- A1
  
  if (output == "filter") {
    return(list(num = num, den = unname(object$theta)))
  }

  x1 <- .filter(object, z, num, A, r)
  x2 <- rev(.filter(object, rev(z), num, A, r))
  
  x1 <- (x1+x2)*(um.uc$sig2/object$sig2)
  x1 <- ts(x1, start = start, frequency = s)
  return(x1) 
}


#' Addition or substraction of univariate (ARIMA) models
#'
#' \code{add_um} creates a univariate (ARIMA) model from the addition or
#' substraction of two univariate (arima) models.
#'
#' @param um1,um2 Two "um" S3 objects.
#' @param add logical. If FALSE, the second model is substracted from the first
#'   one.
#' @param tol tolerance to check if a value is null.   
#' @return A "um" S3 object.
#' @note The + and - operators can also be used to add or substract ARIMA models.
#'
#' @examples
#' um1 <- um(i = "(1 - B)", ma = "(1 - 0.8B)")
#' um2 <- um(i = "(1 - B12)", ma = "(1 - 0.8B^12)")
#' um3 <- add_um(um1, um2)
#' um4 <- um3 - um2
#' @export
add_um <- function(um1, um2, add = TRUE, tol = 1.e-5) {
  stopifnot(is.um(um1) && is.um(um2))
  theta1 <- um1$theta; theta2 <- um2$theta
  if (add) {
    if (um1$d > 0 || um2$d > 0) {
      nabla <- common.factors(um1$nabla, um2$nabla, tol, TRUE) 
      nabla1 <- polydivC(um1$nabla, nabla, FALSE, tol)
      nabla2 <- polydivC(um2$nabla, nabla, FALSE, tol)
      theta1 <- polymultC(theta1, nabla2)
      theta2 <- polymultC(theta2, nabla1)
    } else nabla <- NULL
    if (um1$p > 0 || um2$p > 0) {
      phi <- common.factors(um1$phi, um2$phi, tol, TRUE) 
      phi1 <- polydivC(um1$phi, phi, FALSE, tol)
      phi2 <- polydivC(um2$phi, phi, FALSE, tol)
      theta1 <- polymultC(theta1, phi2)
      theta2 <- polymultC(theta2, phi1)
    } else phi <- NULL
    q <- max(c(length(theta1), length(theta2))) - 1
    g1 <- tacovC(1, theta1, um1$sig2, q)
    g2 <- tacovC(1, theta2, um2$sig2, q)
    g <- g1 + g2
    ma <- wold.pol(g, type = "c", tol = tol)
    if (um1$p > 0 || um2$p > 0)
      phi <- polymultC(polymultC(phi, phi1), phi2)
    if (um1$d > 0 || um2$d > 0)
      nabla <- polymultC(polymultC(nabla, nabla1), nabla2)
  } else {
    stopifnot(um1$d >= um2$d || um1$p >= um2$p)
    phi <- 1; nabla <- 1; theta <- 1
    if (um1$d > 0) {
      nabla <- polydivC(um1$nabla, um2$nabla, FALSE, tol)
      r <- polydivC(um1$nabla, um2$nabla, TRUE, tol)
      if (sum(abs(r)) > tol)
        stop("Incompatible models")
      theta2 <- polymultC(theta2, nabla)
    }
    if (um1$p > 0) {
      phi <- polydivC(um1$phi, um2$phi, FALSE, tol)
      r <- polydivC(um1$phi, um2$phi, TRUE, tol)
      if (sum(abs(r)) > tol)
        stop("Incompatible models")
      theta2 <- polymultC(theta2, phi)
    }
    q <- max(c(length(theta1), length(theta2))) - 1
    g1 <- tacovC(1, theta1, um1$sig2, q)
    g2 <- tacovC(1, theta2, um2$sig2, q)
    g <- g1 - g2
    if (g[1] < 0)
      stop("Incompatible models")
    ma <- wold.pol(g, type = "c", tol = tol)
    nabla <- um1$nabla
    phi <- um1$phi
  }
  oldw <- getOption("warn") 
  options(warn = -1) 
  on.exit(options(warn = oldw))
  sig2 <- ma[1]
  if (length(ma) > 1) ma <- ma[-1]
  else ma <- 1
  if (length(nabla) > 1 && length(ma) > 1) {
    lst <- common.factors(nabla, ma, tol)
    if (length(lst) > 0) {
      nabla1 <- polyexpand(lst)
      nabla <- polydivC(nabla, nabla1, FALSE, tol)
      ma <- polydivC(ma, nabla1, FALSE, tol)
    }
  }
  if (length(phi) > 1 && length(ma) > 1) {
    lst <- common.factors(phi, ma, tol)
    if (length(lst) > 0) {
      phi1 <- polyexpand(lst)
      phi <- polydivC(phi, phi1, FALSE, tol)
      ma <- polydivC(ma, phi1, FALSE, tol)
    } 
  } 
  if (abs(sig2) < tol^2) 
    um3 <- um(sig2 = sig2)
  else
    um3 <- um(ar = as.lagpol(phi, 1, "phi"), i = as.lagpol(nabla, 1, "nabla"),
            ma = as.lagpol(ma, 1, "theta"), sig2 = sig2)
  return(um3)
  
}

#' @export
Ops.um <- function(e1, e2) {
  if (nargs() == 1L) 
    stop(gettextf("unary %s not defined for \"um\" objects", 
                  .Generic), domain = NA)
  switch(.Generic,
              `+` = {
                add_um(e1, e2)  
              },
              `-` = {
                add_um(e1, e2, FALSE)
              },
              `*` = {
                if (is.numeric(e1)) {
                  e2$sig2 <- e2$sig2*e1[1]
                  e2
                } else if (is.numeric(e2)) {
                  e1$sig2 <- e1$sig2*e2[1]
                  e1
                } else stop(gettextf("%s not defined for \"um\" objects",
                                                   .Generic), domain = NA)
              },
              `/` = {
                if (is.numeric(e1)) {
                  e2$sig2 <- e2$sig2/e1[1]
                  e2
                } else if (is.numeric(e2)) {
                  e1$sig2 <- e1$sig2/e2[1]
                  e1
                } else stop(gettextf("%s not defined for \"um\" objects",
                                     .Generic), domain = NA)
              },
              stop(gettextf("%s not defined for \"um\" objects",
                            .Generic), domain = NA))
}

#' Structural form for an ARIMA model
#'
#' \code{sform} finds the structural form for an ARIMA model from its the
#' eventual forecast function.
#'
#' @param mdl an object of class \code{um}.
#'
#' @return An object of class \code{ucm}
#'
#' @export
sform <- function (mdl, ...) { UseMethod("sform") }

#' @rdname sform
#' @param z an optional time series.
#' @param msoe logical, TRUE for multiple source of errors and FALSE for single
#'   source of error.
#' @param cform logical. TRUE for contemporaneous form and FALSE for future
#'   form.
#' @param H an optional matrix to reduce the number of variances.
#' @param tol tolerance to check if the elements of b and C are zero.
#' @param nonadm character, the method to overcome nonadmissibility: non-linear 
#' least squares, quadratic programming or none.   
#' @param envir environment, see "\code{\link{um}}".
#' @param ... other arguments.
#'
#' @examples
#'
#' airl <- um(i = list(1, c(1, 12)), ma = "(1 - 0.8B)(1 - 0.8B12)")
#' sf <- sform(airl, index = c(1, 0, rep(2, 11)))
#' sf
#'
#' @export
#' 
sform.um <- function(mdl, z = NULL, msoe = TRUE, H = NULL, cform = TRUE, 
                     tol = 1.490116e-08, nonadm = c("quadprog", "nnls", "none"), 
                     envir = NULL, ...) {
  if (mdl$p + mdl$d + mdl$q == 0)
    stop("white noise process")
  nonadm <- match.arg(nonadm)
  if (is.null (envir)) envir <- parent.frame ()
  if (!is.null(z)) {
    mdl$z <- deparse(substitute(z))
    z <- eval(parse(text = z), envir)
  } else if (!is.null(mdl$z)) z <- eval(parse(text = mdl$z), envir)

  if (is.null(mdl$mu)) mu <- 0
  else mu <- mdl$mu
  ariroots <- unlist( lapply( c(mdl$ar, mdl$i), function(x) roots(x, FALSE)) )
  R <- sortrootsC(1/ariroots)
  C <- decompFC(R, mu)
  psi <- psi.weights(mdl, lag.max = ncol(C), var.psi = TRUE)[-1]
  C1 <- C[-nrow(C), ]
  C2 <- C[-1, ]
  c1 <- C[1, ]
  iC1 <- solve(C1)
  C <- iC1 %*% C2
  f <- iC1 %*% psi
  d <- as.vector(f)
  if (cform) {
    b <- as.vector(c1 %*% solve(C))
    kappa <- (1 - sum(b*d))
  } else {
    b <- c1
    kappa <- 1
  }
  b[abs(b) < tol] <- 0
  C[abs(C) < tol] <- 0
  
  if (msoe) {
    lf <- .LeverrierFaddeev(b, C)
    k <- length(b)
    g <- as.numeric( tacovC(1, mdl$theta, mdl$sig2, k) )
    g.i <- tacovC(1, lf$p, 1, k)
    A <- apply(lf$A, 2, function(x) {
      tacovC(1, x, 1, k)
    })
    A <- unname(cbind(g.i, A))
    if (is.null(H) && ncol(A) > 2) {
      i <- apply(A[, -1] - A[, -ncol(A)], 2, function(x) all(abs(x) < tol))
      if (any(i)) {
        H <- diag(ncol(A))
        j <- (2:ncol(A))[i]
        H[j, j-1] <- 1
        H <- H[, -j]
      }
    } 
    if (!is.null(H)) {
      s2 <- as.numeric(ginv(A %*% H) %*% g)
      s2 <- rowSums(H %*% diag(s2))
      S <- diag(s2) 
    } else {
      s2 <- as.numeric(ginv(A) %*% g)
      S <- diag(s2)
    }
    if (any(s2 < 0)) {
      warning("Nonadmissible structural form.")
      s2 <- sapply(s2, function(x) {
        if (x < 0 && abs(x)/max(s2) < tol) 0
        else x
      })
      if (nonadm == "nnls") {
        if (!is.null(H)) {
          x <- nnls::nnls(A %*% H, g)
          x <- rowSums(H %*% diag(x$x))
          S <- diag(x) 
        } else {
          x <- abs(nnls::nnls(A, g)$x)
          S <- diag(x)
        }
      } else if (nonadm == "quadprog") {
        if (!is.null(H)) A1 <- A %*% H
        else A1 <- A
        D <- t(A1) %*% A1
        d <- t(A1) %*% g
        x <- abs(quadprog::solve.QP(D, d, diag(1, ncol(A1)), rep(0, ncol(A1)))$solution)
        S <- diag(x)
        if (!is.null(H)) S <- diag(rowSums(H %*% S)) 
      } 
    }
    names(s2) <- paste0("s", 1:length(s2))
    A[abs(A) < tol] <- 0 
    ma <- cbind(lf$p, rbind(lf$A, rep(0,k)))
    ma[abs(ma) < tol] <- 0
    if (!is.null(H))
      aux <- list(C1 = C1, C2 = C2, ma = ma, A = A, H = H, b = g, x = s2)
    else
      aux <- list(C1 = C1, C2 = C2, ma = ma, A = A, b = g, x = s2)
    if (any(s2 < 0)) {
      if (nonadm == "nnls") aux$x.nnls <- x
      else if (nonadm == "quadprog") aux$x.quadprog <- x
    }
  } else {
    g <- c(kappa, f)
    S <- (g %*% t(g))*mdl$sig2
    aux <- list(C1 = C1, C2 = C2, g = g, psi = psi, sig2 = mdl$sig2)
  }
  S <- (S + t(S))/2
  S[abs(S) < .Machine$double.eps] <- 0
  mdl1 <- ssm(z, b, C, S, NULL, bc = mdl$bc, cform = cform)
  mdl1$z.name <- mdl$z
  mdl1$aux <- aux
  return(mdl1)
}

#' Unscramble MA polynomial
#'
#' @param um an object of class \code{um}.
#' 
#' @return 
#' A numeric vector \code{c(1, a1, ..., ad)}
#' 
#' @note 
#' This function returns the member variable \code{um$theta}.
#' 
#' @examples
#' um1 <- um(ma = "(1 - 0.8B)(1 - 0.5B)")
#' theta(um1)
#' 
#' @export
theta <- function (um) { UseMethod("theta") }

#' @rdname theta
#' @export
theta.um <- function(um) {
  stopifnot(inherits(um, "um"))
  as.lagpol(um$theta)
}

#  This function is based on the arima function of the stats package
#  of R. Below the copyright statement of the arima function is reproduced. 
#
#  File src/library/stats/R/arma0.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1999-2019 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

#' Diagnostic Plots for Time-Series Fits Description
#'
#' \code{tsdiag.um} is a wrap of the stats::tsdiag function.
#'
#' @param object a fitted \code{um} object.
#' @param gof.lag	the maximum number of lags for a Portmanteau goodness-of-fit test
#' @param ... additional arguments. 
#' 
#' @seealso stats::tsdiag.
#' 
#' @export
tsdiag.um <- function(object, gof.lag = 10, ...)
{
  ## plot standardized residuals, acf of residuals, Ljung-Box p-values
  oldpar <- par(mfrow = c(3, 1))
  on.exit(par(oldpar))
  rs <- resid(object)
  stdres <- rs/sqrt(object$sig2)
  plot(stdres, type = "h", main = "Standardized Residuals", ylab = "")
  abline(h = 0)
  acf(rs, plot = TRUE, main = "ACF of Residuals",
      na.action = na.pass)
  nlag <- gof.lag
  pval <- numeric(nlag)
  for(i in 1L:nlag) pval[i] <- Box.test(rs, i, type="Ljung-Box")$p.value
  plot(1L:nlag, pval, xlab = "lag", ylab = "p value", ylim = c(0,1),
       main = "p values for Ljung-Box statistic")
  abline(h = 0.05, lty = 2, col = "blue")
}

#' Unobserved components
#'
#' \code{decomp} estimates the unobserved components of a time series (trend,
#' seasonal, cycle, stationary and irregular) from the eventual forecast
#' function.
#'
#' @param mdl an object of class \code{\link{um}} or \code{\link{tfm}}.
#' @param method forward/backward forecasts or a mixture of the two.
#' @param ... additional arguments.
#'
#' @return A matrix with the unobserved components.
#'
#' @examples
#' Z <- AirPassengers
#' um1 <- um(Z, i = list(1, c(1, 12)), ma = list(1, c(1, 12)), bc = TRUE)
#' #uc1 <- decomp(um1)
#' @export
#'
decomp <- function (mdl, ...) { UseMethod("decomp") }

#' @rdname decomp
#' @param z an object of class \code{\link{ts}}.
#' @param envir environment in which the function arguments are evaluated.
#'    If NULL the calling environment of this function will be used.
#' @export
decomp.um <- function(mdl, z = NULL, 
                     method = c("msx", "msx0", "mixed", "forecast", "backcast"),
                     envir = parent.frame(), ...) {
  
  stopifnot(mdl$p+mdl$d > 1)
  method <- match.arg(method)

  if (method == "msx")
    return(.uc_msx(mdl, z, TRUE, envir))
  else if (method == "msx0")
    return(.uc_msx(mdl, z, FALSE, envir))

  z <- z.um(mdl, z, envir)

  if (mdl$bc) z <- log(z)
  if (is.null(z)) {
    if (is.null(mdl$z)) stop("argument z required")
    else {
      z <- eval(parse(text = mdl$z), envir)
    }
  }
    
    if (is.null(mdl$mu)) mu <- 0
    else mu <- mdl$mu
    
    ariroots <- unlist( lapply( c(mdl$ar, mdl$i), function(x) roots(x, FALSE)) )
    R <- sortrootsC(1/ariroots)
    matH <- decompHC(R, mu)
    matF <- decompFC(R, mu)
    if (mdl$bc) z <- log(z)
    if (method == "forecast") {
      matB <-  deceffBC(z, FALSE, mu, mdl$phi, mdl$nabla,
                        mdl$theta, mdl$sig2, matF, 1)
    } else if (method == "backcast") {
      matB <-  deceffBC(z, FALSE, mu, mdl$phi, mdl$nabla,
                        mdl$theta, mdl$sig2, matF, 2)
    } else {
      matB <-  deceffBC(z, FALSE, mu, mdl$phi, mdl$nabla,
                        mdl$theta, mdl$sig2, matF, 3)
    }
    
    s <- frequency(z)
    f <- matF[1, ]
    m <- ncol(matB)
    irreg <- ts( matB[, m], end = end(z), frequency = frequency(z) )
    matB <- matB[, 1:(m-1)]
    if (sum(matH[1, ]) > 0 ) 
      trend <- ts( matB %*% (matH[1, ] * f), end = end(z), frequency = s )
    else
      trend <- NULL
    if (sum(matH[2, ]) > 0 )
      seas <- ts( matB %*% (matH[2, ] * f), end = end(z), frequency = s )
    else
      seas <- NULL
    if (sum(matH[3, ]) > 0 )
      exp <- ts( matB %*% (matH[3, ] * f), end = end(z), frequency = s )
    else
      exp <- NULL
    if (sum(matH[4, ]) > 0 )
      cycle <- ts( matB %*% (matH[4, ] * f), end = end(z), frequency = s )
    else
      cycle <- NULL
    uc <- list(z = z, trend = trend, seas = seas, 
               exp = exp, cycle = cycle, irreg = irreg, 
               F = matF, B = matB, f = f, H = matH, R = R)
    class(uc) <- "uc.um"
    
    return(uc)
    
  }

.uc_msx <- function(mdl, z = NULL, canonical = FALSE, envir = NULL) {
  stopifnot(is.um(mdl))
  z <- z.um(mdl, z, envir)
  if (mdl$bc) X <- cbind(series = log(z))
  else X <- cbind(series = z)
  if (mdl$p > 1) {
    msx1 <- msx(mdl, ar = mdl$phi, canonical = canonical)
    trans <- wkfilter(mdl, msx1$signal, z)
    X <- cbind(X, trans = trans)
  } 
  
  if (mdl$d > 1) {
    i <- factors(as.lagpol(mdl$nabla), full = FALSE)
    hasTrend <- abs(sum(i[[1]]$pol)) < 1e-5
    if (hasTrend) {
      msx2 <- msx(mdl, i = i[[1]], canonical = canonical)
      trend <- wkfilter(mdl, msx2$signal, z)
      X <- cbind(X, trend = trend)
      hasSeas <- length(i) > 1
    } else {
      trend <- NULL 
      hasSeas <- TRUE
    }
  } 
  
  if (hasSeas) {
    if (hasTrend) {
      msx3 <- msx(mdl, i = i[-1], canonical = canonical)
    } else {
      msx3 <- msx(mdl, i = i, canonical = canonical)
    }
    seas <- wkfilter(mdl, msx3$signal, z)
    X <- cbind(X, seas = seas)
  }
  irreg <- X[ ,1]
  for (i in 2:ncol(X))
    irreg <- irreg - X[, i]
  X <- cbind(X, irreg = irreg)
  
  return(X)
}

#' @export
plot.uc.um <- function(x, main = NULL, ...) {
  stopifnot(inherits(x, "ucc.um"))
  k <- 2
  if (!is.null(x$trend)) k <- k + 1
  if (!is.null(x$seas)) k <- k + 1
  if (!is.null(x$exp)) k <- k + 1
  if (!is.null(x$cycle)) k <- k + 1
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  par(mfrow = c(k, 1), oma = c(3.5, 0.5, 2.5, 0.5),
      mar = c(0,2.5,0,2.5), mgp = c(1.5,0.6,0))
  plot.ts(x$z, ylab = "Series", axes = F, frame.plot = TRUE)
  if(!is.null(main)) title(main, line = 1, outer = TRUE)
  axis(side = 2)
  axis(side = 4, labels = F)
  left <- FALSE
  if (!is.null(x$trend)) {
    plot.ts(x$trend, ylab = "Trend", axes = F, frame.plot = TRUE)
    axis(side = 2, labels = left)
    axis(side = 4, labels = !left)
    left <- !left
  }
  if (!is.null(x$seas)) {
    plot.ts(x$seas, ylab = "Seasonal", axes = F, frame.plot = TRUE)
    axis(side = 2, labels = left)
    axis(side = 4, labels = !left)
    left <- !left
  }
  if (!is.null(x$exp)) { 
    plot.ts(x$exp, ylab = "Exp", axes = F, frame.plot = TRUE)
    axis(side = 2, labels = left)
    axis(side = 4, labels = !left)
    left <- !left
  }
  if (!is.null(x$cycle)) { 
    plot.ts(x$cycle, ylab = "Cycle", axes = F, frame.plot = TRUE)
    axis(side = 2, labels = left)
    axis(side = 4, labels = !left)
    left <- !left
  }
  plot.ts(x$irreg, ylab = "Irregular", type = "h",  axes = F, 
          frame.plot = TRUE)
  axis(side = 2, labels = left)
  axis(side = 4, labels = !left)
  axis(side = 1)
  abline(h = 0)
  mtext("Time", side = 1, line = 2)
  invisible(NULL)
  
}

# Non-exported functions
#' @noRd
backcast.um <- function(um, z = NULL, n.back = NULL, envir=NULL) {

  if (is.null (envir)) envir <- parent.frame ()
  if (is.null(z)) {
    if (is.null(um$z)) stop("argument z required")
    else {
      z <- eval(parse(text = um$z), envir)
    }
  }

  if (is.null(n.back) || n.back < 1) stop("argument n.back required")
  
  if (!is.ts(z)) z <- ts(z)
  end <- end(z)
  s <- frequency(z)
  
  if (is.null(um$mu)) mu <- 0
  else mu <- um$mu
  p <- backcastC(z, um$bc, mu, um$phi, um$nabla, um$theta, um$sig2, 1, n.back)  
  ts(as.numeric(p), end = end, frequency = s)
}

#' @noRd
eq.um <-function(um, digits = 2, arima = TRUE) {
  stopifnot(is.um(um))
  
  eq <- "'"  
  if (!is.null(um$ar)) {
    for (i in 1:um$kar) {
      pol <- paste("(", as.character(um$ar[[i]], digits, TRUE, TRUE), ")", sep = "")
      if (um$ar[[i]]$p > 1) pol <- paste(pol, "'^", um$ar[[i]]$p, "*'", sep = "")
      eq <- paste(eq, pol, sep = "")      
    }
  }
  
  if (arima && !is.null(um$i)) {
    for (i in 1:length(um$i)) {
      pol <- paste("(", as.character(um$i[[i]], digits, TRUE, TRUE), ")", sep = "")
      if (um$i[[i]]$p > 1) pol <- paste(pol, "'^", um$i[[i]]$p, "*'", sep = "")
      eq <- paste(eq, pol, sep = "")      
    }
    eq <- paste(eq, "z'[t]*' = ", sep = "")
  } else eq <- paste(eq, "w'[t]*' = ", sep = "")
  if (!is.null(um$ma)) {
    for (i in 1:um$kma) {
      pol <- paste("(", as.character(um$ma[[i]], digits, TRUE, TRUE), ")", sep = "")
      if (um$ma[[i]]$p > 1) pol <- paste(pol, "'^", um$ma[[i]]$p, "*'", sep = "")
      eq <- paste(eq, pol, sep = "")      
    }
  }
  eq <- paste(eq, "a'[t]", sep = "")
  eq
}

#' @noRd
hessian.um <- function(um, z = NULL, method = c("exact", "cond"), envir=NULL) {
  stopifnot(inherits(um, "um"))
  if (is.null (envir)) envir <- parent.frame ()
  if (!is.stationary.um(um)) stop("Non-stationary AR preestimates")
  #if (!is.invertible.um(um)) stop("Non-invertible MA preestimates")

  method <- match.arg(method)
  exact <- method == "exact"

  if (is.null(z)) {
    if (is.null(um$z)) stop("argment z required")
    else z <- eval(parse(text = um$z), envir)
  } else {
    um$z <- deparse(substitute(z))
  }
  
  logLik.arma <- function(b) {
    um <<- .update_um(um, b)
    if (is.null(um$mu)) {
      if (exact) ll <- ellarmaC(w, um$phi,um$theta)
      else ll <- cllarmaC(w, um$phi,um$theta)
    } else {
      if (exact) ll <- ellarmaC(w-um$mu, um$phi,um$theta)
      else ll <- cllarmaC(w-um$mu, um$phi,um$theta)
    }
    return(ll)
  }
  
  w <- diffC(z, um$nabla, um$bc)
  b <- param.um(um)
  H <- numDeriv::hessian(logLik.arma, b)
  H
}

#' @noRd
is.um <- function(um) {
  return(inherits(um, "um"))
}

#' @noRd
is.stationary.um <- function(um) {
  if (um$kar > 0) all(sapply(um$ar, admreg.lagpol))
  else return(TRUE)
}

#' @noRd
is.invertible.um <- function(um) {
  if (um$kma > 0) all(sapply(um$ma, admreg.lagpol))
  else return(TRUE)
}

#' @noRd
param.um <- function(um) {
  unlist(um$param, use.names = TRUE)
}

#' @export
plot.logLik.um <- function(x, z = NULL, par.name, support, 
                           method = "exact", ...) {
  w <- w.um(x, z)
  x$nabla <- 1
  x$bc <- FALSE
  b <- x$param[par.name]
  if (is.null(b)) stop("unknown parameter")
  vll <- sapply(support, function(par) {
    x$param[par.name] <- par
    x <- update(x, x$param)
    if (method == "exact") ll <- ellarmaC(w, x$phi, x$theta)
    else ll <- cllarmaC(w, x$phi, x$theta)
    ll
  })
  return(vll)
}

#' Update an object \code{um} 
#' 
#' \code{.update_um} updates an object of class \code{um} with a new vector 
#' of parameters
#' 
#' @param um an object of class \code{um}.
#' @param param a numeric vector.
#' 
#' @return 
#' An object of class \code{um}.
#' @noRd
.update_um <- function(um, param) {
  if (is.null(names(um$param))) return(um)
  um$param[] <- param[names(um$param)]
  um$is.adm <- TRUE
  
  if (!is.null(um$mu)) um$mu <- param["mu"]
  
  if (um$kar > 0) {
    um$ar <- lapply(um$ar, function(pol) {
      pol <- .update.lagpol(pol, param)
      ok <- admreg.lagpol(pol)
      if (!ok) um$is.adm <<- FALSE
      pol
    })
    um$phi <- polyexpand(um$ar)
  }
  
  if (um$kma > 0) {
    um$ma <- lapply(um$ma, function(pol) {
      pol <- .update.lagpol(pol, param) 
      ok <- admreg.lagpol(pol, FALSE)
      if (!ok) um$is.adm <<- FALSE
      pol
    })
    um$theta <- polyexpand(um$ma)
  }
  
  return(um)
  
}

#' Stationary time series for the ARIMA model
#'
#' \code{w} is the stationary time series for the ARIMA model.
#'
#' @param um an object of class \code{um}.
#' @param z an object of class \code{ts}.
#' @param wtilde logical. If TRUE the mean \code{mu} is subtracted from the
#'   stationary series.
#' @param envir environment in which the function arguments are evaluated.
#'    If NULL the calling environment of this function will be used.
#' @return an object of class \code{ts}.
#' @seealso \code{\link{z.um}}.
#' @examples
#' Z <- AirPassengers
#' um1 <- um(Z, bc = TRUE, i = list(1, c(1, 12)), ma = "(1 - 0.8B)(1 - 0.8B12)")
#' w(um1)
#' @noRd
w.um <- function(um, z = NULL, wtilde = FALSE, envir=NULL) {
  if (is.null (envir)) envir <- parent.frame ()
  z <- z.um(um, z, envir=envir)
  w <- as.vector(diffC(z, um$nabla, um$bc))
  if (wtilde & !is.null(um$mu)) w <- w - um$mu
  if (is.ts(z)) w <- ts(w, end = end(z), frequency = frequency(z))
  return(w)
}

#' Time series for the ARIMA model
#'
#' \code{z} is the time series for the ARIMA model.
#'
#' @param um an object of class \code{um}.
#' @param z an object of class \code{ts}.
#' @param envir environment in which the function arguments are evaluated.
#'    If NULL the calling environment of this function will be used.
#' @return an object of class \code{ts}.
#' @note This function is internally used by many methods of the \code{um} class
#'   to retrive the time series associated with the ARIMA model or to assign one
#'   new.
#' @examples
#' Z <- AirPassengers
#' um1 <- um(Z, bc = TRUE, i = list(1, c(1, 12)), ma = "(1 - 0.8B)(1 - 0.8B^12)")
#' z(um1)
#' z(um1, Z)
#' @noRd
z.um <- function(um, z = NULL, envir=NULL) {
  stopifnot(inherits(um, "um"))
  if (is.null (envir)) envir <- parent.frame ()
  if (is.null(z)) {
    if (is.null(um$z)) stop("argment z required")
    else {
      z <- eval(parse(text = um$z), envir)
    }
  }
  return(z)
}
