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
um <- function(z = NULL, ar = NULL, i = NULL, ma = NULL, mu = NULL,  
               sig2 = 1.0, bc = FALSE, fit = TRUE, ...) {

  call <- match.call()  
  if (is.numeric(z)){
    z <- deparse(substitute(z))
  } 
  
  if (!is.null(ar)) {
    ar <- .lagpol0(ar, "ar", "phi")
    phi <- polyexpand(ar)
  } else {
    phi <- 1.0
  }
  names(phi) <- paste("[phi", 0:(length(phi)-1), "]", sep = "")

  if (!is.null(i)) {
    i <- .lagpol0(i, "i","delta")
    nabla <- polyexpand(i)
  } else {
    nabla <- 1.0
  }
  
  if (!is.null(ma)) {
    ma <- .lagpol0(ma, "ma", "theta")
    theta <- polyexpand(ma)
  } else {
    theta <- 1.0
  }
  names(theta) <- paste("[theta", 0:(length(theta)-1), "]", sep = "")
  param <- unlist( lapply(c(ar, ma), function(x) unlist(x$param)) )
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
  mod <- update.um(mod, unlist(param, use.names = TRUE))
  if (mod$is.adm) {
    if (!is.null(z) & fit) mod <- fit.um(mod, ...)
    else mod$b <- param.um(mod)
  } else {
    warning("non-admisible model")
  }
  
  return(mod)
  
}


#' Convert \code{arima} into \code{um}.
#'
#' \code{as.um} converts an object of class \code{arima} into an object 
#' of class \code{um}.
#'
#' @param arima an object of class \code{arima}.
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
as.um <- function(arima) {
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
    mu <- -unname(arima$coef["intercept"])
    if (p+P > 1) mu <- mu/(1 - sum(arima$model$phi))
  } else mu <- NULL

  um <- um(ar = ar, i = i, ma = ma, mu = mu, sig2 = arima$sigma2)
  um$z <- arima$series
  um
}


#' Theoretical autocovariances of an ARMA model
#'
#' \code{autocov} computes the autocovariances of an ARMA model.
#'
#' @param um an object of class \code{um}.
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
autocov <- function (um, ...) { UseMethod("autocov") }

#' @rdname autocov
#' @export
autocov.um <- function(um, lag.max = 10, ...) {
  stopifnot(inherits(um, "um"))
  g <- as.numeric( tacovC(um$phi, um$theta, um$sig2, lag.max) )
  names(g) <- paste("gamma", 0:lag.max, sep = "")
  g
}


#' Theoretical simple/partial autocorrelations of an ARMA model
#'
#' \code{autocorr} computes the simple/partial autocorrelations of an ARMA model.
#'
#' @param um an object of class \code{um}.
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
autocorr <- function (um, ...) { UseMethod("autocorr") }

#' @rdname autocorr
#' @export
autocorr.um <- function(um, lag.max = 10, par = FALSE, ...) {
  stopifnot(inherits(um, "um"))
  if (!par) {
    g <- as.numeric( tacovC(um$phi, um$theta, um$sig2, lag.max) )
    names(g) <- paste("rho", 0:lag.max, sep = "")
    g <- g/g[1]
    g
  } else {
    g <- as.numeric( pacorrC(um$phi, um$theta, lag.max) )
    names(g) <- paste("phi", 1:lag.max, ".", 1:lag.max, sep = "")
    g
  }
}

#' Calendar effects
#'
#' \code{calendar} extends the ARIMA model \code{um} by estimating a transfer
#' function model with seven deterministic variables to capture the calendar
#' variation in a monthly time series. Two equivalent representations are
#' available: (1) D1, D2, ..., D7, (2) L, D1-D7, ..., D6-D7 where D1, D2, ...,
#' D7 are deterministic variables representing the number of Mondays, Tuesdays,
#' ..., Sundays, L = D1 + D2 + ... + D7 is the of the month. Optionally, a
#' deterministic variable to estimate the Easter effect can also be included.
#'
#' @param um an object of class \code{\link{um}}.
#' @param z a time series.
#' @param form representation for calendar effects: \code{form = td} form (1)
#'   above, \code{form = dif} form (1).
#' @param easter logical. If \code{TRUE} an Easter effect is also estimated.
#' @param n.ahead a positive integer to extend the sample period of the
#'   deterministic variables with \code{n.ahead} observations, which could be
#'   necessary to forecast the output.
#' @param p.value estimates with a p-value greater than p.value are omitted.   
#' @param ... additional arguments.
#' 
#' @return An object of class "\code{\link{tfm}}".
#' 
#' @references W. R. Bell & S. C. Hillmer (1983) Modeling Time Series with
#' Calendar Variation, Journal of the American Statistical Association, 78:383,
#' 526-534, DOI: 10.1080/01621459.1983.10478005
#'
#' @examples
#' Y <- tfarima::rsales
#' um1 <- um(Y, i = list(1, c(1, 12)), ma = list(1, c(1, 12)), bc = TRUE)
#' tfm1 <- calendar(um1)
#' 
#' @export
calendar <- function (um, ...) { UseMethod("calendar") }

#' @rdname calendar
#' @export
calendar.um <- 
function(um, z = NULL, form = c("dif", "td"), easter = FALSE, n.ahead = 0, 
         p.value = 1, ...)
{
  if (is.null(z)) z <- z.um(um)
  if (frequency(z) != 12) stop("function only implemented for monthly ts")
  form <- match.arg(form)
  
  n.ahead <- abs(n.ahead)
  xreg <- CalendarVar(z, form, easter, n.ahead)
  tfm1 <- tfm(xreg = xreg, noise = um)
  if (p.value < 0.999) {
    p <- summary.tfm(tfm1, p.values = TRUE)
    p <- (p[1:ncol(xreg)] <= p.value)
    if (all(p)) return(tfm1)
    if (any(p)) {
      xreg <- xreg[, p]
      tfm1 <- tfm(xreg = xreg, noise = tfm1$noise)
      return(tfm1)
    }
    return(um)
  }
  return(tfm1)
}

#' Diagnostic checking
#'
#' \code{diagchk} displays tools for diagnostic checking.
#'
#' @param mdl an object of class \code{um} or \code{tfm}.
#' @param method exact or conditional residuals.   
#' @param lag.max number of lags for ACF/PACF.
#' @param std logical. If TRUE standardized residuals are shown.
#' @param ... additional arguments.
#' 
#' @export
diagchk <- function (mdl, ...) { UseMethod("diagchk") }


#' @rdname diagchk
#' @param mdl an object of class \code{um}.
#' @param z optional, an object of class \code{ts}.
#' 
#' @examples
#' z <- AirPassengers
#' airl <- um(z, i = list(1, c(1,12)), ma = list(1, c(1,12)), bc = TRUE)
#' diagchk(airl)
#' @export
diagchk.um <- function(mdl, z = NULL, method = c("exact", "cond"), 
                       lag.max = NULL, std = TRUE, ...) {
  u <- residuals.um(mdl, z)
  ide(u, graphs = c("plot", "hist", "acf", "pacf", "cpgram"),
      lag.max = lag.max, std = std, ...)
}


#' Graphs for ARMA models
#'
#' \code{display} shows graphs characterizing one or a list of ARMA models.
#'
#' @param um an object of class \code{um} or a list of these objects.
#' @param lag.max number of lags for ACF/PACF.
#' @param n.freq number of frequencies for the spectrum.
#' @param log.spec logical. If TRUE log spectrum is computed.
#' @param graphs vector of graphs.
#' @param byrow orientation of the graphs.
#' @param eq logical. If TRUE the model equation is used as title.
#' @param ... additional arguments.
#' 
#' @examples
#' um1 <- um(ar = "(1 - 0.8B)(1 - 0.8B^12)")
#' um2 <- um(ma = "(1 - 0.8B)(1 - 0.8B^12)")
#' display(list(um1, um2))
#' 
#' @export 
display <- function (um, ...) { UseMethod("display") }

#' @rdname display
#' @export
display.um <- 
function(um, lag.max = 25, n.freq = 501, log.spec = FALSE,  
         graphs = c("acf", "pacf", "spec"), byrow = FALSE, eq = TRUE, ...) 
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
      text(.5, .5, parse(text = eq.um(um[[i]])), cex = 1.25,
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
        plot(Lag, ACF, type = "h", ylim = c(-1, 1), ...)
        abline(h = 0)
      } else if (graphs[j] == "pacf") {
        PACF <- as.numeric( pacorrC(mod$phi, mod$theta, lag.max) )
        plot(Lag, PACF, type = "h", ylim = c(-1, 1), ...)
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


#' Easter effect
#'
#' \code{easter} extends the ARIMA model \code{um} by including a regression
#' variable to capture the Easter effect.
#'
#' @param um an object of class \code{\link{um}}.
#' @param z a time series.
#' @param n.ahead a positive integer to extend the sample period of the
#'   Easter regression variable with \code{n.ahead} observations, which could be
#'   necessary to forecast the output.
#' @param ... additional arguments.
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
easter.um <- function(um, z = NULL, n.ahead = 0, ...) {
  if (is.null(z)) z <- z.um(um)
  if (frequency(z) != 12) stop("function only implemented for monthly ts")
  
  n.ahead <- abs(n.ahead)
  easter <- EasterVar(z, n.ahead = n.ahead)
  tfm(xreg = easter, noise = um)
  
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
#' @examples
#' z <- AirPassengers
#' airl <- um(i = list(1, c(1, 12)), ma = list(1, c(1, 12)), bc = TRUE)
#' airl <- fit(airl, z)
#' @export
fit.um <- function(mdl, z = NULL, method = c("exact", "cond"), 
                   optim.method = "BFGS", show.iter = FALSE, ... ) {
  
  stopifnot(inherits(mdl, "um"))
  if (!is.stationary.um(mdl)) stop("Non-stationary AR preestimates")
  if (!is.invertible.um(mdl)) stop("Non-invertible MA preestimates")
  if ((mdl$p > 50 || mdl$q > 50) && length(method) > 1) method <- "cond"
  method <- match.arg(method)
  exact <- method == "exact"
  
  if (is.null(z)) {
    if (is.null(mdl$z)) stop("argment z required")
    else z <- eval(parse(text = mdl$z), .GlobalEnv)
  } else {
    mdl$z <- deparse(substitute(z))
  }
  
  logLik.arma <- function(b) {
    mdl <<- update.um(mdl, b)
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
  
  w <- diffC(z, mdl$nabla, mdl$bc)
  b <- param.um(mdl)
  ll0 <- -logLik.arma(b)
  opt <- stats::optim(b, logLik.arma, method = optim.method, hessian = F)  
  if(opt$convergence > 0)
    warning(gettextf("possible convergence problem: optim gave code = %d",
                     opt$convergence), domain = NA)
  
  b <- opt$par
  mdl <- update.um(mdl, b)

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
         bc = NULL, fit = TRUE, ...) 
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
  um(mdl$z, ar = ar, i = i, ma = ma, mu = mu, bc = bc, sig2 = mdl$sig2)
  
}

#' Unscramble I polynomial
#' 
#' \code{nabla} multiplies the I polynomials of an object of 
#' the \code{um} class.
#'
#' @param um an object of class \code{um}.
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
nabla <- function (um) { UseMethod("nabla") }

#' @rdname nabla
#' @export
nabla.um <- function(um) {
  stopifnot(inherits(um, "um"))
  as.lagpol(um$nabla)
}

#' Outliers detection at known/unknown dates
#'
#' \code{outliers} performs a detection of four types of anomalies (AO, TC, LS
#' and IO) in a time series described by an ARIMA model. If the dates of the
#' outliers are unknown, an iterative detection process like that proposed by
#' Chen and Liu (1993) is conducted.
#'
#' @param mdl an object of class \code{\link{um}} or \code{\link{tfm}}.
#' @param z a time series.
#' @param dates a list of dates c(year, season). If \code{dates = NULL}, an
#'   iterative detection process is conducted.
#' @param c a positive constant to compare the z-ratio of the effect of an
#'   observation and decide whether or not it is an outlier. This argument is
#'   only used when \code{dates = NULL}.
#' @param calendar logical; if true, calendar effects are also estimated.
#' @param easter logical; if true, Easter effect is also estimated.
#' @param n.ahead a positive integer to extend the sample period of the
#'   intervation variables with \code{n.ahead} observations, which could be
#'   necessary to forecast the output.
#' @param p.value estimates with a p-value greater than p.value are omitted.   
#' @param ... additional arguments.   
#' @return an object of class "\code{\link{tfm}}" or a table.
#' @export
outliers <- function (mdl, ...) { UseMethod("outliers") }

#' @rdname outliers
#' @examples
#' Y <- rsales
#' um1 <- um(Y, i = list(1, c(1, 12)), ma = list(1, c(1, 12)), bc = TRUE)
#' outliers(um1)
#' @export
outliers.um <- function(mdl, z = NULL, dates = NULL, c = 3, calendar = FALSE, 
                        easter = FALSE, n.ahead = 0, p.value = 1, ...) {
  if (is.null(z)) z <- z.um(mdl)
  else if(!is.ts(z)) z <- ts(z)
  start <- start(z)
  freq <- frequency(z)
  if( freq < 2) {
    calendar <- FALSE
    easter <- FALSE
  }
  
  if (is.null(dates)) indx = 0  
  else if (is.numeric(dates)) {
    if (length(dates) == 2 && (dates[2] >= 1 && dates[2] <= freq)) {
      indx <- (dates[1] - start[1] + 1)*freq - (start[2] - 1) - (freq - dates[2])
    } 
    else {
      indx <- dates
    }
  } else if (is.list(dates)) { 
    indx <- sapply(dates, function(x) {
      if (freq > 1)
        (x[1] - start[1] + 1)*freq - (start[2] - 1) - (freq - x[2])
      else 
        (date[1] - start[1] + 1)*freq
    })
  } else stop("invalid stop argument")
  indx <- unique(indx)
  
  if (is.null(mdl$mu)) mu <- 0
  else mu <- mdl$mu
  
  tfm1 <- NULL
  if (calendar||easter) {
    if (calendar)
      tfm1 <- calendar.um(mdl, z, easter = easter, n.ahead = n.ahead)
    else
      tfm1 <- easter.um(mdl, z, n.ahead)
    N <- noise.tfm(tfm1, diff = FALSE)
    A <- outliersC(N, FALSE,  mu, tfm1$noise$phi, tfm1$noise$nabla, 
                   tfm1$noise$theta, indx, abs(c))
  } else { 
    A <- outliersC(z, mdl$bc,  mu, mdl$phi, mdl$nabla, mdl$theta, indx, abs(c))
  }
  if (ncol(A) == 1) {
    warning("no outlier")
    if (is.null(tfm1)) return(mdl)
    else return(tfm1)
  }
  
  df <- data.frame(obs = A[, 1], date = YearSeason(z, 2)[A[, 1]],
                   type = A[, 2], effect = A[, 3], t.ratio = A[, 4], 
                   w1 = A[, 5], t.w1 = A[, 6])
  df[[3]] <- factor(df[[3]], levels = 1:4, labels = c("IO", "AO", "LS", "TC"))
  
  n <- length(z)
  if (!is.null(n.ahead)) n <- n + n.ahead
  
  tfi <- lapply(1:nrow(df), function(i) {
    p <- double(n)
    p[df[i, 1]] <- 1
    p <- ts(p, start = start, frequency = freq)
    prefix <- paste(df[i, 3], df[i, 2], sep = "")
    if (df[i, 3] == "IO") 
      tf(p, w0 = df[i, 4], ma = mdl$ma, ar = c(mdl$i, mdl$ar), par.prefix = prefix)
    else if (df[i, 3] == "AO") 
      tf(p, w0 = df[i, 4], par.prefix = prefix)
    else if (df[i, 3] == "TC") 
      tf(p, w0 = df[i, 4], ar = "1 - 0.6B", par.prefix = prefix)
    else if (df[i, 3] == "LS") {
      p <- cumsum(p)
      p <- ts(p, start = start, frequency = freq)
      tf(p, w0 = df[i, 4], par.prefix = prefix)
    } else stop("unkown input")
  })
  
  if (calendar||easter) xreg <- tfm1$xreg
  else xreg <- NULL
  tfm1 <- tfm(xreg = xreg, inputs = tfi, noise = mdl)
  if (p.value < 0.999) 
    tfm1 <- varsel.tfm(tfm1, p.value = p.value)
  tfm1
}


#' Unscramble AR polynomial
#' 
#' \code{phi} multiplies the AR polynomials of an object of 
#' the \code{um} class.
#'
#' @param um an object of class \code{um}.
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
phi <- function (um) { UseMethod("phi") }

#' @rdname phi
#' @export
phi.um <- function(um) {
  stopifnot(inherits(um, "um"))
  as.lagpol(um$phi)
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
#' @param ... additional arguments.
#' @return An object of class "\code{\link{tfm}}".
#' @examples
#' Z <- AirPassengers
#' um1 <- um(Z, i = list(1, c(1, 12)), ma = list(1, c(1, 12)), bc = TRUE)
#' p <- predict(um1, n.ahead = 12)
#' p
#' plot(p, n.back = 60)
#' 
#' @export predict.um
#' @export
predict.um <- function(object, z = NULL, ori = NULL, n.ahead = 1, level = 0.95, ...) {
  
  if (is.null(z)) {
    if (is.null(object$z)) stop("argument z required")
    else {
      z <- eval(parse(text = object$z), .GlobalEnv)
    } 
  }
  
  if (is.null(object$mu)) mu <- 0
  else mu <- object$mu
  
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
                            symbol= TRUE, col = c("black", "blue", "red"), ...) {
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
  
  if (symbol) lines(zf, type = "o", pch = 16, cex = 0.6, col = col[2])
  else lines(zf, col = col[2])
  
  k <- ncol(x$low)
  for (i in 1:k) {
    lines(zf.low[, k-i+1], col = col[3], lty = i)
    lines(zf.upp[, k-i+1], col = col[3], lty = i)
  }
  
  invisible(NULL)
  
}


#' @export
print.um <- function(x, ...) {
  stopifnot(inherits(x, "um"))
  if (is.null(names(x$sig2))) names(x$sig2) <- "sig2"
  print( c(unlist(x$param), x$sig2) )
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
residuals.um <- function(object, z = NULL, method = c("exact", "cond"), ...) {

  stopifnot(inherits(object, "um"))
  if (is.null(z)) {
    if (is.null(object$z)) stop("argument z required")
    else {
      z <- eval(parse(text = object$z), .GlobalEnv)
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


#' Roots of the lag polynomials of an ARIMA model
#'
#' \code{roots} compute the roots of the AR, I, MA lag polynomials an ARIMA
#' model.
#'
#' @param x an object of class \code{um}.
#' @param opr character that indicates which operators are selected.
#' @param ... additional arguments.
#' 
#' @return 
#' List of matrices with the roots of each single polynomial.
#' 
#' @export
roots <- function (x, ...) { UseMethod("roots") }

#' @rdname roots
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


#' Time series simulation form an ARIMA or TF model
#' 
#' \code{sim} generates a random time series from an object of class \code{um} or \code{tfm}.   
#' 
#' @param mdl an object of class \code{um} or \code{tfm}.
#' @param n number of observations.
#' @param seed an integer.
#' @param ... additional arguments.
#' @return 
#' An object of class \code{ts}.
#' @export
sim <- function (mdl, ...) { UseMethod("sim") }

#' @rdname sim
#' @param z0 initial conditions for the nonstationary series.
#' @export
sim.um <- function(mdl, n = 100, z0 = NULL, seed = NULL, ...) {
  stopifnot(inherits(mdl, "um"), n > length(mdl$nabla))
  if (is.null(mdl$mu)) mu <- 0
  else mu <- mdl$mu
  if (is.null(z0)) {
    if (is.null(mdl$z)) z0 <- 0
    else z0 <- eval(parse(text = mdl$z), .GlobalEnv)
  }
  if (!is.null(seed)) set.seed(seed)
  a <- rnorm(n, 0, sqrt(mdl$sig2))
  z <- simC(a, mdl$bc, mu, mdl$phi, mdl$nabla, mdl$theta, z0)
  z <- as.numeric(z)
  if (is.ts(z0)) z <- ts(z, start = start(z0), frequency = frequency(z0))  
  return(z)
}  


#' Spectrum of an ARMA model
#'
#' \code{spec} computes the spectrum of an ARMA model.
#'
#' @param um an object of class \code{um}.
#' @param n.freq number of frequencies.
#' @param ... additional parameters.
#' 
#' @return 
#' A matrix with the frequencies and the power spectral densities.
#' 
#' @note 
#' The I polynomial is ignored.
#' 
#' @examples
#' um1 <- um(i = "(1 - B)(1 - B^12)", ma = "(1 - 0.8B)(1 - 0.8B^12)")
#' s <- spec(um1, lag.max = 13)
#' @export
spec <- function (um, ...) { UseMethod("spec") }

#' @rdname spec
#' @export
spec.um <- function(um, n.freq = 501, ...) {
  stopifnot(inherits(um, "um"))
  A <- spectrumC(um$phi, um$theta, um$sig2, n.freq) 
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
                       digits = max(3L, getOption("digits") - 3L), ...) {
  
  stopifnot(inherits(object, "um"))
  model.name <- deparse(substitute(object))
  
  if (is.null(z)) {
    if (is.null(object$z)) stop("argument z required")
    else {
      z <- eval(parse(text = object$z), .GlobalEnv)
    } 
  }
  
  method <- match.arg(method)
  if (!is.null(object$method)) method <- object$method
  
  logLik.arma <- function(b) {
    object <<- update.um(object, b)
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
    object <<- update.um(object, b)
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
  z <- b/se
  p <- pnorm(abs(z), lower.tail = F)*2
  
  X <- cbind(b, g, se, z, p)
  colnames(X) <- c("Estimate", "Gradient", "Std. Error", "z Value", "Pr(>|z|)")
  rownames(X) <- b.names
  
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


#' @export
print.summary.um <- function(x, stats = TRUE, 
                             digits = max(3L, getOption("digits") - 3L), ...) {
  stopifnot(inherits(x, "summary.um")||inherits(x, "summary.tfm"), ...) 
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
        format(x$N, digits = 0, width = w.right, justify = "right"), 
        format("Effective nobs", width = w.left, justify = "left"),
        format(x$n, digits = 0, width = w.right, justify = "right"), "\n")
    
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
#' @inheritParams ucomp
#' 
#' @return \code{seasadj} returns a seasonal adjusted time series.
#' 
#' @export
seasadj <- function (mdl, ...) { UseMethod("seasadj") }

#' @rdname seasadj
#' @examples
#' Y <- AirPassengers
#' um1 <- um(Y, bc = TRUE, i = list(1, c(1,12)), ma = list(1, c(1,12)))
#' Y <- seasadj(um1)
#' ide(Y)
#' @export
seasadj.um <- 
function(mdl, z = NULL, method = c("mixed", "forecast", "backcast"), ...) 
{
  stopifnot(mdl$p+mdl$d > 1)
  
  if (is.null(z)) {
    if (is.null(mdl$z)) stop("argument z required")
    else {
      z <- eval(parse(text = mdl$z), .GlobalEnv)
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


#' Sum of univariate (ARIMA) models
#' 
#' \code{sum_um} creates a univariate (ARIMA) model from the
#' sum of serveral univariate (arima) models.
#' 
#' @param ... List of "um" S3 objects.
#' @return A "um" S3 object.
#' @examples
#' um1 <- um(i = "(1 - B)", ma = "(1 - 0.8B)")
#' um2 <- um(i = "(1 - B12)", ma = "(1 - 0.8B^12)")
#' um3 <- sum_um(um1, um2)
#' @export
sum_um <- function(...) {
  models <- list(...)
  if (length(models) == 1) {
    if (is.um(models[[1]])) return(models[[1]])
    else models <- models[[1]]
  }
  
  if (!all(sapply(models, is.um))) {
    m <- list()
    l <- length(models)
    for (j in 1:l) {
      if (is.um(models[[j]]))
        m <- c(m, list(models[[j]]))
      else if(all(sapply(models[[j]], is.um))) {
        m <- c(m, models[[j]])
      } else
        stop("Arguments must be \"um\" objects")
    }
    l <- length(m)
    if (l < 1) stop("Arguments must be \"um\" objects")
    else if(l == 1) return(m)
    models <- m
  }
  
  group.factors <- function(list.lp, lp) {
    # Group AR or I operators of a "um" object by roots
    # E.g.: (1-B)*(1-B^s) = (1-B)^2*(1+B+...+B^s-1)
    n <- length(list.lp)
    for (j in 1:n) {
      lp1 <- list.lp[[j]]
      gcd <- polygcdC(lp1$pol, lp$pol)
      if (nrow(gcd) > 1) {
        q1 <- polydivC(lp1$pol, gcd, FALSE)
        q <- polydivC(lp$pol, gcd, FALSE)
        if (nrow(q1) == 1) {
          lp1$p <- lp1$p + lp$p
          lp1$Pol <- as.numeric(polyraiseC(lp1$pol, lp1$p))
          list.lp[[j]] <- lp1
          if (nrow(q) == 1) {
            lp <- NULL
            break
          } else lp <- as.lagpol(q, lp$p)
        } else {
          list.lp[[j]] <- as.lagpol(q1, p = lp1$p)
          if (nrow(q) == 1) {
            lp$p <- lp1$p + lp$p
            lp$Pol <- as.numeric(polyraiseC(lp$pol, lp1$p))
            break
          } else {
            list.lp <- c(list.lp, list(as.lagpol(gcd, lp1$p + lp$p)))
            lp <- as.lagpol(q, p = lp$p)
          }
        }
      }
    }
    if (!is.null(lp)) list.lp <- c(list.lp, list(lp))
    return(list.lp)
  }
  
  lcm.factors <- function(list.lp, lp) {
    # Least common multiple of several lag polynomials
    n <- length(list.lp)
    for (j in 1:n) {
      lp1 <- list.lp[[j]]
      gcd <- polygcdC(lp1$pol, lp$pol)
      if (nrow(gcd) > 1 && gcd[1] > 0.5 ) {
        q <- polydivC(lp$pol, gcd, FALSE)
        if (lp$p <= lp1$p) {
          if (nrow(q) > 1) {
            lp <- as.lagpol(q, p = lp$p)
          } else {
            lp <- NULL
            break
          }
        } else {
          q <- polymultC(polyraiseC(gcd, lp$p - lp1$p), polyraiseC(q, lp$p))
          lp <- as.lagpol(q, p = 1)
        }
      }
    }
    
    if (!is.null(lp)) list.lp <- c(list.lp, list(lp))
    
    return(list.lp)
  }
  
  ar <- lapply(models, function(um) {
    n <- length(um$ar)
    if (n > 1) {
      l <- um$ar[1]
      for (j in 2:n) 
        l <- group.factors(l, um$ar[[j]])
      l
    } else um$ar
  })
  
  ar <- unlist(ar, recursive = FALSE)
  n <- length(ar) 
  if (n == 0) list.ar <- list()
  else if (n == 1) list.ar <- ar 
  else {
    list.ar <- ar[1]
    for (j in 2:n) 
      list.ar <- lcm.factors(list.ar, ar[[j]])
  }
  
  # I factors
  i <- lapply(models, function(um) {
    n <- length(um$i)
    if (n > 1) {
      l <- um$i[1]
      for (j in 2:n) 
        l <- group.factors(l, um$i[[j]])
      l
    } else um$i
  })
  i <- unlist(i, recursive = FALSE)
  
  li <- length(i) 
  if (li == 0) list.i <- list()
  else if (li == 1) list.i <- i 
  else {
    list.i <- i[1]
    for (j in 2:li) 
      list.i <- lcm.factors(list.i, i[[j]])
  }
  
  # MA polynomials
  if (length(list.ar) > 0) phi <- polyexpand(list.ar)
  else {
    list.ar <- NULL
    phi <- 1
  }
  
  if (length(list.i) > 0) nabla <- polyexpand(list.i)
  else {
    list.i <- NULL
    nabla <- 1
  }
  
  nlags <- 0
  list.ma <- lapply(models, function(um) {
    phi1 <- polydivC(phi, um$phi, FALSE)
    nabla1 <- polydivC(nabla, um$nabla, FALSE)
    theta <- polymultC(phi1, nabla1)
    theta <- as.numeric( polymultC(theta, um$theta) )
    nlags <<- max(nlags, length(theta)-1 )
    theta*sqrt(um$sig2)
  })
  
  # MA autocovariances
  g <- sapply(list.ma, function(ma) {
    tacovC(1, ma, 1, nlags)
  })
  
  if (is.matrix(g)) {
    g <- rowSums(g)
    ma <- acovtomaC(g)
    sig2 <- ma[1]^2
    ma <- as.lagpol(c(1, -ma[-1]))
  } else {
    sig2 <- sum(g)
    ma <- NULL
  }
  
  um(ar = list.ar, i = list.i, ma = ma, sig2 = sig2)
  
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
#' \code{ucomp} estimates the unobserved components of a time series (trend,
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
#' uc <- ucomp(um1)
#' @export
#' 
ucomp <- function (mdl, ...) { UseMethod("ucomp") }

#' @rdname ucomp
#' @param z an object of class \code{\link{ts}}.
#' @export
ucomp.um <- 
  function(mdl, z = NULL, method = c("mixed", "forecast", "backcast"), ...) 
  {
    stopifnot(mdl$p+mdl$d > 1)
    
    if (is.null(z)) {
      if (is.null(mdl$z)) stop("argument z required")
      else {
        z <- eval(parse(text = mdl$z), .GlobalEnv)
      } 
    }
    method <- match.arg(method)
    
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
    class(uc) <- "ucomp.um"
    
    return(uc)
    
  }

#' @export
plot.ucomp.um <- function(x, main = NULL, ...) {
  stopifnot(inherits(x, "ucomp.um"))
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

backcast.um <- function(um, z = NULL, n.back = NULL) {
  
  if (is.null(z)) {
    if (is.null(um$z)) stop("argument z required")
    else {
      z <- eval(parse(text = um$z), .GlobalEnv)
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
    for (i in 1:um$ki) {
      pol <- paste("(", as.character(um$i[[i]], digits, TRUE, TRUE), ")", sep = "")
      if (um$i[[i]]$p > 1) pol <- paste(pol, "'^", um$i[[i]]$p, "*'", sep = "")
      eq <- paste(eq, pol, sep = "")      
    }
    eq <- paste(eq, "z'[t]' = ", sep = "")
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


hessian.um <- function(um, z = NULL, method = c("exact", "cond")) {
  stopifnot(inherits(um, "um"))
  if (!is.stationary.um(um)) stop("Non-stationary AR preestimates")
  #if (!is.invertible.um(um)) stop("Non-invertible MA preestimates")
  
  method <- match.arg(method)
  exact <- method == "exact"
  
  if (is.null(z)) {
    if (is.null(um$z)) stop("argment z required")
    else z <- eval(parse(text = um$z), .GlobalEnv)
  } else {
    um$z <- deparse(substitute(z))
  }
  
  logLik.arma <- function(b) {
    um <<- update.um(um, b)
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


is.um <- function(um) {
  return(inherits(um, "um"))
}

is.stationary.um <- function(um) {
  if (um$kar > 0) all(sapply(um$ar, admreg.lagpol))
  else return(TRUE)
}

is.invertible.um <- function(um) {
  if (um$kma > 0) all(sapply(um$ma, admreg.lagpol))
  else return(TRUE)
}


param.um <- function(um) {
  unlist(um$param, use.names = TRUE)
}


plot.logLik.um <- function(object, z = NULL, par.name, support, method = "exact") {
  w <- w.um(object, z)
  object$nabla <- 1
  object$bc <- FALSE
  b <- object$param[par.name]
  if (is.null(b)) stop("unknown parameter")
  vll <- sapply(support, function(x) {
    object$param[par.name] <- x
    object <- update(object, object$param)
    if (method == "exact") ll <- ellarmaC(w, object$phi, object$theta)
    else ll <- cllarmaC(w, object$phi, object$theta)
    ll
  })
  return(vll)
}


#' Update an object \code{um} 
#' 
#' \code{update.um} updates an object of class \code{um} with a new vector 
#' of parameters
#' 
#' @param um an object of class \code{um}.
#' @param param a numeric vector.
#' 
#' @return 
#' An object of class \code{um}.
#' @noRd
update.um <- function(um, param) {
  if (is.null(names(um$param))) return(um)
  um$param[] <- param[names(um$param)] 
  
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
#' @return an object of class \code{ts}.
#' @seealso \code{\link{z.um}}.
#' @examples
#' Z <- AirPassengers
#' um1 <- um(Z, bc = TRUE, i = list(1, c(1, 12)), ma = "(1 - 0.8B)(1 - 0.8B12)")
#' w(um1)
#' @noRd
w.um <- function(um, z = NULL, wtilde = FALSE) {
  z <- z.um(um, z)
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
z.um <- function(um, z = NULL) {
  stopifnot(inherits(um, "um"))
  if (is.null(z)) {
    if (is.null(um$z)) stop("argment z required")
    else {
      z <- eval(parse(text = um$z), .GlobalEnv)
    } 
  }
  return(z)
}

