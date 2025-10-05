## tfarima/R/tf.R
## Jose L Gallego (UC)

#' Transfer function for input
#'
#' \code{tf} creates a rational transfer function for an input, V(B) = w0(1 -
#' w_1B - ... - w_qB^q)/(1-d_1B - ... - d_pB^p)B^dX_t. Note that in this
#' specification the constant term of the MA polynomial is factored out so that
#' both polynomials in the numerator and denominator are normalized and can be
#' specified with the \code{lagpol} function in the same way as the operators of
#' univariate models.
#' @param x Input time series. A ts object or numeric vector. If NULL, input
#'   should be provided inside the \code{um} object.
#' @param delay Integer. Number of periods to delay the input (d in the transfer
#'   function).
#' @param w0 Numeric. Constant term of the transfer function polynomial V(B).
#' @param ar Character string or list. Specification of autoregressive
#'   polynomials in the denominator.
#' @param ma Character string or list. Specification of moving average
#'   polynomials in the numerator.
#' @param um Univariate model object. Model for the stochastic component of the
#'   input series.
#' @param n.back Integer. Number of backcasts to compute for extending the input
#'   series backward.
#' @param par.prefix Character. Prefix for parameter names in transfer function.
#' @param envir Environment. Environment for evaluating function arguments. If
#'   NULL, uses the calling environment.
#'
#' @return An object of class "tf" containing the transfer function specification.
#' Key components include:
#' \describe{
#'   \item{x}{Input time series data (possibly extended with backcasts).}
#'   \item{delay}{Number of periods delay in the transfer function.}
#'   \item{w0}{Constant gain parameter.}
#'   \item{phi}{AR polynomial coefficients (denominator).}
#'   \item{theta}{MA polynomial coefficients (numerator).}
#'   \item{p, q}{Orders of AR and MA polynomials respectively.}
#'   \item{param}{Named list of all model parameters.}
#'   \item{um}{Univariate model for the input series.}
#' }
#'
#' @references
#'
#' Box, G.E., Jenkins, G.M., Reinsel, G.C. and Ljung, G.M. (2015) Time Series
#' Analysis: Forecasting and Control. John Wiley & Sons, Hoboken.
#'
#' Wei, W.W.S. (2006) Time Series Analysis Univariate and Multivariate Methods.
#' 2nd Edition, Addison Wesley, New York, 33-59.
#'
#' @seealso \code{\link{um}}.
#'
#' @examples
#'
#' x <- rep(0, 100)
#' x[50] <- 1
#' tfx <- tf(x, w0 = 0.8, ar = "(1 - 0.5B)(1 - 0.7B^12)")
#'
#' @export
tf <- function(x = NULL, delay = 0,  w0 = 0.01, ar = NULL, ma = NULL, um = NULL, 
               n.back = NULL, par.prefix = "", envir = parent.frame ()) {
  
  if (is.null(x)) {
    if (is.null(um)) stop("input x required")
    else if (is.null(um$z)) stop("input x required")
    else { x <- get(um$z, envir); x.name <- um$z }
  } else if (is.character(x)) {
    x.name <- x; x <- get(x, envir)
  } else {
    x.name <- deparse(substitute(x))
  }

  if (par.prefix == "") par.prefix <- x.name
  pos <- regexpr("\\$", par.prefix)[1] 
  if (pos > -1) par.prefix <- substr(par.prefix, pos+1, nchar(par.prefix))
  else if (regexpr("\\[", par.prefix)[1] > -1) stop("wrong prefix for parameters")
  else if (regexpr("\\,", par.prefix)[1] > -1) stop("wrong prefix for parameters")
  
  if (!is.ts(x)) {
    stopifnot(is.numeric(x))
    x <- as.ts(x)
  }
  
  if (!is.numeric(delay) || length(delay) != 1 || delay < 0 ||
      delay != as.integer(delay)) {
    stop("delay must be a non-negative integer")
  }
  if (!is.numeric(w0) || length(w0) != 1) {
    stop("w0 must be a single numeric value")
  }  

  if (is.numeric(ar)) {if (any(ar == 0)) ar <- NULL}
  if (!is.null(ar)) {
    if (is.null(names(ar))) names(ar) <- paste0(par.prefix, ".d")
    ar <- lagpol0(ar, "ar", envir=envir)
    phi <- polyexpand(ar)
  } else {
    phi <- 1.0
  }
  names(phi) <- paste("[phi", 0:(length(phi)-1), "]", sep = "")

  if (is.numeric(ma)) {if (any(ma == 0)) ma <- NULL}
  if (!is.null(ma)) {
    if (is.null(names(ma))) names(ma) <- paste0(par.prefix, ".w")
    ma <- lagpol0(ma, "ma", envir=envir)
    theta <- polyexpand(ma)
  } else {
    theta <- 1.0
  }
  theta <- theta*w0  
  names(theta) <- paste("[theta", delay:(delay+length(theta)-1), "]", sep = "")
  param <- lapply(c(ar, ma), function(x) unlist(x$param))
  param <- unname(param)
  param <- unlist(param)
  if ( is.null(names(w0)) ) names(w0) <- par.prefix
  param <- as.list(c(w0, param))
  param <- param[!duplicated(names(param))]
  w0.expr <- parse(text = names(w0))
  
  if (is.null(um)) um <- um()
  else {
    stopifnot(is.um(um))
    if (is.null(n.back)) 
      n.back <- max(c(frequency(x)*3, length(phi) - 1, length(theta) + delay - 1))
    else
      n.back <- abs(n.back)
    if (n.back > 0) x <- backcast.um(um, x, n.back, envir=envir)
  }

  tf.input <- list(x.name = x.name, x = x, delay = delay, w0 = w0, w0.expr = w0.expr,
                   phi = phi, theta = theta, ar = ar, ma = ma, kar = length(ar),
                   kma = length(ma), p = length(phi)-1, q = length(theta)-1,
                   um = um, param = param, par.prefix = par.prefix, 
                   t.start = 1, t.end = length(x), n.back = n.back, is.adim = TRUE)
  class(tf.input) <- "tf"
  return(tf.input)
  
}



#' Helper function to create a tf object
#'
#' \code{tfest} estimates the transfer function
#'  \deqn{V(B) = w_0 * (theta(B)/phi(B)) * B^d}
#' that relates the input X_t to the output Y_t.
#' 
#' @param y Output series, a ts object or numeric vector.
#' @param x Input series, a ts object or numeric vector.
#' @param delay Integer. Number of periods delay (default 0).
#' @param p Integer. Order of the AR polynomial (default 1).
#' @param q Integer. Order of the MA polynomial (default 2).
#' @param um.y Univariate model for output series, um object or NULL.
#' @param um.x Univariate model for input series, um object or NULL.
#' @param n.back Integer. Number of backcasts to compute.
#' @param par.prefix Character. Prefix for parameter names.
#' @param envir Environment for evaluating arguments. If NULL, uses calling
#'   environment.
#'
#' @return An object of class \code{\link{tf}} containing preestimated transfer
#'   function parameters.
#'
#' @details Uses prewhitening to estimate initial parameter values.
#'
#' @examples
#'
#' Y <- seriesJ$Y - mean(seriesJ$Y)
#' X <- seriesJ$X - mean(seriesJ$X)
#' umx <- um(X, ar = 3)
#' umy <- fit(umx, Y)
#' tfx <- tfest(Y, X, delay = 3, p = 2, q = 2, um.x = umx, um.y = umy)
#'
#' @export
tfest <- function(y, x, delay = 0, p = 1, q = 2, um.y = NULL, um.x = NULL,
                  n.back = NULL, par.prefix = "", 
                  envir = envir <- parent.frame ()) {
  
  stopifnot(p >= 0, q >= 0, delay >= 0)
  if (is.null(um.y)) stop("Univariate model for output required")
  x.name <- deparse(substitute(x))
  if (!is.ts(y)) y <- ts(y)
  if (!is.ts(x)) x <- ts(x)
  s <- frequency(y)
  if (s != frequency(x)) stop("incompatible series")
  x.copy <- x
  if (is.null(n.back)) n.back  <- length(x)/4
  if (par.prefix == "") par.prefix <- x.name
  tf.x <- tf(x, delay = delay, ar = p, ma = q, um = um.x, par.prefix = par.prefix)
  tf.x$x.name <- x.name
  if (is.null(um.y)) um.y <- um()
  else um.y$param <- NULL
  tfm.x <- tfm(y, inputs = tf.x, noise = um.y, envir = envir)
  return(tfm.x$inputs[[1]])
}

#' Impulse response function
#' 
#' Computes and plots the impulse response function (IRF) or step response 
#' function (SRF) of a transfer function.
#' 
#' @param tf An object of class "tf".
#' @param lag.max Integer. Maximum number of lags to compute (default 10).
#' @param cum Logical. If TRUE computes step response function (cumulative), 
#'   if FALSE computes impulse response function (default FALSE).
#' @param plot Logical. If TRUE creates a plot, if FALSE returns values only (default TRUE).
#' 
#' @return If \code{plot = FALSE}, a named numeric vector with IRF/SRF values.
#'   If \code{plot = TRUE}, creates a plot and returns nothing (invisibly).
#'   
#' @examples
#' # Create transfer function
#' x <- rep(0, 100); x[50] <- 1
#' tfx <- tf(x, w0 = 0.8, ar = "(1 - 0.5B)")
#' 
#' # Plot impulse response function
#' irf(tfx, lag.max = 15)
#' 
#' # Get step response values without plot
#' srf_values <- irf(tfx, lag.max = 10, cum = TRUE, plot = FALSE)
#' 
#' @export
irf <- function(tf, lag.max = 10, cum = FALSE, plot = TRUE) {
  stopifnot(is.tf(tf))
  stopifnot(is.numeric(lag.max), length(lag.max) == 1, lag.max > 0)
  stopifnot(is.logical(cum), is.logical(plot))
  psi <- as.numeric( polyratioC(tf$theta, tf$phi, lag.max - tf$delay) )
  if (tf$delay > 0) {
    psi <- c(rep(0, tf$delay), psi)
  }
  if (cum) {
    psi <- cumsum(psi)
    names(psi) <- paste("[NU", 0:lag.max, "]", sep = "")
    ylab <- "SRF"
  } else {
    names(psi) <- paste("[nu", 0:lag.max, "]", sep = "")
    ylab <- "IRF"
  }
  
  if (plot) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    max_psi <- max(abs(psi))
    if (max_psi == 0) max_psi <- 1
    par(mar = c(3.0, 3.0, 1.0, 1.0), mgp = c(1.5, 0.6, 0))
    plot(0:lag.max, psi, type = "h", ylim = c(-max_psi, max_psi), 
         xlab = "Lag", ylab = ylab)
    abline(h = 0, col = "gray")
    invisible(psi)
  } else {
    psi
  }
}



#' Output of a transfer function or a transfer function model
#' 
#' \code{output} filters the input using the transfer function.
#' 
#' @param object an object of the S3 class "tf".
#' @param ... additional arguments.
#' 
#' @return A "ts" object
#' @export
output <- function(object, ...) {
  UseMethod("output")
}

#' @rdname output
#' @export
output.tf <- function(object, ...) {
  stopifnot(is.tf(object))
  x <- filterC(object$x, object$theta, object$phi, object$delay)
  ts(x, end = end(object$x), frequency = frequency(object$x))
}

#' Predict transfer function
#' @param object Transfer function object.
#' @param n.ahead Periods to predict.
#' @param ... Unused.
#' @return Point forecast for input.
#' @export
predict.tf <- function(object, n.ahead, ...) {
  start = start(object$x)
  frequency = frequency(object$x)
  if (object$um$p + object$um$d + object$um$q > 0) {
    ori <- length(object$x)
    if(is.null(object$um$mu)) mu <- 0
    else mu <- object$um$mu
    X <-  forecastC(object$x, object$um$bc, mu, object$um$phi, object$um$nabla,
                    object$um$theta, object$um$sig2, ori, n.ahead)
    object$x <- ts(X[, 1], start = start(object$x), frequency = frequency(object$x))
  } 
  object
}

#' Print method for transfer function objects
#' @rdname print
#' @export
print.tf <- function(x, ...) {
  if (!is.tf(x)) {
    stop("Object must be of class 'tf'")
  }  
  cat("Input:", x$x.name, "\n")
  cat("Sample:", start(x$x), " - ", end(x$x), "\n")
  cat("Freq.:", frequency(x$x), "\n")
  cat("Delay:", x$delay, "\n")
  
  if ( !is.null(x$ar) ) {
    if (length(x$ar) > 1) {
      cat("AR polynomials:\n")
      print(x$ar)
    }
    cat("AR polynomial:\n")
    print.lagpol(x$phi)
  }
  
  cat("w0:", x$w0, "\n")
  if ( !is.null(x$ma) ) {
    cat("Normalized MA polynomials:\n")
    print(x$ma)
    cat("MA polynomial:\n")
    print.lagpol(x$theta)
  }
  invisible(x)
}

#' @rdname roots
#' @param opr character. Operators for which roots are computed. Options: "arma",
#'   "ar" or "ma".
#' @export   
roots.tf <- function(x, opr = c("arma", "ar", "ma"), ...) {
  stopifnot(is.tf(x))
  opr <- match.arg(opr)
  
  t <- list()
  if (startsWith(opr, "ar") & x$kar) {
    t <- lapply(x$ar, roots.lagpol)
  }
  
  if (endsWith(opr, "ma") & x$kma) {
    t <- c(t, lapply(x$ma, roots.lagpol))
  }
  t
}

# Helper functions for tf objects

is.tf <- function(tf) {
  return(inherits(tf, "tf"))
}

is.list.tf <- function(ltf) {
  if(!base::is.list(ltf)) return(FALSE)
  all( sapply(ltf, is.tf) )
}

is.stationary.tf <- function(tf) {
  if (tf$kar > 0) all(sapply(tf$ar, admreg.lagpol))
  else return(TRUE)
}

is.invertible.tf <- function(tf) {
  if (tf$kma > 0) all(sapply(tf$ma, admreg.lagpol))
  else return(TRUE)
}

has.um.tf <- function(tf) {
  stopifnot(is.tf(tf))
  (tf$um$p + tf$um$d + tf$um$q) > 0
}

list.tf <- function(...) {
  names <- as.list(substitute(list(...)))[-1L]
  names <- sapply(names, as.character)
  args <- list(...)
  l <- length(args)
  ltf <- list()
  for (i in 1:l) {
    if (is.tf(args[[i]])) ltf <- c(ltf, list(args[[i]]))
    else if(is.list.tf(args[[i]])) ltf <- c(ltf, args[[i]])
    else if(is.numeric(args[[i]])) ltf <- c(ltf, list(tf(args[[i]], par.prefix = names[i])))
    else stop( paste(class(args[[i]]), " no tf object", sep = ""))
  }
  ltf
}

param.tf <- function(tf) {
  unlist(tf$param, use.names = FALSE)
}

#' Update an object \code{tf} 
#' 
#' \code{.update_tf} updates an object of class \code{tf} with a new vector 
#' of parameters
#' 
#' @param tf an object of class \code{tf}.
#' @param param a numeric vector.
#' 
#' @return 
#' An object of class \code{tf}.
#' @noRd
.update_tf <- function(tf, param){
  tf$param[] <- param[names(tf$param)]
  tf$w0 <- eval(tf$w0.expr, envir = tf$param)
  if (tf$kar > 0) {
    tf$ar <- lapply(tf$ar, .update.lagpol, param = param) 
    tf$phi <- polyexpand(tf$ar)
  }
  
  if (tf$kma > 0) {
    tf$ma <- lapply(tf$ma, .update.lagpol, param = param) 
    tf$theta <- polyexpand(tf$ma)
  } else tf$theta <- 1
  
  tf$theta <- tf$theta*tf$w0  
  
  return(tf)
  
}

var.predict.tf <- function(tf, n.ahead = 10) {
  stopifnot(inherits(tf, "tf"))
  num <- polymultC( tf$theta, tf$um$theta )
  den <- polymultC(tf$phi, polymultC(tf$um$phi, tf$um$nabla))
  psi <- as.numeric( polyratioC(num, den, n.ahead - 1) )
  names(psi) <- paste("psi", 0:(n.ahead - 1), sep = "")
  cumsum(psi^2)*tf$um$sig2
}









#' @export
sim.tf <- function(mdl, N=100, x0=NULL, envir=NULL, ...) {
  if (is.null (envir)) envir <- parent.frame ()
  if (is.null(mdl$um)) {
    x <- mdl$x
  } else {
   x <- sim.um(mdl$um, N, x0, envir=envir)
  }
  output.tf(mdl)
}
