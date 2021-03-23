#' Transfer function for input
#'
#' \code{tf} creates a rational transfer function for an input, V(B) = w0(1 -
#' w_1B - ... - w_qB^q)/(1-d_1B - ... - d_pB^p)B^dX_t. Note that in this
#' specification the constant term of the MA polynomial is factored out so that
#' both polynomials in the numerator and denominator are normalized and can be
#' specified with the \code{lagpol} function in the same way as the operators of
#' univariate models.
#'
#' @param x input, a ts object or a numeric vector.
#' @param delay integer.
#' @param w0 constant term of the polynomial V(B), double.
#' @param ar list of stationary AR polynomials.
#' @param ma list of MA polynomials.
#' @param um univariate model for stochastic input. 
#' @param n.back number of backcasts to extend the input.
#' @param par.prefix prefix name for parameters.
#'
#' @return An object of the class "tf".
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
tf <- function(x = NULL, delay = 0,  w0 = 0, ar = NULL, ma = NULL, 
               um = NULL, n.back = NULL, par.prefix = "") {
  if (is.null(x)) {
    if (is.null(um)) stop("input x required")
    else if (is.null(um$z)) stop("input x required")
    else { x <- get(um$z, .GlobalEnv); x.name <- um$z }
  } else if (is.character(x)) { 
    x.name <- x; x <- get(x, .GlobalEnv) 
  } else {
    x.name <- deparse(substitute(x))
  }

  if (par.prefix == "") par.prefix <- x.name
  pos <- regexpr("\\$", par.prefix)[1] 
  if (pos > -1) par.prefix <- substr(par.prefix, pos+1, nchar(par.prefix))
  else if (regexpr("\\[", par.prefix)[1] > -1) stop("wrong preffix for parameters")
  else if (regexpr("\\,", par.prefix)[1] > -1) stop("wrong preffix for parameters")
  
  if (!is.ts(x)) {
    stopifnot(is.numeric(x))
    x <- ts(x)
  }
  
  if (!is.null(ar)) {
    ar <- .lagpol0(ar, "ar", paste(par.prefix, ".d", sep = ""))
    phi <- polyexpand(ar)
  } else {
    phi <- 1.0
  }
  names(phi) <- paste("[phi", 0:(length(phi)-1), "]", sep = "")

  if (!is.null(ma)) {
    ma <- .lagpol0(ma, "ma", paste(par.prefix, ".w", sep = ""))
    theta <- polyexpand(ma)
  } else {
    theta <- 1.0
  }
  theta <- theta*w0  
  names(theta) <- paste("[theta", delay:(delay+length(theta)-1), "]", sep = "")
  
  param <- unlist( lapply(c(ar, ma), function(x) unlist(x$param)) )
  if ( is.null(names(w0)) ) {
      names(w0) <- par.prefix
  }
  
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
    if (n.back > 0) x <- backcast.um(um, x, n.back)
  }
  
  tf.input <- list(x.name = x.name, x = x, delay = delay, w0 = w0, w0.expr = w0.expr,
                   phi = phi, theta = theta, ar = ar, ma = ma, kar = length(ar),
                   kma = length(ma), p = length(phi)-1, q = length(theta)-1,
                   um = um, param = param, par.prefix = par.prefix, 
                   t.start = 1, t.end = length(x), n.back = n.back, is.adim = TRUE)
  class(tf.input) <- "tf"
  return(tf.input)
  
}







is.tf <- function(tf) {
  return(inherits(tf, "tf"))
}

has.um.tf <- function(tf) {
  stopifnot(is.tf(tf))
  (tf$um$p + tf$um$d + tf$um$q) > 0
}

is.list.tf <- function(ltf) {
  if(!base::is.list(ltf)) return(FALSE)
  all( sapply(ltf, is.tf) )
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

is.stationary.tf <- function(tf) {
  if (tf$kar > 0) all(sapply(tf$ar, admreg.lagpol))
  else return(TRUE)
}

is.invertible.tf <- function(tf) {
  if (tf$kma > 0) all(sapply(tf$ma, admreg.lagpol))
  else return(TRUE)
}

update.tf <- function(tf, param){
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

print.tf <- function(tf) {
  cat("Input:", tf$x.name, "\n")
  cat("Sample:", start(tf$x), " - ", end(tf$x), "\n")
  cat("Freq.:", frequency(tf$x), "\n")
  cat("Delay:", tf$delay, "\n")
  
  if ( !is.null(tf$ar) ) {
    if (length(tf$ar) > 1) {
      cat("AR polynomials:\n")
      print(tf$ar)
    }
    cat("AR polynomial:\n")
    print.lagpol(tf$phi)
  }
  
  cat("w0:", tf$w0, "\n")
  
  if ( !is.null(tf$ma) ) {
    cat("Normalized MA polynomials:\n")
    print(tf$ma)
    cat("MA polynomial:\n")
    print.lagpol(tf$theta)
  }
  
}

irf.tf <- function(tf, lag.max = 10, cum = FALSE, plot = TRUE) {
  stopifnot(inherits(tf, "tf"))
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
    max <- max(abs(psi))
    par(mar = c(3.0, 3.0, 1.0, 1.0), mgp = c(1.5, 0.6, 0))
    plot(0:lag.max, psi, type = "h", ylim = c(-max, max), xlab = "Lag", ylab = ylab)
    abline(h = 0)
    invisible()
  } else {
    psi
  }
}

#' Output of a transfer function
#' 
#' \code{output} filters the input using the transfer function.
#' 
#' @param tf an object of the S3 class "tf".
#' 
#' @return A "ts" object
output.tf <- function(tf) {
  stopifnot(inherits(tf, "tf"))
  x <- filterC(tf$x, tf$theta, tf$phi, tf$delay)
  ts(x, end = end(tf$x), frequency = frequency(tf$x))
}

#' Preestimates of a transfer function
#' 
#' \code{tfest} provides preestimates of the transfer function
#' between an output and an input.
#' 
#' @param y output, a ts object or a numeric vector.
#' @param x input, a ts object or a numeric vector.
#' @param p order of the AR polynomial, double.
#' @param q order of the MA polynomial, double.
#' @param delay integer.
#' @param um.y univariate model for output, um object or NULL.
#' @param um.x univariate model for input, um object or NUL.
#' @param n.back number of backcasts.
#' @return A "tf" S3 object
#' @export
tfest <- function(y, x, delay = 0, p = 1, q = 2, um.y = NULL, um.x = NULL, 
                  n.back = NULL) {
  
  stopifnot(p >= 0, q >= 0, delay >= 0)
  if (is.null(um.y)) stop("Univariate model for output required")
  x.name <- deparse(substitute(x))
  if (!is.ts(y)) y <- ts(y)
  if (!is.ts(x)) x <- ts(x)
  s <- frequency(y)
  if (s != frequency(x)) stop("incompatible series")
  x.copy <- x
  if (is.null(n.back)) n.back  <- length(x)/4
  tf.x <- tf(x, delay = delay, ar = p, ma = q, um = um.x, par.prefix = x.name)
  if (is.null(um.y)) um.y <- um()
  else um.y$param <- NULL
  tfm.x <- tfm(deparse(substitute(y)), inputs = tf.x, noise = um.y)
  return(tfm.x$inputs[[1]])

}


#' Prewhitened cross correlation function
#'
#' \code{pccf} displays cross correlation function between input and output
#' after prewhitening both through a univariate model.
#'
#' @param x input, a 'ts' object or a numeric vector.
#' @param y output, a 'ts' object or a numeric vector.
#' @param um.x univariate model for input.
#' @param um.y univariate model for output.
#' @param lag.max number of lags, integer.
#' @param main title of the graph.
#' @param nu.weights logical. If TRUE the coefficients of the IRF are 
#' computed instead of the cross-correlations. 
#' @param plot logical value to indicate if the ccf graph must be graphed or
#'   computed.
#' @param ... additional arguments.   
#'
#' @return The estimated cross correlations are displayed in a graph or returned
#'   into a numeric vector.
#' @export
pccf <- function(x, y, um.x = NULL, um.y = NULL, lag.max = NULL, plot = TRUE, 
                 main = NULL, nu.weights = FALSE, ...) {
  xlab <- deparse(substitute(x))
  ylab <- deparse(substitute(y))
  
  if (!is.ts(y)) y <- ts(y)
  if (!is.ts(x)) x <- ts(x)
  if (!is.null(ncol(x))) x <- x[, 1]
  if (!is.null(ncol(y))) y <- y[, 1]

  s <- frequency(y)
  if (s != frequency(x)) stop("incompatible series")
  
  if (!is.null(um.x)) {
    x <- residuals.um(um.x, x, method = "cond")
    if (is.null(um.y)) y <- residuals.um(um.x, y, method = "cond")
    else y <- residuals.um(um.y, y, method = "cond")     
  } else if (!is.null(um.y)) {
    y <- residuals.um(um.y, y, method = "cond")
    x <- residuals.um(um.y, x, method = "cond")
  }
  
  x <- window(x, start = start(y), end = end(y))

  cc <- stats::ccf(x, y, lag.max = lag.max, plot = FALSE)
  lag.max <- (dim(cc$lag)[1]-1)/2

  if (plot) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    if (is.null(main))
      main <- substitute(rho*"("*xlab[t+k]*","*ylab[t]*")",
                         list(xlab = xlab, ylab = ylab))
    if (main !=  "")
      par(mar = c(3.0, 3.0, 3.5, 1.0), mgp = c(1.5, 0.6, 0))
    else
      par(mar = c(3.0, 3.0, 1.0, 1.0), mgp = c(1.5, 0.6, 0))
    # stats::ccf(x, y, lag.max = lag.max, ylab = "CCF", ci.col = "gray",
    #            main = main, plot = plot)
    # abline(v = 0, col = "gray", lty = 2)
    cc$lag[,1,1] <- (-lag.max) : lag.max
    plot(cc, main = main, ylab = "CCF", ci.col = "gray")
    invisible()
  } else {
    cc <- stats::ccf(x, y, lag.max = lag.max, plot = FALSE)
    if (nu.weights) cc <- cc*sd(y)/sd(x)
    return(cc)
  }
}

roots.tf <- function(tf, opr = c("arma", "ar", "ma")) {
  stopifnot(inherits(tf, "tf"))
  opr <- match.arg(opr)
  
  t <- list()
  if (startsWith(opr, "ar") & tf$kar) {
    t <- lapply(tf$ar, roots.lagpol)
  }
  
  if (endsWith(opr, "ma") & tf$kma) {
    t <- c(t, lapply(tf$ma, roots.lagpol))
  }
  t
}

sim.tf <- function(tf, N=100, x0=NULL) {
  if (is.null(tf$um)) {
    x <- tf$x
  } else {
   x <- sim.um(tf$um, N, x0) 
  }
  output.tf(tf)
}

predict.tf <- function(tf, n.ahead) {
  start = start(tf$x)
  frequency = frequency(tf$x)
  if (tf$um$p + tf$um$d + tf$um$q > 0) {
    ori <- length(tf$x)
    if(is.null(tf$um$mu)) mu <- 0
    else mu <- tf$um$mu
    X <-  forecastC(tf$x, tf$um$bc, mu, tf$um$phi, tf$um$nabla,
                    tf$um$theta, tf$um$sig2, ori, n.ahead)
    tf$x <- ts(X[, 1], start = start(tf$x), frequency = frequency(tf$x))
  } 
  tf
}

var.predict.tf <- function(tf, n.ahead = 10) {
  stopifnot(inherits(tf, "tf"))
  
  num <- polymultC( tf$theta, tf$um$theta )
  den <- polymultC(tf$phi, polymultC(tf$um$phi, tf$um$nabla))
  psi <- as.numeric( polyratioC(num, den, n.ahead - 1) )
  names(psi) <- paste("psi", 0:(n.ahead - 1), sep = "")
  cumsum(psi^2)*tf$um$sig2
}
