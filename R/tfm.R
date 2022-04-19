#' Transfer function models
#'
#' \code{tfm} creates a multiple input transfer function model.
#'
#' @param output a ts object or a numeric vector.
#' @param xreg a matrix of regressors.
#' @param inputs a list of tf objects.
#' @param noise a um object for the noise.
#' @param fit logical. If TRUE, model is fitted.
#' @param envir environment in which the function arguments are evaluated. If
#'   NULL the calling environment of this function will be used.
#' @param new.name logical. Argument used internally: if TRUE a new name is
#'   assigned to the output, otherwise it keeps its name saved in noise$z.
#' @param ... additional arguments.
#'
#' @return An object of the class \code{tfm}.
#'
#' @seealso \code{\link{tf}} and \code{\link{um}}.
#'
#' @references
#'
#' Box, G.E., Jenkins, G.M., Reinsel, G.C. and Ljung, G.M. (2015) Time Series
#' Analysis: Forecasting and Control. John Wiley & Sons, Hoboken.
#'
#' @export
tfm <- function(output = NULL, xreg = NULL, inputs = NULL, noise, fit = TRUE,
                envir = NULL, new.name = TRUE,  ...) {
  call <- match.call()
  stopifnot(is.um(noise))
  if (is.null (envir)) envir <- parent.frame ()

  if (!is.null(inputs)) {
    if (inherits(inputs, "tf")) {
      inputs <- list(inputs)
    } else {
      inputs <- inputs[!sapply(inputs, is.null)]
    }
    k <- length(inputs)
  } else k <- 0
  
  if (is.null(output)) {
    if (is.null(noise$z)) stop("missing output")
    else output <- eval(parse(text = noise$z), envir)
  } else if (is.character(output)) {
    noise$z <- output
    output <- eval(parse(text = noise$z), envir)
  } else if (is.numeric(output)) {
    if (new.name)
      noise$z <- deparse(substitute(output))
  } else stop("invalid output")

  N <- length(output)
  stopifnot(N > noise$d)
  if (!is.ts(output)) output <- as.ts(output)
  start <- start(output)
  end <- end(output)
  s <- frequency(output)
  
  y <- diffC(output, noise$nabla, noise$bc)
  n <- length(y)
  
  param <- list()
  if (k > 0) {
    for (i in 1:k) {
      stopifnot(inherits(inputs[[i]], "tf"))
      if (!is.ts(inputs[[i]]$x)) 
        inputs[[i]]$x <- ts(inputs[[i]]$x, start = start, frequency = s)
      if (frequency(inputs[[i]]$x) != s) stop("invalid frequency for input")
      # Check sample period for input: start1 <= start and end1 >= end
      start1 <- start(inputs[[i]]$x)
      t0 <- (start[1] - start1[1])*s + (start1[2] - start[2] + 1)
      end1 <- end(inputs[[i]]$x)
      t1 <- (end[1] - start1[1])*s + (end1[2] - start[2] + 1)
      if (t0 < 1 || (t1 - t0 + 1) != N) {
        stop("incompatible samples")
      }
      t1 <- t0 + N - 1
      inputs[[i]]$t.start <- t0
      inputs[[i]]$t.end <- t1
      if (inputs[[i]]$param[1] == 0) {
        x <- diffC(inputs[[i]]$x[t0:t1], noise$nabla, FALSE)
        inputs[[i]]$param[1] <- lm(y ~ x)$coefficients[2]
      }
    }
    
    param <- lapply(inputs, function(x) x$param)
    param <- unlist(param)
    param <- param[!duplicated(names(param))]
    param <- as.list(param)
  }
  
  if (!is.null(xreg)) {
    xn <- deparse(substitute(xreg))
    xreg <- as.matrix(xreg)
    kx <- ncol(xreg)
    cn <- colnames(xreg)
    if (is.null(cn)) {
      if (kx == 1) cn <- xn
      else cn <- paste(xn, 1:kx, sep = "")
    }
    if (!is.ts(xreg)) xreg <- ts(xreg, start = start, frequency = s)
    if (frequency(xreg) != s) stop("invalid frequency for xreg inputs")
    # Check sample period for X: start1 <= start and end1 >= end
    xreg <- window(xreg, start = start)
    start1 <- start(xreg)
    if (start1[1] != start[1] || start1[2] != start[2])
      stop("incompatible samples")
    if (nrow(xreg) < N)
      stop("insufficient obs. for xreg")
    colnames(xreg) <- cn
    X <- apply(xreg, 2, function(x) {
      x <- diffC(x, noise$nabla, FALSE)
      if (length(x) > n) x[1:n]
      else x
    })
    reg <- lm(y ~ X - 1)
    b <- reg$coefficients
    names(b) <- cn
    param <- c(b, param)
  } else kx <- 0
  
  if (is.um(noise)) {
    param <- param[!names(param) %in% names(noise$param)]
    param <- c(param, noise$param)
  }
  
  mod <- list(output = output, xreg = xreg, inputs = inputs, noise = noise,
              param = param, kx = kx, k = k, optim = NULL, method = NULL,
              call = call)
  class(mod) <- "tfm"
  if (fit) mod <- fit.tfm(mod, envir = envir, ...)

  return(mod)

}

#' Coefficients of a transfer function model
#'
#' \code{coef} extracts the "coefficients" from a TF model.
#'
#' @param object a \code{tfm} object.
#' @param ... other arguments.
#'
#' @return A numeric vector.
#'
#' @export
coef.tfm <- function(object, ...) {
  param.tfm(object)
}


#' @rdname diagchk
#' @param y an object of class \code{ts}.
#' @param envir environment in which the function arguments are evaluated.
#'    If NULL the calling environment of this function will be used.
#' @export
diagchk.tfm <- function(mdl, y = NULL, method = c("exact", "cond"),
                        lag.max = NULL, lags.at = NULL, freq.at = NULL,
                        std = TRUE, envir = NULL, ...) {
  if (is.null (envir)) envir <- parent.frame ()
  u <- residuals.tfm(mdl, y, method, envir = envir)
  ide(u, graphs = c("plot", "hist", "acf", "pacf", "cpgram"), ylab = "u",
      lag.max = lag.max, lags.at = lags.at, freq.at = freq.at,
      std = std, envir=envir, ...)
}


.diff_tfm <- function(mdl, nabla) {
  if (mdl$kx > 0) {
    xreg <- apply(mdl$xreg, 2, function(x) {
      x <- diffC(x, nabla, FALSE)
    })
    mdl$xreg <- xreg
  }
  
  if (mdl$k > 0) {
    d <- length(nabla) - 1
    for (j in 1:mdl$k) {
      end <- end(mdl$inputs[[j]]$x)
      s <- frequency(mdl$inputs[[j]]$x)
      x <- diffC(mdl$inputs[[j]]$x, nabla, mdl$inputs[[j]]$um$bc)
      mdl$inputs[[j]]$x <- ts(x, end = end, frequency = s)
      n <- mdl$inputs[[j]]$t.end - mdl$inputs[[j]]$t.start + 1
      mdl$inputs[[j]]$t.end <- mdl$inputs[[j]]$t.end - d
      n <- n - d
      mdl$inputs[[j]]$t.start <- mdl$inputs[[j]]$t.end - n + 1
      mdl$inputs[[j]]$um$bc <- FALSE
    }
  }
  
  nabla <- polydivC(mdl$noise$nabla, nabla, FALSE)
  mdl$noise$nabla <- nabla
  mdl$noise$i <- list(as.lagpol(nabla))
    
  mdl
}


#' @rdname fit
#' @param y a \code{ts} object.
#' @param fit.noise logical. If TRUE parameters of the noise model are fixed.
#' @param envir environment in which the function arguments are evaluated.
#'    If NULL the calling environment of this function will be used.
#' @param ... additional arguments.
#'
#' @return A \code{tfm} object.
#' @export
fit.tfm <- function(mdl, y = NULL, method = c("exact", "cond"),
                    optim.method = "BFGS", show.iter = FALSE,
                    fit.noise = TRUE, envir = NULL, ...){

  stopifnot(inherits(mdl, "tfm"))
  if (is.null (envir)) envir <- parent.frame ()

  if (!fit.noise) {
    par.noise <- names(mdl$noise$param)
    mdl$noise$param <- unname(mdl$noise$param)
  }
  
  # if (!all(sapply(mdl$inputs, is.stationary.tf)))
  #   stop("Non-stationary AR preestimates for inputs")
  if (!is.stationary.um(mdl$noise))
    stop("Non-stationary AR preestimates")
  if (!is.invertible.um(mdl$noise))
    stop("Non-invertible MA preestimates")
  if (length(method) == 2 && !is.null(mdl$noise$method)) 
    method <- mdl$noise$method
  method <- match.arg(method)
  mdl$noise$method <- method
  exact <- method == "exact"

  noiseTFM <- function() {
    if (is.null(mdl$noise$mu)) ystar <- w
    else ystar <- w - mdl$noise$mu
    if (mdl$kx > 0) {
      ystar <- ystar - xreg %*% unlist(mdl$param[1:mdl$kx])
    }
    if (mdl$k > 0) {
      for (i in 1:mdl$k) {
        x <- filterC(X[[i]], mdl$inputs[[i]]$theta,
                     mdl$inputs[[i]]$phi, mdl$inputs[[i]]$delay)
        if (t0[i] > 1 || t1[i] > 1) ystar <- ystar - x[t0[i]:t1[i]]
        else ystar <- ystar - x
      }
    }
    ystar
  }
  
  logLikTFM <- function(b) {
    mdl <<- update.tfm(mdl, b)
    if (mdl$noise$is.adm) {
      wstar <- noiseTFM()
      if (exact) ll <- -ellarmaC(wstar, mdl$noise$phi, mdl$noise$theta)
      else ll <- -cllarmaC(wstar, mdl$noise$phi, mdl$noise$theta)
    } else {
      ll <- ll0
    }
    if (show.iter) print(c(loglik = ll, b))
    mdl$noise$is.adm <<- TRUE
    return(ll)
  }

  y <- output.tfm(mdl, y, envir)
  start <- start(y)
  end <- end(y)
  s <- frequency(y)
  
  w <- diffC(y, mdl$noise$nabla, mdl$noise$bc)
  N <- length(y)
  n <- length(w)
  d <-  N - n
  
  if (mdl$kx > 0) {
    xreg <- apply(mdl$xreg, 2, function(x) {
      x <- diffC(x, mdl$noise$nabla, FALSE)
      if (length(x) > n) x[1:n]
      else x
    })
  }
  
  if (mdl$k > 0) {
    t0 <- rep(1, mdl$k)
    t1 <- rep(1, mdl$k)
    i <- 1
    X <- lapply(mdl$inputs, function(tf) {
      x <- diffC(tf$x, mdl$noise$nabla, tf$um$bc)
      t0[i] <<- tf$t.start
      t1[i] <<- tf$t.start + n - 1
      x
    })
  }

  b <- param.tfm(mdl)
  ll0 <- 1.797693e+308
  ll0 <- logLikTFM(b)
  opt <- optim(b, logLikTFM, method = optim.method, hessian = F)
  if(opt$convergence > 0)
    warning(gettextf("possible convergence problem: optim gave code = %d",
                     opt$convergence), domain = NA)

  mdl <- update.tfm(mdl, opt$par)
  wstar <- noise.tfm(mdl, y, diff = TRUE, envir=envir)
  if (!is.null(mdl$noise$mu)) wstar <- wstar - mdl$noise$mu
  if (method == "cond")
    res <- condresC(wstar, mdl$noise$phi, mdl$noise$theta, TRUE)
  else res <- gresC(wstar, mdl$noise$phi, mdl$noise$theta)
  mdl$noise$sig2 <- sum(res^2)/length(res)
  mdl$noise$optim <- opt
  mdl$noise$method <- method
  
  if (!fit.noise) names(mdl$noise$param) <- par.noise 
  
  return(mdl)
  
}


#' @export
logLik.tfm <-function(object, y = NULL, method = c("exact", "cond"), envir=NULL, ...) {
  method <- match.arg(method)
  if (is.null (envir)) envir <- parent.frame ()
  w <- noise.tfm(object, y, TRUE, envir=envir)
  if (!is.null(object$noise$mu))
    w <- w - object$noise$mu
  if (method == "exact") ll <- ellarmaC(w, object$noise$phi, object$noise$theta)
  else ll <- cllarmaC(w, object$noise$phi, object$noise$theta)
  return(ll)
  
}


#' @rdname modify
#' @export
modify.tfm <- function(mdl, ...) {
  args <- list(...)
  ar <- args$ar
  i <- args$i
  ma <- args$ma
  bc <- args$bc
  sig2 <- args$sig2
  fit <- args$fit
  um <- modify.um(mdl$noise, ar = ar, i = i, ma = ma, bc = bc, sig2 = sig2,
                  fit = fit)
  tfm(xreg = mdl$xreg, inputs = mdl$inputs, noise = um)
}


#' Noise of a transfer function model
#' 
#' \code{noise} computes the noise of a linear transfer function model.     
#'
#' @param tfm an object of the class \code{tfm}.
#' @param y output of the TF model if it is different to that of the 
#' \code{tfm} object.
#' @param diff logical. If TRUE, the noise is differenced with
#' the "i" operator of the univariate model of the noise.
#' @param exp logical. If TRUE, the antilog transformation is applied.
#' @param envir environment in which the function arguments are evaluated.
#'    If NULL the calling environment of this function will be used.
#' @param ... additional arguments.
#'
#' @return A "ts" object.
#' 
#' @export
noise <- function (tfm, ...) { UseMethod("noise") }

#' @rdname noise
#' @export
noise.tfm <- function(tfm, y = NULL, diff = TRUE, exp = FALSE, envir = NULL, ...) {
  y <- output.tfm(tfm, y, envir = envir)
  start <- start(y)
  end <- end(y)
  s <- frequency(y)
  N <- length(y)
  if (tfm$noise$bc) y <- log(y)
  
  if (diff & tfm$noise$d > 0) {
    y <- diffC(y, tfm$noise$nabla, FALSE)
    n <- length(y)
    if (tfm$kx > 0) {
      for (i in 1:tfm$kx) {
        x <- diffC(tfm$xreg[, i], tfm$noise$nabla, FALSE)
        if (length(x) > n) 
          y <- y - x[1:n, ]*tfm$param[[i]]
        else
          y <- y - x*tfm$param[[i]]
      }
    }
    if (tfm$k > 0) {
      for (i in 1:tfm$k) {
        x <- diffC(tfm$inputs[[i]]$x, tfm$noise$nabla, tfm$inputs[[i]]$um$bc)
        x <- filterC(x, tfm$inputs[[i]]$theta,
                     tfm$inputs[[i]]$phi, tfm$inputs[[i]]$delay)
        t0 <- tfm$inputs[[i]]$t.start
        t1 <- tfm$inputs[[i]]$t.end
        if (t0 > 1 || length(tfm$inputs[[i]]$x) > t1) 
          y <- y - x[t0:(t0+n-1)]
        else
          y <- y - x
      }
    }
  }
  else {
    if (tfm$kx > 0) {
      if (nrow(tfm$xreg) > N) {
        y <- y - as.matrix(tfm$xreg[1:N, ]) %*% unlist(tfm$param[1:tfm$kx])
      } else
        y <- y - tfm$xreg %*% unlist(tfm$param[1:tfm$kx])
    }
    if (tfm$k > 0) {
      for (i in 1:tfm$k) {
        x <- filterC(tfm$inputs[[i]]$x, tfm$inputs[[i]]$theta,
                     tfm$inputs[[i]]$phi, tfm$inputs[[i]]$delay)
        t0 <- tfm$inputs[[i]]$t.start
        t1 <- tfm$inputs[[i]]$t.end
        if (t0 > 1 | length(tfm$inputs[[i]]$x) > t1)
          y <- y - x[t0:t1]
        else
          y <- y - x
      }
    }
    if (tfm$noise$bc && exp) y <- exp(y)
  }
  
  y <- ts(y, end = end, frequency = s)
  if (!is.null(ncol(y))) y <- y[, 1]
  y
}

#' @rdname setinputs
#' @export
setinputs.tfm <- function(mdl, xreg = NULL, inputs = NULL, y = NULL, 
                          envir = parent.frame (), ...) {
  stopifnot(is.tfm(mdl))
  if (!is.null(y)) mdl$noise$z <- deparse(substitute(y))
  y <- output.tfm(mdl, y, envir)
  if (!is.null(xreg)) {
    if (mdl$kx > 0) {
      name <- c(colnames(mdl$xreg), colnames(xreg))
      xreg <- cbind(mdl$xreg, xreg)
      colnames(xreg) <- name
    }
  } else xreg <- mdl$xreg
  
  if (!is.null(inputs) && mdl$k > 0) inputs <- c(mdl$inputs, inputs)  
  else if(mdl$k > 0) inputs <- mdl$inputs
  tfm(output = y, xreg = xreg, inputs = inputs, noise = mdl$noise, 
      new.name = FALSE, envir = envir, ...)
}


#' Summarizing Transfer Function models
#' 
#' \code{summary} method for class "tfm".
#' 
#' @param object a \code{tfm} object.
#' @param y a "ts" object.
#' @param method exact or conditional maximum likelihood.
#' @param digits number of significant digits to use when printing.
#' @param envir environment in which the function arguments are evaluated.
#'    If NULL the calling environment of this function will be used.
#' @param ... additional arguments.
#' @return A \code{tfm} object.
#' @export
summary.tfm <- function(object, y = NULL, method = c("exact", "cond"),
                        digits = max(3L, getOption("digits") - 3L), 
                        envir=NULL, ...) {

  stopifnot(inherits(object, "tfm"))
  if (is.null (envir)) envir <- parent.frame ()
  model.name <- deparse(substitute(object))

  args <- list(...)
  if(is.null(args[["p.values"]])) p.values <- FALSE
  else p.values <- args[["p.values"]]
  if(is.null(args[["table"]])) table <- FALSE
  else table <- args[["table"]]

  method <- match.arg(method)
  if (!is.null(object$noise$method)) method <- object$noise$method

  logLikTFM <- function(b){
    object <<- update.tfm(object, b)
    wstar <- noise.tfm(object, y, diff = TRUE, envir=envir)
    if (!is.null(object$noise$mu)) wstar <- wstar - object$noise$mu
    if (method == "cond") ll <- cllarmaC(wstar, object$noise$phi, object$noise$theta)
    else ll <- ellarmaC(wstar, object$noise$phi, object$noise$theta)
    return(ll)
  }

  resTFM <- function(b) {
    object <<- update.tfm(object, b)
    wstar <- noise.tfm(object, y, diff = TRUE, envir=envir)
    if (!is.null(object$noise$mu)) wstar <- wstar - object$noise$mu
    if (method == "cond")
      res <- condresC(wstar, object$noise$phi, object$noise$theta, TRUE)
    else res <- gresC(wstar, object$noise$phi, object$noise$theta)
    as.vector(res)
  }
  
  y <- output.tfm(object, y, envir)
  w <- diffC(y, object$noise$nabla, object$noise$bc)
  N <- length(y)
  n <- length(w)
  b <- param.tfm(object)
  ll <- logLikTFM(b)
  aic <- (-2.0*ll+2*length(b))/n
  bic <- (-2*ll+log(n)*length(b))/n
  res <- resTFM(b)
  J <- numDeriv::jacobian(resTFM, b)
  g <- t(J) %*% res
  ssr <- sum(res^2)
  object$noise$sig2 <- ssr/length(res)
  H <- t(J) %*% J
  varb <- solve(H)*object$noise$sig2
  b <- param.tfm(object)
  se <- sqrt(diag(varb))
  z <- b/se
  p <- pnorm(abs(z), lower.tail = F)*2

  if (p.values) return(p)
  X <- cbind(b, g, se, z, p)
  colnames(X) <- c("Estimate", "Gradient", "Std. Error", "z Value", "Pr(>|z|)")
  rownames(X) <- names(b)
  
  if (table) return(X)

  res <- noise.tfm(object, y, envir = envir)
  if (!is.null(object$noise$mu)) res <- res - object$noise$mu
  if (method == "cond") {
    res <- as.numeric(condresC(res, object$noise$phi, object$noise$theta, TRUE))
  } else{
    res <- as.numeric(exactresC(res, object$noise$phi, object$noise$theta))
  }
  
  tss <- var(y)*(N-1)
  mean.resid <- mean(res)
  rss <- var(res)*(length(res)-1)
  sd.resid <- sd(res)
  z.mean.resid <- mean.resid*sqrt(length(res))/sd.resid
  p.mean.resid <- pnorm(abs(z.mean.resid), lower.tail = F)*2
  if (is.ts(y)) {
    start <- start(y)
    end <- end(y)
    s <- frequency(y)
  } else {
    start <- NULL
    end <- NULL
    s <- NULL
  }
  
  Q1 <- Box.test(res, lag = object$noise$p+object$noise$q+1, 
                 type = "Ljung-Box", fitdf = object$noise$p+object$noise$q)
  Q2 <- Box.test(res, lag = as.integer(n/4)+object$noise$p+object$noise$q,
                 type = "Ljung-Box", fitdf = object$noise$p+object$noise$q)
  Q <- c(Q1 = Q1, Q2 = Q2)
  
  g <- rep(3, length(res))
  g[1:floor(length(res)/3)] <- 1
  g[(floor(length(res)/3)+1):floor(2*length(res)/3)] <- 2
  h.stat <- bartlett.test(res, g = g)  
  
  if (is.ts(y))
    res <- ts(res, end = end(y), frequency = frequency(y))
  out <- list(call = object$call, model.name = model.name, 
              ml.method = object$noise$method,
              z = object$noise$z, aic = aic, bic = bic, gradient = g, 
              var.coef = varb, logLik = ll, N = N, n = n, 
              mean.resid = mean.resid, sd.resid = sd.resid, 
              z.mean.resid = z.mean.resid, p.mean.resid = p.mean.resid,
              resid = res, rss = ssr, sig2 = object$noise$sig2, Q = Q, 
              h.stat = h.stat, tss = tss, table = X, start = start,
              end = end, s = s)
  class(out) <- "summary.tfm"
  
  out
  
}

#' @export
print.summary.tfm <- function(x, stats = TRUE, short = FALSE,
                              digits = max(3L, getOption("digits") - 3L), ...) {
  print.summary.um(x, stats, digits, short = short, ...)
}

#' @rdname calendar
#' @export
calendar.tfm <-
  function(mdl, y = NULL, form = c("dif", "td", "td7", "td6", "wd"),
           ref = 0, lom = TRUE, lpyear = TRUE, easter = FALSE, len = 4, 
           easter.mon = FALSE, n.ahead = 0, p.value = 1, envir = NULL, ...)
{
  if (is.null (envir)) envir <- parent.frame ()
  if (is.null(y)) y <- output.tfm(mdl, y, envir = envir)
  else mdl$noise$z <- deparse(substitute(y))
  if (frequency(y) != 12) stop("function only implemented for monthly ts")

  if (is.null(n.ahead)) n.ahead <- 0  
  n.ahead <- abs(n.ahead)
  xreg <- CalendarVar(y, form, ref, lom, lpyear, easter, len, easter.mon, n.ahead)
  tfm1 <- setinputs.tfm(mdl, xreg, y = y)
  if (p.value < 0.999) {
    p <- summary.tfm(tfm1, p.values = TRUE)
    p <- (p[1:ncol(xreg)] <= p.value)
    if (all(p)) return(tfm1)
    if (any(p)) {
      xreg <- xreg[, p]
      tfm1 <- tfm(y, xreg = xreg, noise = tfm1$noise, envir = envir, new.name = FALSE, ...)
      return(tfm1)
    }
    return(tfm1$noise)
  }
  return(tfm1)
}


#' @rdname outliers
#' @param y an object of class \code{ts}, optional.
#' @param envir environment in which the function arguments are evaluated.
#'    If NULL the calling environment of this function will be used.
#' @export
outliers.tfm <- function(mdl, y = NULL, types = c("AO", "LS", "TC", "IO"), 
                         dates = NULL,  c = 3, calendar = FALSE, easter = FALSE, 
                         resid = c("exact", "cond"), n.ahead = NULL, 
                         p.value = 1, tc.fix = TRUE, envir=NULL, ...) {
  if (is.null (envir)) envir <- parent.frame ()
  if (is.null(n.ahead)) n.ahead <- 0L
  else n.ahead <- abs(n.ahead)
  if (!is.null(y)) mdl$noise$z <- deparse(substitute(y))
  y <- output.tfm(mdl, y, envir)
  if (mdl$kx > 0)
    n.ahead <- length(y) - nrow(mdl$xreg)
  N <- noise.tfm(mdl, y, diff = FALSE, envir=envir)
  start <- start(N)
  freq <- frequency(N)
  if( freq < 2) {
    calendar <- FALSE
    easter <- FALSE
  }
  if ((mdl$noise$p > 50 || mdl$noise$q > 50) && length(resid) > 1)
    resid <- "cond"
  resid <- match.arg(resid)
  eres <- resid == "exact"
  if (is.numeric(types)) {
    if (length(types) == 1) {
      if (types == 1)types <- "AO"
      else if (types == 2) types <- c("AO", "LS")
      else if (types == 3) types <- c("AO", "LS", "TC")
      else types <- c("AO", "LS", "TC", "IO")
    } else {
      types <- c("AO", "LS", "TC", "IO")[types]
    }    
  }
  types <- toupper(types)
  types <- match.arg(c("AO", "LS", "TC", "IO"), types, several.ok = TRUE)
  types <- sapply(c("AO", "LS", "TC", "IO"), function(x) (x %in% types)*1L)
  
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
  
  bc <- mdl$noise$bc
  mdl$noise$bc <- FALSE
  if (is.null(mdl$noise$mu)) mu <- 0
  else mu <- mdl$noise$mu
  
  tfm1 <- NULL
  if (calendar||easter) {
    if (calendar)
      tfm1 <- calendar.um(mdl$noise, N, easter = easter, 
                          n.ahead = n.ahead, envir=envir, ...)
    else
      tfm1 <- easter.um(mdl$noise, N, n.ahead, envir = envir, ...)
    N <- noise.tfm(tfm1, diff = FALSE, envir=envir)
    A <- outliersC(N, FALSE,  mu, tfm1$noise$phi, tfm1$noise$nabla, 
                   tfm1$noise$theta, types, indx, eres, abs(c))
  } else {
    A <- outliersC(N, FALSE,  mu, mdl$noise$phi, mdl$noise$nabla, 
                   mdl$noise$theta, types, indx, eres, abs(c))
  }
  
  if (ncol(A) == 1) {
    warning("no outlier")
    if (is.null(tfm1)) return(mdl)
    else return(tfm1)
  }
  
  df <- data.frame(obs = A[, 1], date = time(N)[A[, 1]], type = A[, 2],
                   effect = A[, 3], t.ratio = A[, 4], w1 = A[, 5], t.w1 = A[, 6])
  df[[3]] <- factor(df[[3]], levels = 1:4, labels = c("IO", "AO", "LS", "TC"))
  df$date <- sapply(df$date, function(x) {
    year <- trunc(x)
    if (freq > 1) {
      m <- round( (x - year)*freq + 1)
      if (freq == 12 && m < 10) m <- paste(0, m, sep = "")
      year <- as.character(year)
      if (nchar(year) == 4) year <- substr(year, 3, 4)
      paste(year, m, sep = ".")
    } else as.character(year)
  })
  
  n <- length(N) + n.ahead
  if (tc.fix)
    i <- df[, 3] == "AO" | df[, 3] == "LS" | df[, 3] == "TC" 
  else
    i <- df[, 3] == "AO" | df[, 3] == "LS"
  if (any(i)) {
    names <- c()
    X <- sapply((1:nrow(df))[i], function(k) {
      p <- double(n)
      p[df[k, 1]] <- 1
      names <<- c(names, paste0(df[k, 3], df[k, 2]))
      if (df[k, 3] == "AO") p
      else if (df[k, 3] == "LS") cumsum(p)
      else {
        p <- cumsum(p)
        p*(0.7^(cumsum(p)-1))
      }
    })
    if (is.vector(X)) X <- matrix(X, ncol = 1)
    X <- ts(X, start = start, frequency = freq)
    colnames(X) <- names
  } else X <- NULL

  if (any(!i)) {
    tfi <- lapply((1:nrow(df))[!i], function(k) {
      p <- double(n)
      p[df[k, 1]] <- 1
      p <- ts(p, start = start, frequency = freq)
      prefix <- paste(df[k, 3], df[k, 2], sep = "")
      if (df[k, 3] == "IO") tf(p, w0 = df[k, 4], ma = mdl$noise$ma, ar = c(mdl$noise$i, mdl$noise$ar), par.prefix = prefix)
      else if (df[k, 3] == "TC") tf(p, w0 = df[k, 4], ar = 1, par.prefix = prefix)
      else stop("unkown input")
    })
  } else tfi <- NULL 
  
  if (mdl$kx > 0) {
    names <- colnames(mdl$xreg)
    xreg <- mdl$xreg
    if (calendar||easter) {
      names <- c(names, colnames(tfm1$xreg))
      xreg <- cbind(xreg, tfm1$xreg)
    }
    if (any(i)) {
      names <- c(names, colnames(X))
      xreg <- cbind(xreg, X)
    }
    colnames(xreg) <- names
  } else if (calendar||easter) {
    names <- colnames(tfm1$xreg)
    xreg <- tfm1$xreg
    if (any(i)) {
      names <- c(names, colnames(X))
      xreg <- cbind(xreg, X)
    }
    colnames(xreg) <- names
  } else if (any(i))
    xreg <- X
  else
    xreg <- NULL
  
  if (is.null(tfi)) tfi <- mdl$inputs
  else if (!is.null(mdl$inputs)) tfi <- list.tf(mdl$inputs, tfi)
  mdl$noise$bc <- bc
  tfm1 <- tfm(y, xreg = xreg, inputs = tfi, noise = mdl$noise, fit = TRUE, 
              new.name = FALSE, envir = envir)
  if (p.value < 0.999) 
    tfm1 <- varsel.tfm(tfm1, NULL, p.value = p.value, envir = envir)
  
  return(tfm1)
  
}


#' @rdname intervention
#' @export
intervention.tfm <- function(mdl, y = NULL, type, time, n.ahead = 0, 
                             envir = parent.frame(), ...) {
  stopifnot(is.tfm(mdl))
  type <- toupper(type)
  y <- output.tfm(mdl, y, envir = envir)
  if (any(type == "IO")) 
    u <- residuals(mdl, y, envir = envir)
  n <- noise.tfm(mdl, y, diff = FALSE, exp = TRUE, envir = envir)
  s <- frequency(y)
  k1 <- length(type)
  # convert time to list of dates
  if (is.list(time))
    k2 <- length(time)
  else if (is.matrix(time)) {
    if (ncol(time) == 1 || nrow(time) == 1) {
      if (s > 1) stop("invalid time argument")
      time <- lapply(as.vector(time), function(x) c(x, 1))
    } else if (ncol(time) == 2)
      time <- lapply(1:nrow(time), function(x) time[x, ])
    else if (nrow(time) == 2) # ncol(time) > 2
      time <- lapply(1:ncol(time), function(x) time[, x])
    else
      stop("invalid time argument")
    k2 <- length(time)
  } else if (is.vector(time)) {
    k2 <- length(time)
    if (k2 == 1) {
      if (s != 1) stop("invalid time argument")
      time <- list(c(time, 1))
    } else if (k2 == 2) {
      stopifnot(time[2] >= 1 && time[2] <=s)
      time <- list(time)
      k2 <- 1
    } else {
      stopifnot(s == 1)
      time <- lapply(time, function(x) c(x, 1))
    }
  } else stop("invalid time argument")

  testing <- FALSE
  if (k1 == 1 && k2 > k1) { 
    type <- rep(type, k2)
    k1 <- k2
  } else if (k2 == 1 && k1 > 1) {
    time <- time[[1]]
    time <- lapply(1:k1, function(x) time)
    k2 <- k1
    testing <- TRUE
  } else if (k1 != k2) stop("type and time are incompatible arguments")
  
  X <- sapply(1:k1, function(k) {
    if (type[k] == "AO" || type[k] == "P" || type[k] == "IO" || type[k] == "TC")
      InterventionVar(y, time[[k]], "P", n.ahead)
    else if (type[k] == "LS" || type[k] == "S")
      InterventionVar(y, time[[k]], "S", n.ahead)
    else if (type[k] == "R")
      InterventionVar(y, time[[k]], "R", n.ahead)
    else
      stop("type of intervention/outlier unknown")
  })
  
  names <- sapply(1:k1, function(k) {
    if (s > 1) paste0(type[k], time[[k]][1], ".", time[[k]][2])
    else paste0(type[k], time[[k]][1])
  })
  colnames(X) <- names

  i <- type != "IO" & type != "TC" # xreg inputs
  if (testing) {
    tbl <- sapply(1:k1, function(k) {
      x <- X[, k]
      if (i[k]) {
        x <- matrix(x, ncol = 1)
        colnames(x) <- names[k]
        tfm1 <- setinputs(mdl, xreg = x, ...)
      } else {
        if (type[k] == "TC") {
          tf <- tfest(n, x, p = 1, q = 0, um.y = mdl$noise, 
                      par.prefix = names[k], envir = envir)
          tf$x.name <- names[k]
        } else { # "IO"
          tf <- tf(x, w0 = tsvalue(u, time[[k]]), ma = mdl$noise$ma,
                   ar = c(mdl$noise$i, mdl$noise$ar), par.prefix = names[k], 
                   envir = envir)
        }
        tfm1 <- setinputs(mdl, inputs = tf, ...)
      }
      if (type[k] == "TC") {
        s <- summary(tfm1, table = TRUE)
        d <- s[paste0(names[k], ".d1"), 1]
        name <- names[k]
        names[k] <<- paste0(names[k], "(", format(d, digits = 2), ")")
        s[name, ]
      } else summary(tfm1, table = TRUE)[names[k], ]
    })
    colnames(tbl) <- names
    return(tbl)
  } else {
      if (all(i)) {
        return(setinputs(mdl, xreg = X, ...))
      } else {
        if (any(i)) {
          X1 <- X[, i]
          if (is.vector(X1))
            X1 <- matrix(X1, ncol = 1)
          colnames(X1) <- names[i]
        } else {
          X1 <- NULL
        }
        ltf <- lapply( (1:k1)[!i], function(k) {
          x <- X[, k]
          x <- ts(x, start = start(y), frequency = frequency(y))
          if (type[k] == "TC") {
            tf <- tfest(n, x, p = 1, q = 0, um.y = mdl$noise, 
                        par.prefix = names[k]) 
            tf$x.name <- names[k]
          } else { # "IO"
            tf <- tf(x, w0 = tsvalue(u, time[[k]]), ma = mdl$noise$ma,
                     ar = c(mdl$noise$i, mdl$noise$ar), par.prefix = names[k])
          }
          return(tf)
        })
        return(setinputs(mdl, xreg = X1, inputs = ltf, y = y, envir = envir, ...))
      }
  }
}

#' Forecasting with transfer function models
#'
#' \code{predict} computes point and interval predictions for a time series
#' based on a \code{tfm} object.
#' 
#' @param object an object of class \code{\link{um}}.
#' @param y an object of class \code{\link{ts}}.
#' @param ori the origin of prediction. By default, it is the last observation.
#' @param n.ahead number of steps ahead.
#' @param level confidence level.
#' @param i transformation of the series \code{y} to be forecasted. It is a
#' lagpol as those of a \code{\link{um}} object.  
#' @param envir environment in which the function arguments are evaluated.
#'    If NULL the calling environment of this function will be used.
#' @param newdata new data for the predictors for the forecast period. This is
#'   a matrix if there is more than one predictor. The number of columns is
#'   equal to the number of predictors, the number of rows equal to
#'   \code{n.ahead}. If there is one predictor only the data may be provided
#'   alternatively as a vector.    
#' @param ... additional arguments.
#' 
#' @details Forecasts for the inputs of a \code{tfm} object can be provided
#' in tree ways: (1) extending the time series with forecasts so that the length
#' of the intput is greater than the length of the output, (2) computed 
#' internally from the \code{um} object associated to the input and (3) with 
#' the \code{newdata} argument.  
#'
#' @export predict.tfm
#' @export
predict.tfm <- function(object, newdata=NULL, y = NULL, ori = NULL, n.ahead = NULL,
                        level = 0.95, i = NULL,  envir=NULL, ...) {
  stopifnot(is.tfm(object))
  if (is.null (envir)) envir <- parent.frame ()
  y <- output.tfm(object, y, envir=envir)
  if (!is.null(i)) {
    i <- .lagpol0(i, "i","delta", envir=envir)
    nabla <- polyexpand(i)
    object <- .diff_tfm(object, nabla)
    end <- end(y)
    s <- frequency(y)
    y <- as.vector(diffC(y, nabla, object$noise$bc))
    if (object$noise$bc){
      y <- y*100
      object$noise$sig2 <- object$noise$sig2*100^2
    }
    object$noise$bc <- FALSE
    y <- ts(y, end = end, frequency = s)
  }
  
  z <- noise.tfm(object, y, FALSE, envir=envir)
  if (is.null(object$noise$mu)) mu <- 0
  else mu <- object$noise$mu

  if (is.null(ori)) ori <- length(z)
  if (is.null(n.ahead)) n.ahead <- frequency(z)
  
  X <-  forecastC(z, FALSE, mu, object$noise$phi, object$noise$nabla,
                  object$noise$theta, object$noise$sig2, ori, n.ahead)
  t <- (ori+1):(ori+n.ahead)
  z <- ts(X[, 1], start = start(z), frequency = frequency(z))
  start <- start(z)
  end <- end(z)
  n <- length(z)
  s <- frequency(z)
  
  if (!is.null(newdata)) {
    newdata <- as.matrix(newdata)
    stopifnot(nrow(newdata) >= n.ahead)
  }
  
  if (object$kx > 0) {
    if (nrow(object$xreg) < n) {
      if (is.null(newdata)) 
        stop("insufficient number of forecasts for input") 
      if (col(newdata) < object$kx) 
        stop("wrong object 'newdata'")
      Xf <- newdata[1:n.ahead, 1:object$kx]
      if (ncol(newdata) > object$kx) 
        newdata <- newdata[,object$kx:ncol(newdata)]
      else 
        newdata <- NULL
    } else Xf <- as.matrix(object$xreg[t, ])
    z[t] <- z[t] +  Xf %*% unlist(object$param[1:object$kx])
  }
  
  if (object$k > 0) {
    for (i in 1:object$k) {
      start1 <- start(object$inputs[[i]]$x)
      if ( any(colnames(newdata) == object$inputs[[i]]$x.name) ) {
        x <- newdata[, object$inputs[[i]]$x.name]
        t1 <- object$inputs[[i]]$t.end
        object$inputs[[i]]$x <- c(object$inputs[[i]]$x[1:t1], x)
      } else if (length(object$inputs[[i]]$x) - object$inputs[[i]]$t.start + 1 < n) {      
        if (has.um.tf(object$inputs[[i]])) {
          end1 <- end(object$inputs[[i]]$x)
          nahead <- (end[1] - end1[1])*s + end[2] - end1[2]
          object$inputs[[i]] <- 
            predict.tf(object$inputs[[i]], n.ahead = nahead)
        } else stop("insufficient number of forecasts for input")
      }
      x <- object$inputs[[i]]$x
      if (object$inputs[[i]]$um$bc) x <- log(x)
      x <- filterC(x, object$inputs[[i]]$theta,
                   object$inputs[[i]]$phi, object$inputs[[i]]$delay)
      x <- ts(as.numeric(x), start = start1, frequency = s)
      x <- window(x, start = start, end = end) 
      z[t] <- z[t] + x[t] 
      if (has.um.tf(object$inputs[[i]])) {
        v <- var.predict.tf(object$inputs[[i]], n.ahead)
        X[t, 4] <- X[t, 4] + v
      }
    }
  }
  
  dates <- time(zoo::as.zoo(z))
  if (any(level <= 0 || level >= 1)) level[level <= 0 || level >= 1] <- 0.95
  level <- unique(level)
  cv <- qnorm((1-level)/2, lower.tail = F)
  se <- sqrt(X[t, 4])
  z[1:ori] <- y
  if (object$noise$bc) {
    z[t] <- exp(z[t])
    low <- sapply(cv, function(x) z[t]*exp(-x*se)) 
    upp <- sapply(cv, function(x) z[t]*exp(x*se))
  } else {
    low <- sapply(cv, function(x) z[t] - x*se) 
    upp <- sapply(cv, function(x) z[t] + x*se)
  }
  
  out <- list(z = z, rmse = se, low = low, upp = upp, level = level,
              dates = dates, ori = ori, n.ahead = n.ahead, ori.date = dates[ori])
  class(out) <- c("predict.tfm", "predict.um") 
  out  
}


#' @export
print.tfm <- function(x, ...) {
  print(summary(x), short = TRUE, ...)
}


#' Residuals of a transfer function model
#'
#' \code{residuals} computes the exact or conditional residuals of a TF model.
#'
#' @param object a \code{tfm} object.
#' @param y output of the TF model (if it is different to that of the "tfm" 
#' object).
#' @param method a character string specifying the method to compute the
#' residuals, exact or conditional.
#' @param envir environment in which the function arguments are evaluated.
#'    If NULL the calling environment of this function will be used.
#' @param ... additional arguments.
#'
#' @return A "ts" object.
#'
#' @export
residuals.tfm <- function(object, y = NULL, method = c("exact", "cond"), envir=NULL, ...) {

  stopifnot(inherits(object, "tfm"))
  if (is.null (envir)) envir <- parent.frame ()
  method <- match.arg(method)

  w <- noise.tfm(object, y, envir=envir)
  if (!is.null(object$noise$mu)) w <- w - object$noise$mu
  if (method == "cond")
    res <- condresC(w, object$noise$phi, object$noise$theta, TRUE)
  else 
    res <- exactresC(w, object$noise$phi, object$noise$theta)
  
  if (is.ts(w))
    res <- ts(res[, 1], end = end(w), frequency = frequency(w))
  else
    res <- res[, 1]
  
  return(res)
  
}

#' Signal component of a TF model
#'
#' \code{signal} extracts the signal of a TF model.
#'
#' @param mdl an object of the class \code{tfm}.
#' @param y output of the TF model if it is different to that of the 
#' \code{tfm} object.
#' @param diff logical. If TRUE, the noise is differenced with
#' the "i" operator of the univariate model of the noise.
#' @param envir environment in which the function arguments are evaluated.
#'    If NULL the calling environment of this function will be used.
#' @param ... additional arguments.
#'
#' @return A "ts" object.
#' @export
signal <- function (mdl, ...) { UseMethod("signal") }

#' @rdname signal
#' @export
signal.tfm <- function(mdl, y = NULL, diff = TRUE, envir=NULL, ...) {
  if (is.null (envir)) envir <- parent.frame ()
  y <- output.tfm(mdl, y, envir=envir)
  start <- start(y)
  end <- end(y)
  s <- frequency(y)
  N <- length(y)
  
  y <- y*0
  if (diff && mdl$noise$d > 0) {
    y <- diffC(y, mdl$noise$nabla, FALSE)
    n <- length(y)
    if (mdl$kx > 0) {
      for (i in 1:mdl$kx) {
        x <- diffC(mdl$xreg[, i], mdl$noise$nabla, FALSE)
        if (length(x) > n) 
          y <- y + x[1:n, ]*mdl$param[[i]]
        else
          y <- y + x*mdl$param[[i]]
      }
    }
    if (mdl$k > 0) {
      for (i in 1:mdl$k) {
        x <- filterC(mdl$inputs[[i]]$x, mdl$inputs[[i]]$theta,
                     mdl$inputs[[i]]$phi, mdl$inputs[[i]]$delay)
        t0 <- mdl$inputs[[i]]$t.start
        t1 <- mdl$inputs[[i]]$t.end
        if (t0 > 1 || length(mdl$inputs[[i]]$x) > t1) 
          y <- y + diffC(x[t0:t1], mdl$noise$nabla, mdl$inputs[[i]]$um$bc)
        else
          y <- y + diffC(x, mdl$noise$nabla, mdl$inputs[[i]]$um$bc)
      }
    }
  }
  else {
    if (mdl$kx > 0) {
      if (nrow(mdl$xreg) > N)
        y <- y + as.matrix(mdl$xreg[1:N, ]) %*% unlist(mdl$param[1:mdl$kx])
      else
        y <- y + mdl$xreg %*% unlist(mdl$param[1:mdl$kx])
    }
    if (mdl$k > 0) {
      for (i in 1:mdl$k) {
        x <- filterC(mdl$inputs[[i]]$x, mdl$inputs[[i]]$theta,
                     mdl$inputs[[i]]$phi, mdl$inputs[[i]]$delay)
        t0 <- mdl$inputs[[i]]$t.start
        t1 <- mdl$inputs[[i]]$t.end
        if (t0 > 1 || length(mdl$inputs[[i]]$x) > t1)
          y <- y + x[t0:t1]
        else
          y <- y + x
      }
    }
  }
  
  ts(y, end = end, frequency = s)
  
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
#' \code{tsdiag.tfm} is a wrap of the stats::tsdiag function.
#'
#' @param object a fitted \code{um} object.
#' @param gof.lag	the maximum number of lags for a Portmanteau goodness-of-fit test
#' @param ... additional arguments. 
#' 
#' @seealso stats::tsdiag.
#' 
#' @export
tsdiag.tfm <- function(object, gof.lag = 10, ...)
{
  ## plot standardized residuals, acf of residuals, Ljung-Box p-values
  oldpar <- par(mfrow = c(3, 1))
  on.exit(par(oldpar))
  rs <- residuals(object)
  stdres <- rs/sqrt(object$noise$sig2)
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


#' Variable selection
#'
#' \code{varsel} omits non-significant inputs from a transfer function model.
#'
#' @param tfm a \code{tfm} object.
#' @param y a "ts" object.
#' @param p.value probability value to decide whether or not to omit an input.
#' @param envir environment in which the function arguments are evaluated.
#'    If NULL the calling environment of this function will be used.
#' @param ... other arguments.
#' @return A \code{tfm} object or a "um" if no input is significant at that level.
#' @export
varsel <- function (tfm, ...) { UseMethod("varsel") }

#' @rdname varsel
#' @export
varsel.tfm <- function(tfm, y = NULL, p.value = 0.10, envir = NULL, ...) {
  if (is.null (envir)) envir <- parent.frame()
  if (is.null(y)) y <- output.tfm(tfm, y, envir)
  else tfm$noise$z <- deparse(substitute(y))
  
  p <- summary.tfm(tfm, y, p.values = TRUE)
  b <- param.tfm(tfm)
  nms <- names(b)
  names(p) <- nms
  if (!is.null(tfm$xreg)) {
    P <- (p[1:ncol(tfm$xreg)] <= p.value)
    if (all(P)) xreg <- tfm$xreg
    else if (any(P)) {
      xreg <- tfm$xreg[, P]
    } else xreg <- NULL
  } else xreg <- NULL
  if (tfm$k > 0) {
    tfi <- lapply(tfm$inputs, function(tf) {
      if (p[as.character(tf$w0.expr)] < p.value) {
        if (tf$par.prefix == "TC" && tf$p == 1 && length(tf$param) == 2 && tf$delay == 0) {
          if ( p[names(tfm$inputs[[3]]$param)[2]] < p.value ) tf
          else {
            tf1 <- tf(tf$x, w0 = tf$w0, par.prefix = tf$par.prefix)
            tf1$x.name <- tf$x.name
            tf1
          }
        } else tf
      }
      else NULL
    })
    tfi[sapply(tfi, is.null)] <- NULL
  } else tfi <- NULL
  
  if (is.null(xreg) && is.null(tfi)) return(tfm$noise)
  
  tfm(y, xreg = xreg, inputs = tfi, noise = tfm$noise, new.name = FALSE, ...)
}

is.tfm <- function(tfm) {
  inherits(tfm, "tfm")
}


param.tfm <- function(tfm) {
  unlist(tfm$param, use.names = TRUE)
}

output.tfm <- function(tfm, y = NULL, envir = NULL) {
  if (is.null(y)) {
    if (is.ts(tfm$output)) return(tfm$output)
    if (is.null(tfm$output)) stop("argment y required")
    else {
      if (is.null (envir)) envir <- parent.frame ()
      y <- eval(parse(text = tfm$output), envir)
    }
  }  
  if (is.ts(y)) return(y)
  if (!is.numeric(y)) stop("output must be a numeric vector")
  ts(y)
} 

update.tfm <- function(tfm, param) {
  tfm$param[] <- param
  tfm$inputs <- lapply(tfm$inputs, update.tf, param = param)
  tfm$noise <- update.um(tfm$noise, param)
  return(tfm)
}


#' Cross-correlation check
#'
#' \code{ccf} displays ccf between prewhitened inputs and residuals.
#'
#' @param tfm a \code{tfm} object.
#' @param lag.max number of lags.
#' @param method Exact/conditional residuals.
#' @param envir environment in which the function arguments are evaluated.
#'    If NULL the calling environment of this function will be used.
#' @param ... additional arguments.
#'
#'
ccf.tfm <- function(tfm, lag.max = NULL, method = c("exact", "cond"), envir=NULL, ...) {

  if (is.null (envir)) envir <- parent.frame ()
  if (tfm$k < 1) stop("no stochastic input")

  j <- c()
  for (i in 1:tfm$k) 
    if (tfm$inputs[[i]]$um$k > 0)
      j <- c(j, i)
  
  k <- length(j)
  if (k < 1) stop("no stochastic input")
  
  
  if (k > 1) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow = c(k, 1))
  }

  u <- residuals.tfm(tfm, method = method, envir = envir)
  for (i in j) {
    end <- end(tfm$inputs[[i]]$x)
    s <- frequency(tfm$inputs[[i]]$x)
    x <- tfm$inputs[[i]]$x[tfm$inputs[[i]]$n.back:length(tfm$inputs[[i]]$x)]
    x <- ts(x, end = end, frequency = s)
    x <- residuals.um(tfm$inputs[[i]]$um, x, method, envir=envir)
    x <- window(x, start = start(u), end = end(u))
    pccf(x, u, lag.max = lag.max)
  }
  
  invisible(NULL)
  
}


#' @rdname ucomp
#' @param y an object of class \code{\link{ts}}.
#' @param envir environment in which the function arguments are evaluated.
#'    If NULL the calling environment of this function will be used.
#' @export
ucomp.tfm <- function(mdl, y = NULL,
                      method = c("mixed", "forecast", "backcast"), envir=envir, ...) {

  if (is.null (envir)) envir <- parent.frame ()
  y <- noise.tfm(mdl, y, FALSE, envir=envir)
  mdl$noise$bc <- FALSE
  uc <- ucomp.um(mdl$noise, y, method, envir=envir)
  return(uc)

}


#' @rdname sim
#' @param y0 initial conditions for the nonstationary series.
#' 
#' @export
sim.tfm <- function(mdl, n = 100, y0 = NULL, seed = NULL, ...) {
  
  z <- signal.tfm(mdl, diff = TRUE)
  end <- end(z)
  s <- frequency(z)
  
  d <- mdl$noise$d
  if (!is.null(seed)) set.seed(seed)
  a <- rnorm(n - d, 0, sqrt(mdl$noise$sig2))
  if (is.null(mdl$noise$mu)) mu <- 0
  else mu <- mdl$noise$mu
  
  z <- z + simC(a, FALSE, mu, mdl$noise$phi, 1, mdl$noise$theta, 0)
  if (d > 0) {
    y <- double(n)
    b <- -mdl$noise$nabla[2:(d+1)]
    if (!is.null(y0)) {
      if (length(y0) >= d) y[1:d] <- y0[1:d]
    } 
    
    for (i in (d+1):n) 
      y[i] <- z[i-d] + sum(y[(i-1):(i-d)] * b) 
  } else y <- z
  
  if (mdl$noise$bc) ts(exp(y), end = end, frequency = s)
  else ts(y, end = end, frequency = s)
  
}




