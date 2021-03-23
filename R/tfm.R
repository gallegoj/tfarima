#' Transfer function models
#' 
#' \code{tfm} creates a multiple input transfer function model.     
#'
#' @param output a ts object or a numeric vector.
#' @param xreg a matrix of regressors.
#' @param inputs a list of tf objects.
#' @param noise a um object for the noise.
#' @param fit logical. If TRUE, model is fitted. 
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
                ...) {
  stopifnot(is.um(noise))

  if (!is.null(inputs)) {  
    if (inherits(inputs, "tf")) { 
      inputs <- list(inputs)
    } else {
      inputs <- inputs[!sapply(inputs, is.null)]
    }
    k <- length(inputs)
  } else k <- 0
  
  if (is.null(output)) {
    if (is.um(noise)) {
      if (is.null(noise$z)) stop("missing output")
      else {
        txt <- noise$z
        output <- eval(parse(text = txt), .GlobalEnv)
      }
    }
  } else if (is.character(output)) {
    txt <- output
    output <- eval(parse(text = txt), .GlobalEnv)
  } else if (is.numeric(output)) {
    txt <- deparse(substitute(output))
    noise$z <- txt
  }
  
  N <- length(output)
  stopifnot(N > noise$d)
  start <- start(output)
  end <- end(output)
  s <- frequency(output)
  
  y <- diffC(output, noise$nabla, noise$bc)
  n <- length(y)
  
  param <- list()
  if (k > 0) {
    for (i in 1:k) {
      stopifnot(inherits(inputs[[i]], "tf"))
      if (frequency(inputs[[i]]$x) != s) stop("invalid frequency for input")
      start1 <- start(inputs[[i]]$x)
      end1 <- end(inputs[[i]]$x)
      t0 <- (start[1] - start1[1])*s + (start1[2] - start[2]) + 1
      t1 <- (end1[1] - start[1])*s + end1[2] - start[2] + 1 
      if (t0 < 1 || t1 < t0 || t1 < N) {
        stop("incompatible samples")
      }
      inputs[[i]]$t.start <- t0
      inputs[[i]]$t.end <- t1
      if (inputs[[i]]$param[1] == 0) {
        x <- diffC(inputs[[i]]$x[t0:(t0+N-1)], noise$nabla, FALSE)
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
    if (is.ts(xreg)) {
      if (frequency(xreg) != s) stop("invalid frequency for xreg inputs")
      start1 <- start(xreg)
      end1 <- end(xreg)
      t0 <- (start[1] - start1[1])*s + (start1[2] - start[2]) + 1
      t1 <- (end1[1] - start[1])*s + end1[2] - start[2] + 1 
      if (t0 < 1 || t1 < t0 || t1 < N) {
        stop("incompatible samples")
      }
      if (t0 > 1) xreg <- xreg[t0:t1, ]
    } else {
      xreg <- ts(xreg, start = start, frequency = s)
    }
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
  } else if (k == 0) stop("missing inputs")
  else kx <- 0
  
  if (is.um(noise)) {
    param <- param[!names(param) %in% names(noise$param)]
    param <- c(param, noise$param)
  }
  
  mod <- list(output = txt, xreg = xreg, inputs = inputs, noise = noise, 
              param = param, kx = kx, k = k, optim = NULL, method = NULL)
  class(mod) <- "tfm"
  if (fit) mod <- fit.tfm(mod, output, ...)
  
  return(mod)
  
}


#' @rdname diagchk
#' @param y an object of class \code{ts}.
#' @export
diagchk.tfm <- function(mdl, y = NULL, method = c("exact", "cond"), 
                        lag.max = NULL, std = TRUE, ...) {
  u <- residuals.tfm(mdl, y, method)
  ide(u, graphs = c("plot", "hist", "acf", "pacf", "cpgram"), lag.max = lag.max, std = std, ...)
}


#' @rdname fit
#' @param y a \code{ts} object.
#' @param fit.noise logical. If TRUE parameters of the noise model are fixed.
#' @param ... additional arguments.
#' 
#' @return A \code{tfm} object.
#' @export
fit.tfm <- function(mdl, y = NULL, method = c("exact", "cond"), 
                    optim.method = "BFGS", show.iter = FALSE, 
                    fit.noise = TRUE, ...){
  
  stopifnot(inherits(mdl, "tfm"))
  
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
  if (is.null(y)) {
    if (is.null(mdl$output)) stop("output required")
    else {
      y <- eval(parse(text = mdl$output), .GlobalEnv)
    }
  } 

  if (!is.ts(y)) y <- as.ts(y)
  
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
        ystar <- ystar - x
      }
    }
    ystar
  }
  
  logLikTFM <- function(b) {
    mdl <<- update.tfm(mdl, b)
    if (mdl$noise$is.adm) {
      wstar <- noiseTFM()
      if (exact) ll <- ellarmaC(wstar, mdl$noise$phi, mdl$noise$theta)
      else ll <- cllarmaC(wstar, mdl$noise$phi, mdl$noise$theta)
    } else {
      ll <- ll0
    }
    
    if (show.iter) print(c(ll, b))
    mdl$noise$is.adm <<- TRUE
    return(-ll)
  }
  
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
    X <- lapply(mdl$inputs, function(tf) {
      x <- diffC(tf$x, mdl$noise$nabla, tf$um$bc)
      t0 <- tf$t.start
      t1 <- tf$t.end
      if (t0 > 1 || t1 > N) x[t0:(t0+n-1)]
      else x
    })
  }

  b <- param.tfm(mdl)
  ll0 <- -logLikTFM(b)
  opt <- optim(b, logLikTFM, method = optim.method, hessian = F)
  if(opt$convergence > 0)
    warning(gettextf("possible convergence problem: optim gave code = %d",
                     opt$convergence), domain = NA)
  
  mdl <- update.tfm(mdl, opt$par)
  wstar <- noise.tfm(mdl, diff = TRUE)
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
logLik.tfm <-function(object, y = NULL, method = c("exact", "cond"), ...) {
  method <- match.arg(method)
  w <- noise.tfm(object, y, TRUE)
  if (!is.null(object$noise$mu))
    w <- w - object$noise$mu
  if (method == "exact") ll <- ellarmaC(w, object$noise$phi, object$noise$theta)
  else ll <- cllarmaC(w, object$noise$phi, object$noise$theta)
  
  return(c(ll, param.tfm(object)))
  
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
#' @param ... additional arguments.
#' 
#' @return A "ts" object.
#' 
#' @export
noise <- function (tfm, ...) { UseMethod("noise") }

#' @rdname noise
#' @export
noise.tfm <- function(tfm, y = NULL, diff = TRUE, exp = FALSE, ...) {
  y <- output.tfm(tfm, y)
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
        x <- filterC(tfm$inputs[[i]]$x, tfm$inputs[[i]]$theta,
                     tfm$inputs[[i]]$phi, tfm$inputs[[i]]$delay)
        t0 <- tfm$inputs[[i]]$t.start
        t1 <- tfm$inputs[[i]]$t.end
        if (t0 > 1 || t1 > N) 
          y <- y - diffC(x[t0:(t0+N-1)], tfm$noise$nabla, tfm$inputs[[i]]$um$bc)
        else
          y <- y - diffC(x, tfm$noise$nabla, tfm$inputs[[i]]$um$bc)
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
        if (t0 > 1 || t1 > N)
          y <- y - x[t0:(t0+N-1)]
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


#' \code{setinputs} adds new inputs into a transfer function model.     
#'
#' @param tfm a \code{tfm} object.
#' @param xreg a matrix of inputs.
#' @param inputs a list of tf objects.
#' @param ... further arguments to be passed to particular methods
#' 
#' @return A \code{tfm} object.
#' 
#' @export
setinputs <- function (tfm, ...) { UseMethod("setinputs") }

#' @rdname setinputs
#' @export
setinputs.tfm <- function(tfm, xreg = NULL, inputs = NULL, ...) {
  stopifnot(is.tfm(tfm))
  if (!is.null(xreg)) {
    if (!is.null(tfm$xreg))
      xreg <- cbind(xreg, tfm$xreg)
  } else xreg <- tfm$xreg
  
  if (!is.null(inputs)) {
    if (!is.null(tfm$inputs) && length(tfm$inputs) > 1) 
      inputs <- c(inputs, tfm$inputs)  
  } else inputs <- tfm$inputs
  tfm(tfm$output, xreg, inputs, tfm$noise)
  
}


#' Summarizing Transfer Function models
#' 
#' \code{summary} method for class "tfm".
#' 
#' @param object a \code{tfm} object.
#' @param y a "ts" object.
#' @param method exact or conditional maximum likelihood.
#' @param digits number of significant digits to use when printing.
#' @param ... additional arguments.
#' @return A \code{tfm} object.
#' @export
summary.tfm <- function(object, y=NULL, method = c("exact", "cond"),
                        digits = max(3L, getOption("digits") - 3L), ...) {
  
  stopifnot(inherits(object, "tfm"))
  model.name <- deparse(substitute(object))
  
  args <- list(...)
  if(is.null(args[["p.values"]])) p.values <- FALSE
  else p.values <- args[["p.values"]]
  
  if (is.null(y)) {
    if (is.null(object$output)) stop("argument y required")
    else y <- eval(parse(text = object$output), .GlobalEnv)
  }
  
  method <- match.arg(method)
  if (!is.null(object$noise$method)) method <- object$noise$method
  
  logLikTFM <- function(b){
    object <<- update.tfm(object, b)
    wstar <- noise.tfm(object, diff = TRUE)
    if (!is.null(object$noise$mu)) wstar <- wstar - object$noise$mu 
    if (method == "cond") ll <- cllarmaC(wstar, object$noise$phi, object$noise$theta)
    else ll <- ellarmaC(wstar, object$noise$phi, object$noise$theta)
    return(ll)
  }
  
  resTFM <- function(b) {
    object <<- update.tfm(object, b)
    wstar <- noise.tfm(object, diff = TRUE)
    if (!is.null(object$noise$mu)) wstar <- wstar - object$noise$mu 
    if (method == "cond") 
      res <- condresC(wstar, object$noise$phi, object$noise$theta, TRUE)
    else res <- gresC(wstar, object$noise$phi, object$noise$theta)
    as.vector(res)
  }
  
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
  
  res <- noise.tfm(object)
  if (!is.null(object$noise$mu)) res <- res - object$noise$mu 
  if (method == "cond") {
    res <- as.numeric(condresC(res, object$noise$phi, object$noise$theta, TRUE))
  } else{
    res <- as.numeric(exactresC(res, object$noise$phi, object$noise$theta))
  }
  
  
  X <- cbind(b, g, se, z, p)
  colnames(X) <- c("Estimate", "Gradient", "Std. Error", "z Value", "Pr(>|z|)")
  rownames(X) <- names(b)
  
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
              z = object$output, aic = aic, bic = bic, gradient = g, 
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
print.summary.tfm <- function(x, stats = TRUE, 
                              digits = max(3L, getOption("digits") - 3L), ...) {
  print.summary.um(x, stats, digits, ...)
}


#' @rdname outliers
#' @param y an object of class \code{ts}
#' @export
outliers.tfm <- function(mdl, y = NULL, dates = NULL, c = 3, calendar = FALSE, 
                         easter = FALSE, n.ahead = NULL, p.value = 1, ...) {
  
  if (is.null(n.ahead)) n.ahead = 0
  else n.ahead <- abs(n.ahead)
  N <- noise.tfm(mdl, y, diff = FALSE)
  start <- start(N)
  freq <- frequency(N)
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
  
  bc <- mdl$noise$bc
  mdl$noise$bc <- FALSE
  if (is.null(mdl$noise$mu)) mu <- 0
  else mu <- mdl$noise$mu
  
  if (calendar||easter) {
    if (calendar)
      tfm1 <- calendar.um(mdl$noise, N, easter = easter, n.ahead = n.ahead)
    else
      tfm1 <- easter.um(mdl$noise, N, n.ahead)
    N <- noise.tfm(tfm1, diff = FALSE)
    A <- outliersC(N, FALSE,  mu, tfm1$noise$phi, tfm1$noise$nabla, tfm1$noise$theta, indx, abs(c))
  } else { 
    A <- outliersC(N, 1,  mu, mdl$noise$phi, mdl$noise$nabla, mdl$noise$theta, indx, abs(c))
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
  
  n <- length(N)
  if (!is.null(n.ahead)) n <- n + n.ahead
  tfi <- lapply(1:nrow(df), function(i) {
    p <- double(n)
    p[df[i, 1]] <- 1
    p <- ts(p, start = start, frequency = freq)
    prefix <- paste(df[i, 3], df[i, 2], sep = "")
    if (df[i, 3] == "IO") tf(p, w0 = df[i, 4], ma = mdl$noise$ma, ar = c(mdl$noise$i, mdl$noise$ar), par.prefix = prefix)
    else if (df[i, 3] == "AO") tf(p, w0 = df[i, 4], par.prefix = prefix)
    else if (df[i, 3] == "TC") tf(p, w0 = df[i, 4], ar = 1, par.prefix = prefix)
    else if (df[i, 3] == "LS") {
      p <- cumsum(p)
      p <- ts(p, start = start, frequency = freq)
      tf(p, w0 = df[i, 4], par.prefix = prefix)
    } else stop("unkown input")
  })
  
  if (calendar||easter) {
    if (is.null(mdl$xreg)) xreg <- tfm1$xreg
    else xreg <- merge(mdl$xreg, tfm1$xreg)
  } else xreg <- mdl$xreg
  
  if (is.null(tfi)) tfi <- mdl$inputs
  else if (!is.null(mdl$inputs)) tfi <- list.tf(mdl$inputs, tfi)
  
  mdl$noise$bc <- bc
  tfm1 <- tfm(xreg = xreg, inputs = tfi, noise = mdl$noise, fit = TRUE)
  if (p.value < 0.999) 
    tfm1 <- varsel.tfm(tfm1, p.value = p.value)
  
  return(tfm1)
  
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
#' @param ... additional arguments.
#'
#' @export predict.tfm
#' @export
predict.tfm <- function(object, y = NULL, ori = NULL, n.ahead = NULL, 
                        level = 0.95, ...) {
  stopifnot(is.tfm(object))
  y <- output.tfm(object, y)
  z <- noise.tfm(object, y, FALSE)  
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
  
  if (object$kx > 0) {
    if (nrow(object$xreg) < n) 
      stop("insufficient number of forecasts for input")
    z[t] <- z[t] + as.matrix(object$xreg[t, ]) %*% unlist(object$param[1:object$kx])
  }
  
  if (object$k > 0) {
    for (i in 1:object$k) {
      start1 <- start(object$inputs[[i]]$x)
      if (length(object$inputs[[i]]$x) < n) {
        if (has.um.tf(object$inputs[[i]])) {
          end1 <- end(object$inputs[[i]]$x)
          nahead <- (end[1] - end1[1])*s + end[2] - end1[2]
          object$inputs[[i]] <- 
            predict.tf(object$inputs[[i]], n.ahead = nahead)
          v <- var.predict.tf(object$inputs[[i]], n.ahead)
          X[t, 4] <- X[t, 4] + v
        } else stop("insufficient number of forecasts for input")
      }
      x <- object$inputs[[i]]$x
      if (object$inputs[[i]]$um$bc) x <- log(x)
      x <- filterC(x, object$inputs[[i]]$theta,
                   object$inputs[[i]]$phi, object$inputs[[i]]$delay)
      x <- ts(as.numeric(x), start = start1, frequency = s)
      x <- window(x, start = start, end = end) 
      z[t] <- z[t] + x[t] 
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
  print(c(param.tfm(x), sig2 = x$noise$sig2), ...)
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
#' @param ... additional arguments.
#'       
#' @return A "ts" object.
#' 
#' @export
residuals.tfm <- function(object, y = NULL, method = c("exact", "cond"), ...) {
  
  stopifnot(inherits(object, "tfm"))
  method <- match.arg(method)
  
  w <- noise.tfm(object, y)
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
#' @param ... additional arguments.
#' 
#' @return A "ts" object.
#' @export
signal <- function (mdl, ...) { UseMethod("signal") }

#' @rdname signal
#' @export
signal.tfm <- function(mdl, y = NULL, diff = TRUE, ...) {
  y <- output.tfm(mdl, y)
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
        if (t0 > 1 || t1 > N) 
          y <- y + diffC(x[t0:(t0+N-1)], mdl$noise$nabla, mdl$inputs[[i]]$um$bc)
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
        if (t0 > 1 || t1 > N)
          y <- y + x[t0:(t0+N-1)]
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
#' @param p.value probability value to decide whether or not to omit an input.
#' @param ... additional arguments.
#' @return A \code{tfm} object or a "um" if no input is significant at that level.
#' @export
varsel <- function (tfm, ...) { UseMethod("varsel") }

#' @rdname varsel
#' @export
varsel.tfm <- function(tfm, p.value = 0.10, ...) {
  p <- summary.tfm(tfm, p.value = TRUE)
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
  
  tfm(xreg = xreg, inputs = tfi, noise = tfm$noise)
  
}






is.tfm <- function(tfm) {
  inherits(tfm, "tfm")
}


param.tfm <- function(tfm) {
  unlist(tfm$param, use.names = TRUE)
}

output.tfm <- function(tfm, y = NULL) {
  if (is.null(y)) {
    if (is.null(tfm$output)) stop("argment y required")
    y <- eval(parse(text = tfm$output), .GlobalEnv)
  } else {
    if (!is.numeric(y)) stop("output must be a numeric vector")
  }
  
  if (!is.ts(y)) y <- ts(y)
  y
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
#' @param ... additional arguments.
#'
#'
ccf.tfm <- function(tfm, lag.max = NULL, method = c("exact", "cond"), ...) {

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

  u <- residuals.tfm(tfm, method = method)
  for (i in j) {
    end <- end(tfm$inputs[[i]]$x)
    s <- frequency(tfm$inputs[[i]]$x)
    x <- tfm$inputs[[i]]$x[tfm$inputs[[i]]$n.back:length(tfm$inputs[[i]]$x)]
    x <- ts(x, end = end, frequency = s)
    x <- residuals.um(tfm$inputs[[i]]$um, x, method)
    x <- window(x, start = start(u), end = end(u))
    pccf(x, u, lag.max = lag.max)
  }
  
  invisible(NULL)
  
}


#' @rdname ucomp
#' @param y an object of class \code{\link{ts}}.
#' @export
ucomp.tfm <- function(mdl, y = NULL, 
                      method = c("mixed", "forecast", "backcast"), ...) {

  y <- noise.tfm(mdl, y, FALSE)  
  mdl$noise$bc <- FALSE
  uc <- ucomp.um(mdl$noise, y, method)
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




