## tfarima/R/tfm.R
## Jose L Gallego (UC)

#' Transfer Function Model Constructor
#'
#' Creates and optionally fits a multiple-input transfer function model. 
#' A transfer function model relates an output time series to one or more 
#' input series (transfer functions), exogenous regressors, and a noise model.
#'
#' @param output A numeric vector or \code{ts} object representing the 
#'   dependent (output) time series. If \code{NULL}, it is taken from 
#'   \code{noise$z}.
#' @param xreg A numeric matrix or \code{ts} object of exogenous regressors. 
#'   Columns correspond to different regressors. Defaults to \code{NULL}.
#' @param inputs A list of transfer function objects of class \code{tf}. 
#'   Each element represents one stochastic input. Can also be a single 
#'   \code{tf} object. Defaults to \code{NULL}.
#' @param noise An object of class \code{um} describing the univariate 
#'   noise model. This defines the ARIMA-type structure for the residuals.
#' @param fit Logical. If \code{TRUE} (default), the model parameters are 
#'   estimated by maximum likelihood after construction.
#' @param new.name Logical. Internal use. If \code{TRUE} (default), a new 
#'   name is assigned to the output series. Otherwise, the name stored in 
#'   \code{noise$z} is preserved.
#' @param envir Environment in which the function arguments are evaluated. 
#'   If \code{NULL}, the calling environment is used.   
#' @param ... Additional arguments passed to \code{\link{fit.tfm}} when 
#'   \code{fit = TRUE}.
#'
#' @details
#' All series must have the same frequency. Input series must span at least 
#' the same period as output. The function applies differencing and Box-Cox 
#' transformation as specified in \code{noise}.
#'
#' @return Object of class \code{tfm} with components: output, xreg, inputs, 
#'   noise, param, kx, k, optim, method, and call.
#'
#' @seealso \code{\link{tf}}, \code{\link{um}}, \code{\link{fit.tfm}}
#'
#' @references
#' Box, G. E., Jenkins, G. M., Reinsel, G. C., & Ljung, G. M. (2015).
#' \emph{Time Series Analysis: Forecasting and Control} (5th ed.). Wiley.
#'
#' @examples
#' \dontrun{
#' Y <- seriesJ$Y - mean(seriesJ$Y)
#' X <- seriesJ$X - mean(seriesJ$X)
#' umx <- um(X, ar = 3)
#' umy <- fit(umx, Y)
#' tfx <- tfest(Y, X, delay = 3, p = 2, q = 2, um.x = umx, um.y = umy)
#' tfmy <- tfm(Y, inputs = tfx, noise = um(ar = 2))
#' }
#'
#' @export
tfm <- function(output = NULL, xreg = NULL, inputs = NULL, noise, fit = TRUE,
                new.name = TRUE, envir = parent.frame (), ...) {
  call <- match.call()
  stopifnot(is.um(noise))

  if (!is.null(inputs)) {
    if (inherits(inputs, "tf")) {
      inputs <- list(inputs)
    } else {
      inputs <- inputs[!sapply(inputs, is.null)]
      if (any(!sapply(inputs, is.tf)))
        stop("inputs must be 'tf' objects")
    }
    k <- length(inputs)
  } else k <- 0
  
  if (is.null(output)) {
    if (is.null(noise$z)) stop("missing output series")
    else {
      if (exists(noise$z, envir = envir)) {
        output <- get(noise$z, envir = envir)
      } else stop("Output not found in environment: ", output)    
    }
  } else if (is.character(output)) {
    noise$z <- output
    if (exists(output, envir = envir)) {
      output <- get(output, envir = envir)
    } else stop("Output not found in environment: ", output)    
  } else if (is.numeric(output)) {
    if (new.name) noise$z <- deparse(substitute(output))
  } else {
    stop("Invalid output type: must be ts, character or numeric")
  }
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
    param <- unlist(unname(param))
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
      else cn <- paste0(xn, 1:kx)
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
      if (length(x) > n) x[1:n] else x
    })
    reg <- tryCatch(
      lm(y ~ X - 1),
      error = function(e) {
        stop(sprintf("Failed to estimate 'xreg' parameters: %s", e$message))
      }
    )    
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
  
  if (fit) {
    mod <- tryCatch(
      fit.tfm(mod, envir = envir, ...),
      error = function(e) {
        warning("Model fitting failed: ", e$message)
        mod
      }
    )
  }  

  return(mod)

}

#' Fit a Transfer Function Model
#'
#' Estimates the parameters of a transfer function model of class \code{tfm} by
#' (conditional or exact) maximum likelihood.
#'
#' @param mdl An object of class \code{tfm} created with \code{\link{tfm}}.
#' @param y Optional \code{ts} object containing the output series. If
#'   \code{NULL}, the output stored in \code{noise} is used.
#' @param method Character string specifying likelihood method: "exact" for
#'   exact maximum likelihood or "cond" for conditional maximum likelihood.
#'   Default is "exact".
#' @param optim.method Character. Optimization method passed to
#'   \code{\link{optim}}. Default is \code{"BFGS"}. Other options:
#'   "Nelder-Mead", "CG", "L-BFGS-B".
#' @param show.iter Logical. If \code{TRUE}, prints iteration progress of the
#'   likelihood optimization.
#' @param fit.noise Logical. If \code{TRUE} (default), the parameters of the
#'   noise model are estimated. If \code{FALSE}, noise parameters are fixed at
#'   their current values.
#' @param envir Environment in which the function arguments are evaluated. If
#'   \code{NULL}, the calling environment is used.
#' @param ... Additional arguments.
#'
#' @return An updated object of class \code{tfm} containing fitted parameters,
#'   estimated innovation variance, and optimization details.
#'
#' @seealso \code{\link{tfm}}
#'
#' @examples
#' \dontrun{
#' Y <- seriesJ$Y - mean(seriesJ$Y)
#' X <- seriesJ$X - mean(seriesJ$X)
#' umx <- um(X, ar = 3)
#' umy <- fit(umx, Y)
#' tfx <- tfest(Y, X, delay = 3, p = 2, q = 2, um.x = umx, um.y = umy)
#' tfmy <- tfm(Y, inputs = tfx, noise = um(ar = 2), fit = FALSE)
#' tfmy_fit <- fit(tfmy)
#' }
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
  if (!is.invertible.um(mdl$noise) && !is.null(mdl$noise$param))
    stop("Non-invertible MA preestimates")
  if (length(method) == 2 && !is.null(mdl$noise$method)) 
    method <- mdl$noise$method
  method <- match.arg(method)
  mdl$noise$method <- method
  exact <- method == "exact"
  
  noiseTFM <- function(par) {
    mdl <- .update_tfm(mdl, par)
    if (is.null(mdl$noise$mu)) wstar <- w
    else wstar <- w - mdl$noise$mu
    if (mdl$kx > 0) {
      wstar <- wstar - xreg %*% unlist(mdl$param[1:mdl$kx])
    }
    if (mdl$k > 0) {
      for (i in 1:mdl$k) {
        x <- filterC(X[[i]], mdl$inputs[[i]]$theta,
                     mdl$inputs[[i]]$phi, mdl$inputs[[i]]$delay)
        if (t0[i] > 1 || t1[i] > 1) wstar <- wstar - x[t0[i]:t1[i]]
        else wstar <- wstar - x
      }
    }
    list(mdl = mdl, wstar = wstar)
  }
  
  logLikTFM <- function(par) {
    lst <- noiseTFM(par)
    mdl <- lst$mdl
    wstar <- lst$wstar
    if (mdl$noise$is.adm) {
      if (exact) ll <- -ellarmaC(wstar, mdl$noise$phi, mdl$noise$theta)
      else ll <- -cllarmaC(wstar, mdl$noise$phi, mdl$noise$theta)
    } else {
      ll <- ll0
    }
    if (show.iter) print(c(loglik = ll, par))
    return(ll)
  }
  
  y <- .output_tfm(mdl, y, envir)
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
      if (length(x) > n)
        t1[i] <<- tf$t.start + n - 1
      i <<- i + 1
      x
    })
  }
  
  b0 <- param.tfm(mdl)
  ll0 <- logLikTFM(b0)
  opt <- optim(b0, logLikTFM, method = optim.method, ...)
  if(opt$convergence > 0)
    warning(gettextf("possible convergence problem: optim gave code = %d",
                     opt$convergence), domain = NA)
  
  mdl <- .update_tfm(mdl, opt$par)
  wstar <- noise.tfm(mdl, y, diff = TRUE, envir=envir)
  if (!is.null(mdl$noise$mu)) wstar <- wstar - mdl$noise$mu
  if (method == "cond")
    res <- condresC(wstar, mdl$noise$phi, mdl$noise$theta, TRUE)
  else res <- gresC(wstar, mdl$noise$phi, mdl$noise$theta)
  mdl$noise$sig2 <- sum(res^2)/length(res)
  mdl$optim <- mdl$noise$optim <- opt
  mdl$noise$method <- method
  
  if (!fit.noise) names(mdl$noise$param) <- par.noise 
  
  return(mdl)
  
}

#' Summarize Transfer Function Model
#' 
#' Produces summary statistics for a fitted transfer function model including
#' parameter estimates, standard errors, and diagnostic tests.
#' 
#' @param object A fitted \code{tfm} object.
#' @param y Optional \code{ts} object for alternative output series.
#' @param method Character: "exact" or "cond" for residual calculation.
#' @param digits Number of significant digits for printing.
#' @param envir Environment for evaluation. NULL uses calling environment.
#' @param ... Additional arguments:
#'   \code{p.values} (logical) returns only p-values;
#'   \code{table} (logical) returns only coefficient table.
#'
#' @details
#' Computes parameter estimates with standard errors (from Jacobian),
#' z-statistics, p-values, AIC, BIC, log-likelihood, Ljung-Box tests
#' (at lags p+q+1 and n/4+p+q), and Bartlett heteroscedasticity test.
#'
#' @return Object of class \code{summary.tfm} containing: call, coefficient
#'   table, variance-covariance matrix, residuals, diagnostic statistics,
#'   information criteria, and time series attributes.
#'
#' @seealso \code{\link{print.summary.tfm}}
#'
#' @examples
#' \dontrun{
#' Y <- seriesJ$Y - mean(seriesJ$Y)
#' X <- seriesJ$X - mean(seriesJ$X)
#' umx <- um(X, ar = 3)
#' umy <- fit(umx, Y)
#' tfx <- tfest(Y, X, delay = 3, p = 2, q = 2, um.x = umx, um.y = umy)
#' tfmy <- tfm(Y, inputs = tfx, noise = um(ar = 2))
#' sm <- summary(tfmy)
#' print(sm)
#' }
#'
#' @export
summary.tfm <- function(object, y = NULL, method = c("exact", "cond"),
                        digits = max(3L, getOption("digits") - 3L), 
                        envir = parent.frame(), ...) {
  
  stopifnot(inherits(object, "tfm"))
  model.name <- deparse(substitute(object))
  
  args <- list(...)
  if(is.null(args[["p.values"]])) p.values <- FALSE
  else p.values <- args[["p.values"]]
  if(is.null(args[["table"]])) table <- FALSE
  else table <- args[["table"]]
  
  method <- match.arg(method)
  if (!is.null(object$noise$method)) method <- object$noise$method
  
  resTFM <- function(b) {
    mdl <- .update_tfm(object, b)
    wstar <- noise.tfm(mdl, y, diff = TRUE, envir=envir)
    if (!is.null(mdl$noise$mu)) wstar <- wstar - mdl$noise$mu
    if (method == "cond")
      res <- condresC(wstar, mdl$noise$phi, mdl$noise$theta, TRUE)
    else res <- gresC(wstar, mdl$noise$phi, mdl$noise$theta)
    as.vector(res)
  }
  
  y <- .output_tfm(object, y, envir)
  w <- diffC(y, object$noise$nabla, object$noise$bc)
  N <- length(y)
  n <- length(w)
  b <- param.tfm(object)
  
  ll <- logLik(object, y, method, envir)
  aic <- (-2.0*ll+2*length(b))/n
  bic <- (-2*ll+log(n)*length(b))/n
  
  res <- resTFM(b)
  J <- numDeriv::jacobian(resTFM, b)
  g <- t(J) %*% res
  ssr <- sum(res^2)
  object$noise$sig2 <- ssr/length(res)
  if (!is.null(object$optim$hessian)) H <- object$optim$hessian
  else H <- t(J) %*% J
  varb <- tryCatch ({ 
    solve(H)*object$noise$sig2
  }, error = function(e) {
    message("Error inverting matrix: ", e$message)
    return(NULL)
  })
  b <- param.tfm(object)
  if (is.null(varb)) {
    warning("Could not compute standard errors: Hessian is singular")
    se <- rep(NA_real_, length(b))
    z <- rep(NA_real_, length(b))
    p <- rep(NA_real_, length(b))
  } else {
    se <- sqrt(diag(varb))
    z <- b/se
    p <- pnorm(abs(z), lower.tail = FALSE)*2
  }  

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

#' Print Summary of Transfer Function Model
#'
#' Print method for objects of class \code{summary.tfm}.
#'
#' @param x A \code{summary.tfm} object.
#' @param stats Logical. If TRUE, prints diagnostic statistics.
#' @param short Logical. If TRUE, prints abbreviated output.
#' @param digits Number of significant digits.
#' @param ... Additional arguments.
#'
#' @seealso \code{\link{summary.tfm}}
#'
#' @export
print.summary.tfm <- function(x, stats = TRUE, short = FALSE,
                              digits = max(3L, getOption("digits") - 3L), ...) {
  print.summary.um(x, stats, digits, short = short, ...)
}

#' Print Transfer Function Model
#'
#' Print method for objects of class \code{tfm}.
#'
#' @param x A \code{tfm} object.
#' @param ... Additional arguments passed to \code{\link{summary.tfm}} 
#'   and \code{\link{print.summary.tfm}}.
#'
#' @details
#' Prints a summary of the transfer function model by calling 
#' \code{summary(x)} with \code{short = TRUE}.
#'
#' @seealso \code{\link{tfm}}, \code{\link{summary.tfm}}
#'
#' @export
print.tfm <- function(x, ...) {
  print(summary(x), short = TRUE, ...)
}

#' Signal component of a TF model
#'
#' \code{signal} extracts the signal of a TF model.
#'
#' @param mdl an object of the class \code{tfm}.
#' @param y output of the TF model if it is different to that of the \code{tfm}
#'   object.
#' @param diff logical. If TRUE, the signal is differenced with the "i" operator
#'   of the univariate model of the noise.
#' @param type Character vector specifying signal components to extract: "xreg"
#'   for exogenous regressors only,"inputs" for transfer function inputs only.
#'   If both provided (default), includes all components.
#' @param envir Environment in which the function arguments are evaluated. By
#'   default, the calling environment is used.
#' @param ... additional arguments.
#'
#' @return A "ts" object.
#' @export
signal <- function (mdl, ...) { UseMethod("signal") }

#' @rdname signal
#' @export
signal.tfm <- function(mdl, y = NULL, diff = FALSE, type = c("xreg", "inputs"),
                       envir = parent.frame(), ...) {
  stopifnot(is.tfm(mdl))
  type <- match.arg(type, several.ok = TRUE)
  y <- .output_tfm(mdl, y, envir=envir)
  start <- start(y)
  end <- end(y)
  s <- frequency(y)
  N <- length(y)
  y <- y*0
  if ("xreg" %in% type) {
    if (mdl$kx > 0) {
      if (nrow(mdl$xreg) > N)
        y <- y + as.matrix(mdl$xreg[1:N, ]) %*% unlist(mdl$param[1:mdl$kx])
      else
        y <- y + mdl$xreg %*% unlist(mdl$param[1:mdl$kx])
    }
  }
  
  if ("inputs" %in% type) {
    if (mdl$k > 0) {
      for (i in 1:mdl$k) {
        x <- filterC(mdl$inputs[[i]]$x, mdl$inputs[[i]]$theta,
                     mdl$inputs[[i]]$phi, mdl$inputs[[i]]$delay)
        t0 <- mdl$inputs[[i]]$t.start
        t1 <- mdl$inputs[[i]]$t.end
        if (t0 > 1 || length(mdl$inputs[[i]]$x) > t1) y <- y + x[t0:t1]
        else y <- y + x
      }
    }
  }
  if (diff && N > mdl$noise$d) y <- diffC(y, mdl$noise$nabla, FALSE)
  
  ts(y, end = end, frequency = s)
}

#' Extract Noise Component from Transfer Function Model
#' 
#' Computes the noise series (output minus fitted signal) from a transfer
#' function model.
#'
#' @param mdl A \code{tfm} object.
#' @param y Optional \code{ts} object for alternative output series.
#' @param diff Logical. If TRUE (default), returns differenced noise series
#'   (stationary). If FALSE, returns noise in original scale.
#' @param exp Logical. If TRUE, applies exponential transformation (inverse
#'   of log). Only relevant when \code{diff = FALSE} and Box-Cox transformation
#'   was used (\code{bc = TRUE}).
#' @param envir Environment for evaluation. NULL uses calling environment.
#' @param ... Additional arguments.
#'
#' @return A \code{ts} object containing the noise series, computed as
#'   output minus all transfer function and regressor effects.
#' 
#' @details
#' The noise represents the component of the output not explained by the
#' transfer functions and exogenous regressors. When \code{diff = TRUE},
#' the differencing operator from the noise model is applied, resulting in
#' a stationary series suitable for ARMA modeling.
#'
#' @seealso \code{\link{signal.tfm}}, \code{\link{residuals.tfm}}, \code{\link{tfm}}
#' 
#' @export
noise <- function (mdl, ...) { UseMethod("noise") }

#' @rdname noise
#' @export
noise.tfm <- function(mdl, y = NULL, diff = TRUE, exp = FALSE,
                      envir = parent.frame(), ...) {
  stopifnot(is.tfm(mdl))
  y <- .output_tfm(mdl, y, envir = envir)
  if (diff && exp) {
    warning("'exp' is ignored when 'diff = TRUE'", call. = FALSE)
  }  
  start <- start(y)
  end <- end(y)
  s <- frequency(y)
  N <- length(y)
  if (mdl$noise$bc) y <- log(y)
  
  if (diff && mdl$noise$d > 0) {
    y <- diffC(y, mdl$noise$nabla, FALSE)
    n <- length(y)
    if (mdl$kx > 0) {
      for (i in 1:mdl$kx) {
        x <- diffC(mdl$xreg[, i], mdl$noise$nabla, FALSE)
        if (length(x) > n) x <- x[1:n]
        y <- y - x*mdl$param[[i]]
      }
    }
    if (mdl$k > 0) {
      for (i in 1:mdl$k) {
        x <- diffC(mdl$inputs[[i]]$x, mdl$noise$nabla, mdl$inputs[[i]]$um$bc)
        x <- filterC(x, mdl$inputs[[i]]$theta,
                     mdl$inputs[[i]]$phi, mdl$inputs[[i]]$delay)
        t0 <- mdl$inputs[[i]]$t.start
        t1 <- mdl$inputs[[i]]$t.end
        y <- y - x[t0:(t0+n-1)]
      }
    }
  }
  else {
    if (mdl$kx > 0) {
      if (nrow(mdl$xreg) > N) {
        y <- y - as.matrix(mdl$xreg[1:N, ]) %*% unlist(mdl$param[1:mdl$kx])
      } else
        y <- y - mdl$xreg %*% unlist(mdl$param[1:mdl$kx])
    }
    if (mdl$k > 0) {
      for (i in 1:mdl$k) {
        if (mdl$inputs[[i]]$um$bc) {
          x <- filterC(log(mdl$inputs[[i]]$x), mdl$inputs[[i]]$theta,
                       mdl$inputs[[i]]$phi, mdl$inputs[[i]]$delay)
        } else {
          x <- filterC(mdl$inputs[[i]]$x, mdl$inputs[[i]]$theta,
                       mdl$inputs[[i]]$phi, mdl$inputs[[i]]$delay)
        }
        t0 <- mdl$inputs[[i]]$t.start
        t1 <- mdl$inputs[[i]]$t.end
        y <- y - x[t0:t1]
      }
    }
    if (mdl$noise$bc && exp) y <- exp(y)
  }
  ts(as.vector(y), end = end, frequency = s)
}

#' Extract Residuals from Transfer Function Model
#'
#' Computes exact or conditional residuals from a fitted transfer function model.
#'
#' @param object A fitted \code{tfm} object.
#' @param y Optional \code{ts} object for alternative output series.
#' @param method Character: "exact" estimates presample values  to
#'   compute residuals; "cond" fixes presample values at zero.
#' @param envir Environment for evaluation. NULL uses calling environment.
#' @param ... Currently unused.
#'
#' @return A \code{ts} object containing model residuals with the same
#'   time series attributes as the output series.
#'
#' @export
residuals.tfm <- function(object, y = NULL, method = c("exact", "cond"),
                          envir=parent.frame(), ...) {
  
  stopifnot(inherits(object, "tfm"))
  method <- match.arg(method)
  
  w <- noise.tfm(object, y, envir=envir)
  if (!is.null(object$noise$mu)) w <- w - object$noise$mu
  if (method == "cond")
    res <- condresC(w, object$noise$phi, object$noise$theta, TRUE)
  else 
    res <- exactresC(w, object$noise$phi, object$noise$theta)
  res <- as.vector(res)
  if (is.ts(w))
    res <- ts(res, end = end(w), frequency = frequency(w))
  else
    res <- res
  
  return(res)
  
}

#' Diagnostic Checking for Transfer Function Models
#'
#' Produces diagnostic plots for residuals of a fitted transfer function model.
#'
#' @param mdl A fitted \code{tfm} object.
#' @param y Optional \code{ts} object for alternative output series.
#' @param method Character: "exact" or "cond" for residual calculation.
#' @param lag.max Maximum lag for ACF/PACF plots.
#' @param lags.at Specific lags to display in ACF/PACF.
#' @param freq.at Specific frequencies for cumulative periodogram.
#' @param std Logical. If TRUE, standardizes residuals.
#' @param envir Environment for evaluation. NULL uses calling environment.
#' @param ... Additional arguments passed to \code{\link{ide}}.
#'
#' @details
#' Generates five diagnostic plots: time series plot of residuals, 
#' histogram, ACF, PACF, and cumulative periodogram. Uses the \code{\link{ide}}
#' function for plotting.
#'
#' @seealso \code{\link{tfm}}, \code{\link{tsdiag.tfm}}
#'
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

#' Cross-correlation check
#'
#' \code{ccf} displays ccf between prewhitened inputs and residuals.
#'
#' @param x a \code{tfm} object.
#' @param lag.max number of lags.
#' @param method Exact/conditional residuals.
#' @param envir environment in which the function arguments are evaluated.
#' @param ... additional arguments.
#'
#' @export
ccf.tfm <- function(x, lag.max = NULL, method = c("exact", "cond"), 
                    envir = parent.frame (), ...) {
  stopifnot(is.tfm(x))
  if (x$k < 1) 
    stop("Model has no transfer function inputs")
  j <- c()
  for (i in 1:x$k) 
    if (x$inputs[[i]]$um$k > 0)
      j <- c(j, i)
  
  k <- length(j)
  if (k < 1) 
    stop("Model has no transfer function inputs")
  
  if (k > 1) {
    tryCatch({
      oldpar <- par(mfrow = c(k, 1), mar = c(3, 4, 2, 1))
      on.exit(par(oldpar))
    }, error = function(e) {
      stop(
        "Graphics device too small for ", k, " plots.\n",
        "Resize window or use: dev.new(width=7, height=", 3*k, ")"
      )
    })
  }  
  u <- residuals.tfm(x, method = method, envir = envir)
  for (i in j) {
    end <- end(x$inputs[[i]]$x)
    s <- frequency(x$inputs[[i]]$x)
    a <- x$inputs[[i]]$x[x$inputs[[i]]$n.back:length(x$inputs[[i]]$x)]
    a <- ts(a, end = end, frequency = s)
    a <- residuals.um(x$inputs[[i]]$um, a, method, envir=envir)
    a <- window(a, start = start(u), end = end(u))
    pccf(a, u, lag.max = lag.max)
  }
  
  invisible(NULL)
}

#' Forecast Transfer Function Model
#'
#' Computes point forecasts and prediction intervals for transfer function models.
#' 
#' @param object A fitted \code{tfm} object.
#' @param newdata Optional matrix or vector of future values for exogenous
#'   regressors and inputs. Rows correspond to forecast horizon, columns to
#'   predictors.
#' @param y Optional \code{ts} object for alternative output series.
#' @param ori Forecast origin (observation index). Default is last observation.
#' @param n.ahead Number of steps ahead to forecast. Default is series frequency.
#' @param level Confidence level(s) for prediction intervals (0-1). Default is 0.95.
#'   Can be a vector for multiple intervals.
#' @param i Optional differencing operator (lagpol) to apply before forecasting.
#' @param envir Environment for evaluation. NULL uses calling environment.
#' @param ... Additional arguments (currently unused).
#' 
#' @details 
#' Future values for transfer function inputs can be provided in three ways:
#' (1) extending input series beyond output length, (2) automatic forecasting
#' from associated \code{um} models, or (3) via the \code{newdata} argument.
#' 
#' If Box-Cox transformation was used, forecasts are back-transformed and
#' intervals adjusted accordingly.
#'
#' @return Object of class \code{predict.tfm} containing:
#' \item{z}{Complete series including forecasts}
#' \item{rmse}{Root mean square error for each forecast}
#' \item{low, upp}{Lower and upper prediction interval bounds (matrices)}
#' \item{level}{Confidence level(s) used}
#' \item{dates}{Time points for all observations}
#' \item{ori, ori.date}{Forecast origin (index and date)}
#' \item{n.ahead}{Number of forecasts}
#'
#' @seealso \code{\link{tfm}}, \code{\link{fit.tfm}}
#'
#' @export
predict.tfm <- function(object, newdata=NULL, y = NULL, ori = NULL, n.ahead = NULL,
                        level = 0.95, i = NULL,  envir=NULL, ...) {
  stopifnot(is.tfm(object))
  if (is.null (envir)) envir <- parent.frame ()
  y <- .output_tfm(object, y, envir=envir)
  if (!is.null(i)) {
    i <- lagpol0(i, "i", envir=envir)
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
    Xf <- NULL
    if (!is.null(newdata)) {
      nms <- colnames(object$xreg)
      nm1 <- colnames(newdata)
      if (all(nms %in% nm1)) {
        Xf <- as.matrix(newdata[1:n.ahead, nms, drop = FALSE])
      } else if(any(nms %in% nm1)) {
        if (nrow(object$xreg) < n) 
          stop("insufficient number of forecasts for input") 
        Xf <- as.matrix(object$xreg[t, , drop = FALSE])
        Xf[, nms %in% nm1] <- newdata[1:n.ahead, nms[nms %in% nm1]]
      }
    }  
    
    if (is.null(Xf)) {
      if (nrow(object$xreg) < n) 
        stop("insufficient number of forecasts for input") 
      Xf <- as.matrix(object$xreg[t, ])  
    }
    
    z[t] <- z[t] +  Xf %*% unlist(object$param[1:object$kx])
  }

  if (object$k > 0) {
    for (i in 1:object$k) {
      start1 <- start(object$inputs[[i]]$x)
      if ( any(colnames(newdata) %in% object$inputs[[i]]$x.name) ) {
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
  z[1:ori] <- y[1:ori]
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

#' @rdname modify
#' @export
modify.tfm <- function(mdl, ar = NULL, i = NULL, ma = NULL, mu = NULL, 
                       sig2 = NULL, bc = NULL, ...) {
  um1 <- modify.um(mdl$noise, ar = ar, i = i, ma = ma, bc = bc, sig2 = sig2,
                   fit = FALSE, ...)
  um1$z <- mdl$noise$z
  tfm(xreg = mdl$xreg, inputs = mdl$inputs, noise = um1, ...)
}

#' @rdname setinputs
#' 
#' @details
#' For \code{tfm} objects: If the model already has inputs of the same type,
#' new ones are appended (combined). The model is re-fitted by default unless
#' \code{fit = FALSE}.
#'
#' @export
setinputs.tfm <- function(mdl, xreg = NULL, inputs = NULL, y = NULL, 
                          envir = parent.frame (), ...) {
  stopifnot(is.tfm(mdl))
  if (!is.null(y)) mdl$noise$z <- deparse(substitute(y))
  y <- .output_tfm(mdl, y, envir)
  if (!is.null(xreg)) {
    xreg <- as.matrix(xreg)
    if (mdl$kx > 0) {
      stopifnot(nrow(mdl$xreg) == nrow(xreg))
      nms <- colnames(mdl$xreg)
      nms1 <- colnames(xreg)
      if (any(nms1 %in% nms))
        stop("Input name in use")
      xreg <- cbind(mdl$xreg, xreg)
      colnames(xreg) <- c(nms, nms1)
    }
  } else xreg <- mdl$xreg
  
  if (!is.null(inputs) && mdl$k > 0) {
    if (is.tf(inputs)) inputs <- list(inputs)
    if (!all(sapply(inputs, is.tf)))
      stop("'inputs' must be a 'tf' object or list of 'tf' objects")
    inputs <- c(mdl$inputs, inputs)  
  } else if(mdl$k > 0) inputs <- mdl$inputs
  tfm(output = y, xreg = xreg, inputs = inputs, noise = mdl$noise, 
      new.name = FALSE, envir = envir, ...)
}

#' @rdname calendar
#' @export
calendar.tfm <- function(mdl, y = NULL, 
                         form = c("dif", "td", "td7", "td6", "wd"), 
                         ref = 0, lom = TRUE, lpyear = TRUE, easter = FALSE,
                         len = 4, easter.mon = FALSE, n.ahead = 0, p.value = 1, 
                         envir = parent.frame (), ...) {
  stopifnot(is.tfm(mdl))
  if (!is.null(y)) y.name <- deparse(substitute(y))
  else y.name <- NULL
  y <- .output_tfm(mdl, y, envir = envir)
  if (frequency(y) != 12) stop("function only implemented for monthly ts")
  
  n.ahead <- abs(n.ahead)
  xreg <- CalendarVar(y, form, ref, lom, lpyear, easter, len, easter.mon, n.ahead)
  tfm1 <- setinputs.tfm(mdl, xreg, y = y)
  if (!is.null(y.name)) mdl$noize$z <- y.name
  if (p.value > 0.999) return(tfm1)
  
  p <- summary.tfm(tfm1, p.values = TRUE)
  p <- (p[1:ncol(xreg)] <= p.value)
  if (all(p)) return(tfm1)
  if (any(p)) {
    xreg <- tfm1$xreg[, p, drop = FALSE]
    tfm1 <- tfm(y, xreg = xreg, inputs = tfm1$inputs, noise = tfm1$noise,
                envir = envir, new.name = FALSE, ...)
    return(tfm1)
  }
  return(mdl)
}

#' @rdname outliers
#' @param y an object of class \code{ts}, optional.
#' @export
outliers.tfm <- function(mdl, y = NULL, types = c("AO", "LS", "TC", "IO"), 
                         dates = NULL,  c = 3, calendar = FALSE, easter = FALSE, 
                         resid = c("exact", "cond"), n.ahead = 0, 
                         p.value = 1, tc.fix = TRUE, 
                         envir = parent.frame (), ...) {
  stopifnot(is.tfm(mdl))
  n.ahead <- abs(n.ahead)
  if (!is.null(y)) mdl$noise$z <- deparse(substitute(y))
  y <- .output_tfm(mdl, y, envir)
  if (mdl$kx > 0) {
    if (nrow(mdl$xreg) - length(y) < n.ahead)
      stop(paste0("'n.ahead' incompatible with xreg"))
    n.ahead <- nrow(mdl$xreg) - length(y)
  }
  N <- noise.tfm(mdl, y, diff = FALSE, envir=envir)
  start <- start(N)
  freq <- frequency(N)
  if( freq < 12) {
    calendar <- FALSE
    easter <- FALSE
  }
  if ((mdl$noise$p > 50 || mdl$noise$q > 50) && length(resid) > 1)
    resid <- "cond"
  resid <- match.arg(resid)
  eres <- resid == "exact"
  types <- match.arg(types, several.ok = TRUE)
  types <- c("AO", "LS", "TC", "IO") %in% types
  if (!is.null(dates)) {
    if (is.numeric(dates) && length(dates) == 2) dates <- list(dates)
    if (is.list(dates)) {
      indx <- sapply(dates, function(x) {
        if (length(x) != 2) stop("invalid date format")
        if (freq > 1) (x[1] - start[1] + 1)*freq - (start[2] - 1) - (freq - x[2])
        else (date[1] - start[1] + 1)
      })
      if (any(indx) < 1||indx > length(y)) stop("invalid date format")
    } else if (is.numeric(dates)) indx <- dates
    else indx <- 0
    indx <- unique(indx)
  }  else indx <- 0

  bc <- mdl$noise$bc
  mdl$noise$bc <- FALSE
  if (is.null(mdl$noise$mu)) mu <- 0
  else mu <- mdl$noise$mu
  
  tfm1 <- NULL
  if (calendar||easter) {
    if (calendar) tfm1 <- calendar.um(mdl$noise, N, easter = easter, 
                          n.ahead = n.ahead, envir=envir, ...)
    else tfm1 <- easter.um(mdl$noise, N, n.ahead, envir = envir, ...)
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
      prefix <- paste0(df[k, 3], df[k, 2])
      if (df[k, 3] == "IO") {
        tf(p, w0 = df[k, 4], ma = mdl$noise$ma,
           ar = c(mdl$noise$i, mdl$noise$ar), par.prefix = prefix)
      } else if (df[k, 3] == "TC") {
        tf(p, w0 = df[k, 4], ar = 1, par.prefix = prefix)
      } else stop("unkown input")
    })
  } else tfi <- NULL 
  
  nms <- c()
  xreg <- cbind()
  if (mdl$kx > 0) {
    nms <- c(nms, colnames(mdl$xreg))
    xreg <- cbind(xreg, mdl$xreg)
  }
  if (calendar||easter) {
    nms <- c(nms, colnames(tfm1$xreg))
    xreg <- cbind(xreg, tfm1$xreg)
  }
  if (any(i)) {
    nms <- c(nms, colnames(X))
    xreg <- cbind(xreg, X)
  }
  if (!is.null(nms)) colnames(xreg) <- nms

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
  if (!all(type %in% c("AO", "P", "LS", "S", "TC", "IO", "R")))
    stop("invalid intervention type")
  
  y <- .output_tfm(mdl, y, envir = envir)
  s <- frequency(y)
  
  N <- noise.tfm(mdl, y, diff = FALSE, exp = TRUE, envir = envir)
  k1 <- length(type)
  
  # convert time to list of dates
  if (!is.list(time)) {
    if (is.vector(time)) {
      if (s == 1) time <- lapply(1:nrow(time), function(x) time[x, ])
      else if (length(time) %% 2 == 0) {
        time <- matrix(time , ncol = 2)
      } else stop("invalid date format")
    } else if(is.matrix(time)) {
      x <- dim(time)
      if (x[1] != 2 && x[2] == 2) stop("invalid date format")
      if (x[1] == 2 && x[2] != 2) time <- t(time)
      else if(x[1] == 2 && x[2] != 2 && all(time[2, ] <= s)) time <- t(time)
    }
    time <- lapply(1:nrow(time), function(x) time[x, ])
  }
  k2 <- length(time)
  
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
    type1 <- switch(
      type[k],
      "AO" = , "P" = , "IO" = , "TC" = "P",
      "LS" = , "S" = "S",
      "R" = "R",
      stop("Unknown type: ", type[k])
    )
    InterventionVar(y, time[[k]], type1, n.ahead)
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
        x <- ts(x, end = end(N), frequency = frequency(N))
        if (type[k] == "TC") {
          tf <- tfest(N, x, p = 1, q = 0, um.y = mdl$noise, 
                      par.prefix = names[k], envir = envir)
          tf$x.name <- names[k]
        } else { # "IO"
          u <- residuals(mdl, y, envir = envir)
          w0 <- tsvalue(u, time[[k]])
          tf <- tf(x, w0 = w0, ma = mdl$noise$ma, 
                   ar = c(mdl$noise$i, mdl$noise$ar), par.prefix = names[k], 
                   envir = envir)
        }
        tfm1 <- setinputs(mdl, inputs = tf, ...)
      }
      if (type[k] == "TC") {
        s <- summary(tfm1, table = TRUE)
        d <- s[names(tf$param)[2], 1]
        name <- names[k]
        names[k] <<- paste0(names[k], "(", format(d, digits = 2), ")")
        s[name, ]
      } else summary(tfm1, table = TRUE)[names[k], ]
    })
    colnames(tbl) <- names
    tbl <- tbl[, order(tbl[5, ])]
    return(tbl)
  } else {
    if (all(i)) {
      return(setinputs(mdl, xreg = X, ...))
    } else {
      if (any(i)) {
        X1 <- X[, i, drop = FALSE]
        colnames(X1) <- names[i]
      } else {
        X1 <- NULL
      }
      ltf <- lapply( (1:k1)[!i], function(k) {
        x <- X[, k]
        x <- ts(x, start = start(y), frequency = frequency(y))
        if (type[k] == "TC") {
          tf <- tfest(N, x, p = 1, q = 0, um.y = mdl$noise, 
                      par.prefix = names[k]) 
          tf$x.name <- names[k]
        } else { # "IO"
          u <- residuals(mdl, y, envir = envir)
          tf <- tf(x, w0 = tsvalue(u, time[[k]]), ma = mdl$noise$ma,
                   ar = c(mdl$noise$i, mdl$noise$ar), par.prefix = names[k])
        }
        return(tf)
      })
      return(setinputs(mdl, xreg = X1, inputs = ltf, y = y, envir = envir, ...))
    }
  }
}

#' Variable selection
#'
#' \code{varsel} omits non-significant inputs from a transfer function model.
#'
#' @param tfm a \code{tfm} object.
#' @param y a "ts" object.
#' @param p.value probability value to decide whether or not to omit an input.
#' @param envir environment in which the function arguments are evaluated.
#'    By default, the calling environment of this function will be used.
#' @param ... other arguments.
#' @return A \code{tfm} object or a "um" if no input is significant at that level.
#' @export
varsel <- function (tfm, ...) { UseMethod("varsel") }

#' @rdname varsel
#' @export
varsel.tfm <- function(tfm, y = NULL, p.value = 0.10,
                       envir = parent.frame(), ...) {
  stopifnot(is.tfm(tfm))
  if (is.null(y)) y <- .output_tfm(tfm, y, envir)
  else tfm$noise$z <- deparse(substitute(y))
  
  p <- tryCatch(
    summary(tfm, y, p.values = TRUE, envir = envir),
    error = function(e) {
      stop("Failed to compute p-values: ", e$message)
    }
  )  
  if (is.null(p) || length(p) == 0) return(tfm)
  b <- param.tfm(tfm)
  nms <- names(b)
  names(p) <- nms
  xreg <- NULL
  if (!is.null(tfm$xreg)) {
    P <- (p[1:ncol(tfm$xreg)] <= p.value)
    if (all(P)) xreg <- tfm$xreg
    else if (any(P)) {
      xreg <- tfm$xreg[, P, drop = FALSE]
      if (!is.matrix(xreg)) {
        xreg <- as.matrix(xreg)
        colnames(xreg) <- colnames(tfm$xreg)[P]
      }
    } 
  } 
  
  tfi <- NULL
  if (tfm$k > 0) {
    tfi <- lapply(tfm$inputs, function(tf) {
      nmw0 <- as.character(tf$w0.expr)
      if (!nmw0 %in% names(p)) {
        warning(sprintf(
          "No p-value found for parameter '%s'. Keeping input.",
          nmw0
        ))
        return(tf)
      }
      if (p[nmw0] <= p.value) return(tf)
      else return(NULL)
    })  
    tfi[sapply(tfi, is.null)] <- NULL
  } 
  
  if (is.null(xreg) && is.null(tfi)) return(tfm$noise)
  tfm(y, xreg = xreg, inputs = tfi, noise = tfm$noise, new.name = FALSE, ...)
}

#' Coefficients of a Transfer Function Model
#'
#' Extracts the estimated coefficients from a fitted transfer function model 
#' of class \code{tfm}. This is a method for the generic \code{\link{coef}} 
#' function.
#'
#' @param object An object of class \code{tfm}.
#' @param ... Further arguments (currently unused).
#'
#' @return A named numeric vector with the estimated coefficients of the model, 
#' including regression coefficients, transfer function parameters, and noise 
#' model parameters.
#'
#' @seealso \code{\link{tfm}}
#'
#' @examples
#' \dontrun{
#' mdl <- tfm(y, xreg = X, noise = um())
#' coef(mdl)
#' }
#' @export
coef.tfm <- function(object, ...) {
  param.tfm(object)
}

#' Log-Likelihood of Transfer Function Model
#'
#' Computes the log-likelihood for a fitted transfer function model.
#'
#' @param object A fitted \code{tfm} object.
#' @param y Optional \code{ts} object for alternative output series.
#' @param method Character: "exact" estimates presample values; "cond" fixes
#'   presample values at zero.
#' @param envir Environment for evaluation. 
#' @param ... Additional arguments (currently unused).
#'
#' @return Numeric value of the log-likelihood.
#'
#' @seealso \code{\link{tfm}}, \code{\link{fit.tfm}}, \code{\link{AIC.tfm}}
#'
#' @export
logLik.tfm <-function(object, y = NULL, method = c("exact", "cond"), 
                      envir = parent.frame (), ...) {
  method <- match.arg(method)
  w <- noise.tfm(object, y, TRUE, envir=envir)
  if (!is.null(object$noise$mu))
    w <- w - object$noise$mu
  if (method == "exact") ll <- ellarmaC(w, object$noise$phi, object$noise$theta)
  else ll <- cllarmaC(w, object$noise$phi, object$noise$theta)
  return(ll)
}

#' AIC and BIC for Transfer Function Models
#'
#' Computes Akaike's Information Criterion (AIC) and Bayesian Information
#' Criterion (BIC) for transfer function models.
#'
#' @param object A fitted \code{tfm} object.
#' @param ... Additional \code{tfm} objects for model comparison.
#' @param k Numeric. Penalty per parameter. Default is 2 for AIC.
#'   Use \code{k = log(n)} for BIC where n is sample size.
#'
#' @return If one model: numeric value of AIC/BIC. If multiple models:
#'   data frame with columns df (degrees of freedom) and AIC for each model.
#'
#' @details
#' AIC = -2*logLik + k*npar, where npar is the number of parameters.
#' Lower values indicate better fit penalized for complexity.
#'
#' @seealso \code{\link{logLik.tfm}}, \code{\link{BIC.tfm}}
#'
#' @examples
#' \dontrun{
#' model1 <- tfm(output, inputs = tf1, noise = noise1)
#' model2 <- tfm(output, inputs = tf2, noise = noise2)
#' 
#' # Single model AIC
#' AIC(model1)
#' 
#' # Compare models
#' AIC(model1, model2)
#' 
#' # BIC
#' BIC(model1)
#' }
#'
#' @export
AIC.tfm <- function(object, ..., k = 2) {
  
  # Obtener lista de modelos
  models <- list(object, ...)
  
  # Si solo hay un modelo
  if (length(models) == 1) {
    ll <- logLik(object)
    npar <- length(param.tfm(object))
    aic <- -2 * ll + k * npar
    return(aic)
  }
  
  # Múltiples modelos: crear tabla comparativa
  model_names <- as.character(match.call()[-1])
  model_names <- model_names[model_names != "k"]
  
  results <- data.frame(
    df = integer(length(models)),
    AIC = numeric(length(models)),
    row.names = model_names
  )
  
  for (i in seq_along(models)) {
    if (!inherits(models[[i]], "tfm")) {
      stop(sprintf("Object %d is not a 'tfm' object", i))
    }
    ll <- logLik(models[[i]])
    npar <- length(param.tfm(models[[i]]))
    results$df[i] <- npar
    results$AIC[i] <- -2 * ll + k * npar
  }
  
  return(results)
}

#' @rdname AIC.tfm
#' @export
BIC.tfm <- function(object, ...) {
  
  # BIC usa k = log(n)
  y <- .output_tfm(object)
  w <- diffC(y, object$noise$nabla, object$noise$bc)
  n <- length(w)
  AIC.tfm(object, ..., k = log(n))
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
  stopifnot(is.tfm(object))
  rs <- residuals(object)
  stdres <- rs/sqrt(object$noise$sig2)
  n <- length(rs)
  if (gof.lag >= n || gof.lag < 1) {
    warning(sprintf(
      "gof.lag (%d) is too large for %d residuals. Using gof.lag = %d",
      gof.lag, n, floor(n / 2)
    ))
    gof.lag <- floor(n/4)
  }
  
  oldpar <- par(
    mfrow = c(3, 1),
    mar = c(3, 4, 2, 1) + 0.1,  # bottom, left, top, right
    oma = c(0, 0, 0, 0)
  )
  on.exit(par(oldpar))  
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
  
  
  invisible(list(
    std.residuals = stdres,
    acf = acf(rs, plot = FALSE, na.action = na.pass),
    ljung.box.pvalues = pval
  ))
  
}

#' @rdname sim
#' @param n Number of observations to simulate.
#' @param z0 Initial conditions for nonstationary series. Default is \code{NULL} (zero initial conditions).
#' @param n0 Number of initial observations to discard as burn-in. Default is \code{0}.
#' @param a Optional vector of innovations with length \code{n + n0}. If \code{NULL}, 
#'   innovations are drawn from \eqn{N(0, \sigma^2)}.
#' @param seed Random seed for reproducibility.
#' @param envir Environment for argument evaluation. Default is \code{parent.frame()}.
#' @export
sim.tfm <- function(mdl, n = 100, z0 = NULL, n0 = 0, a = NULL, seed = NULL, 
                    envir = parent.frame(), ...) {
  stopifnot(is.tfm(mdl))
  N <- sim.um(mdl$noise, n = n, z0 = z0, n0 = 0, a = a, seed = seed, 
           envir = envir,  ...)
  S <- signal.tfm(mdl, diff = FALSE, envir = envir)
  end <- end(N)
  freq <- frequency(N)
  
  if (mdl$noise$bc) unname(ts(exp(S)*N, end = end, frequency = freq))
  else unname(ts(S+N, end = end, frequency = freq))
}

#' @noRd
is.tfm <- function(tfm) {
  inherits(tfm, "tfm")
}

#' @noRd
param.tfm <- function(tfm) {
  unlist(tfm$param, use.names = TRUE)
}

# Update tfm object parameters from vector
# @param object tfm object
# @param param numeric vector of parameter values
# @return tfm object with updated parameters in all components
# @noRd
.update_tfm <- function(object, param) {
  object$param[] <- param
  object$inputs <- lapply(object$inputs, .update_tf, param = param)
  object$noise <- .update_um(object$noise, param)
  return(object)
}

# Remove differencing from tfm to predict transformed output
# Adjusts model by dividing out nabla from noise and applying it to inputs
# @param mdl tfm object
# @param nabla differencing operator to remove from model
# @return modified tfm for predicting differenced series
# @noRd
.diff_tfm <- function(mdl, nabla) {
  stopifnot(is.tfm(mdl))
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
  
  nabla <- polydivC(mdl$noise$nabla, nabla, FALSE, 1e-5)
  mdl$noise$nabla <- nabla
  mdl$noise$i <- list(as.lagpol(nabla))
  mdl
}

#' Extract output series from tfm object
#' @param object tfm object
#' @param y optional ts object (if NULL, extracts from object)
#' @param envir environment for evaluation
#' @return ts object containing the output series
#' @noRd
.output_tfm <- function(object, y = NULL, envir = parent.frame (), ...) {
  if (is.null(y)) {
    if (is.ts(object$output)) return(object$output)
    if (is.null(object$output)) stop("argment y required")
    else {
      y <- eval(parse(text = object$output), envir)
    }
  }  
  if (!is.numeric(y)) stop("output must be a numeric vector")
  as.ts(y)
} 







#' @rdname decomp
#' @param y an object of class \code{\link{ts}}.
#' @param envir environment in which the function arguments are evaluated.
#'    If NULL the calling environment of this function will be used.
#' @export
decomp.tfm <- function(mdl, y = NULL,
                      method = c("mixed", "forecast", "backcast"), 
                      envir = NULL, ...) {

  if (is.null (envir)) envir <- parent.frame ()
  y <- noise.tfm(mdl, y, FALSE, TRUE, envir = envir)
  uc <- decomp.um(mdl$noise, y, method)
  return(uc)

}
