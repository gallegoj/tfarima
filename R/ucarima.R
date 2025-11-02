## tfarima/R/ucarima.R
## Jose L Gallego (UC)

#' Unobserved components ARIMA models
#'
#' \code{ucarima} creates an S3 object that combines two or more ARIMA models
#' (objects of class \code{\link{um}}).
#'
#' @param z an object of class \code{ts} or NULL to specify theoretical models.
#' @param bc logical. If TRUE logs are taken.
#' @param ucm a list of \code{um} objects specifying the ARIMA models for
#'   components such as trend, seasonality, and irregular terms. Alternatively,
#'   can be a character string: \code{"ht"} for Harvey-Todd seasonal
#'   decomposition, or \code{NULL}/missing for Harvey-Durbin seasonal
#'   decomposition (default when only time series \code{z} is provided).
#' @param ar list of stationary AR lag polynomials (\code{lagpol} objects).
#' @param xreg matrix of explanatory variables.
#' @param fit logical. If TRUE, model is fitted.
#' @param envir environment.
#' @param ... additional arguments.
#'
#' @return An object of class \code{ucarima}.
#'
#' @references Harvey, A.C. (1989) Forecasting, Structural Time Series Models
#' and the Kalman Filter. Cambridge University Press, Cambridge.
#'
#' @examples
#' trend <- um(i = "(1 - B)", sig2 = c(s2t = 1))
#' seas <- um(i = "(1+B)", sig2 = c(s2s = 0.05))
#' irreg <- um(sig2 = c(s2i = 0.75))
#' uca1 <- ucarima(ucm = list(trend = trend, seas = seas, irreg = irreg))
#' uca1
#'
#' # Trigonometric seasonality
#' uca2 <- ucarima(AirPassengers, bc = TRUE)
#' uca2
#' 
#' # Dummy seasonality
#' uca3 <- ucarima(AirPassengers, bc = TRUE, ucm = "ht")
#' uca3
#'
#' @export
ucarima <- function(z = NULL, bc = FALSE, ucm = NULL, ar = NULL, 
                    xreg = NULL, fit = TRUE, envir = parent.frame(), ...) {
  call <- match.call()
  z.name <- if (is.numeric(z)) deparse(substitute(z)) else NULL
  
  # Construct default UCM if not provided
  if (is.null(ucm) || is.character(ucm)) {
    if (is.null(z)) stop("argument 'z' is missing when 'ucm' is not specified") 
    stopifnot("'z' must be a time series object" = is.ts(z))
    stopifnot("time series must have frequency > 1" = frequency(z) > 1)
    s <- frequency(z)
    
    # Trend: IMA(2, 1) with positive MA parameter
    lp <- lagpol(param = c(theta = 0.8), coef = "1/(1+exp(-theta))")
    trend <- um(i = 2, ma = lp, sig2 = c(s2t = 1))
    # Irregular: white noise
    irreg <- um(sig2 = c(s2i = 0.75))
    # Seasonal
    if (s > 1) {
      if (!is.null(ucm) && ucm == "ht") {
        seas <- um(i = as.character(s), sig2 = c(s2s = 0.25))
      } else {
        if (s %% 2 == 0) seas <- um(i = "1/2", sig2 = 0.5)
        else seas <- um(sig2 = 0)
        if (s > 2) {
          k <- as.integer((s-1)/2)
          for (j in 1:k) {
            i <- paste0(j, "/", s)
            f <- 2 * pi * j / s 
            s2 <- (1 + sin(f))
            ma <- lagpol(coef = cos(f)/s2)
            seas <- seas + um(i = i, ma = ma, sig2 = s2)
          }
          ma <- lagpol(coef = -seas$theta[-1])
        } else ma <- NULL
        seas <- um(i = seas$i, ma = ma, sig2 = c(s2s = 0.25))
      }
      ucm <- list(trend = trend, seas = seas, irreg = irreg)
    } else ucm <- list(trend = trend, irreg = irreg)
  }

  # Validate UCM structure
  stopifnot("all UCM components must be 'um' objects" = all(sapply(ucm, is.um)))
  stopifnot("UCARIMA must contain at least 2 components" = length(ucm) > 1)    
  # Check for common factors between components
  k <- length(ucm)
  for (h in seq_len(k-1)) {
    if (length(ucm[[h]]$phi) > 1) {
      for (j in (h+1):k) {
        if (length(ucm[[j]]$phi) > 1) {
          p <- common.factors(ucm[[h]]$phi, ucm[[j]]$phi)
          if (length(p) > 0)
            stop("common AR factors are not allowed between components")
          p <- common.factors(ucm[[h]]$nabla, ucm[[j]]$nabla)
          if (length(p) > 0)
            stop("common I factors are not allowed between components")
        }
      }
    }
  }

  param <- lapply(ucm, function(x) unlist(x$param))
  param <- unname(param)
  if (!is.null(ar)) {
    ar <- lagpol0(ar, "ar")
    param1 <- lapply(ar, function(x) unlist(x$param))
    param1 <- unname(param1)
    param <- c(param1, param)
  }
  param <- unlist(param)
  param <- param[!duplicated(names(param))]
  param <- as.list(param)
   
  if (is.null(names(ucm)))
    names(ucm) <- paste0("uc", 1:k)

  if (!is.null(xreg) && !is.null(z))
    a <- solve(t(xreg) %*% xreg, t(xreg) %*% z)
  else
    a <- NULL
  
  mdl <- structure(list(ucm = ucm, k = k, ar = ar, xreg = xreg, a = a, z = z, 
                        z.name = z.name, bc = bc, call = call, nabla = NULL, 
                        p = 0, d = 0, q = 0, param = param, is.adm = TRUE),
                   class = "ucarima")
  mdl$nabla <- nabla(mdl, lp = FALSE)
  mdl$d <- length(mdl$nabla) - 1
  p <- sapply(1:k, function(x) mdl$ucm[[x]]$p)
  d <- sapply(1:k, function(x) mdl$ucm[[x]]$d)
  mdl$q <- max(sapply(1:k, function(x) sum(p[-x]) + sum(d[-x]) + mdl$ucm[[x]]$q))
  mdl$p <- sum(p)
  if (!is.null(ar)) {
    p <- sapply(ar, function(x) length(x$Pol) - 1)
    mdl$p <- mdl$p + sum(p)
  }
  mdl <- .update_ucarima(mdl, c(a, unlist(param, use.names = TRUE)))
  if (mdl$is.adm) {
    if (!is.null(z) && fit) {
      mdl <- tryCatch({
        mdl <- fit.ucarima(mdl, envir = envir, ...)
        mdl
      }, error = function(e) {
        warning("Model fitting failed: ", conditionMessage(e))
        mdl$is.adm <- FALSE
        mdl
      })      
    } else mdl$b <- c(a, unlist(param, use.names = TRUE))
  } else {
    warning("non-admissible model")
  }
  return(mdl)
}

#' Estimation of UCARIMA models
#'
#' Estimates the parameters of a UCARIMA model using maximum likelihood.
#'
#' @param mdl An object of class \code{\link{ucarima}}.
#' @param z Optional. A time series object (\code{ts}) used for estimation if
#'   not already included in \code{mdl}.
#' @param method Character string specifying the optimization method passed to 
#'   \code{\link[stats]{optim}}. Default is \code{"BFGS"}.
#' @param show.iter Logical. If \code{TRUE}, displays parameter estimates during
#'   each iteration of the optimization process. Default is \code{FALSE}.
#' @param envir Environment where the estimation is evaluated. If \code{NULL}
#'   (default), uses the parent environment.
#' @param ... Additional arguments passed to \code{\link[stats]{optim}}.
#'
#' @return An updated \code{ucarima} object with estimated parameters,
#'   log-likelihood, and convergence information.
#'
#' @examples
#' # Define a local level model with trend and irregular components
#' trend <- um(i = 1, sig2 = c(s2t = 0.5))
#' irreg <- um(sig2 = c(s2i = 1))
#'
#' # Create and estimate the UCARIMA model
#' uca <- ucarima(ucm = list(trend = trend, irreg = irreg))
#' uca_fitted <- fit(uca, Nile)
#'
#' @export
fit.ucarima <- function(mdl, z = NULL, method = "BFGS", show.iter = FALSE,
                        envir = NULL, ...) {
  
  stopifnot("mdl must be a ucarima object" = inherits(mdl, "ucarima"))
  if (!is.environment(envir)) envir <- parent.frame()

  # Log-likelihood function for UCARIMA model
  .lluca <- function(par, type = 0, ...) {
    if (kX + k1 > 0) mdl <<- .update_ucarima(mdl, par[1:(kX+k1)])
    if (mdl$is.adm) {
      # Construct MA polynomial matrix A      
      if (maxp > 0) {
        A <- sapply(seq_len(mdl$k), function(x) {
          ma <- phi.ucarima(mdl, x, FALSE, FALSE)
          ma <- polymultC(ma, D[[x]])
          ma <- polymultC(ma, mdl$ucm[[x]]$theta)
          c(ma, rep(0, mdl$q - length(ma) + 1))
        })
      } else {
        A <- sapply(seq_len(mdl$k), function(x) {
          ma <- as.numeric(polymultC(mdl$ucm[[x]]$theta, D[[x]]))
          c(ma, rep(0, mdl$q - length(ma) + 1))
        })
      }
      A <- t(A)
      # AR polynomial
      if (mdl$p > 0) phi <- phi.ucarima(mdl, lp = FALSE)
      else phi <- 1
      if (kX + k1 > 0) sigmas[!s2fix] <<- exp(par[-(1:(kX + k1))])^2
      else sigmas[!s2fix] <<- exp(par)^2
      S <- diag(colSums(E * sigmas))
      if (!is.null(mdl$xreg)) {
        if (type == 0)
          ll <- llucaC(w - X %*% par[1:kX], phi, A, S, s2, type)
        else {
          w <<- w - X %*% par[1:kX] # w is changed by llucaC function
          ll <- llucaC(w, phi, A, S, s2, type)
        }
      } else {
        ll <- llucaC(w, phi, A, S, s2, type)
      }
      if (show.iter)
        print(c(ll = ll, s2, par))
    } else ll <- ll0
    mdl$is.adm <<- TRUE
    return(-ll)
  }
  
  # Residuals function
  .resuca <- function(par, type = 1, ...) {
    w1 <- w*1
    ll <- .lluca(par, type)
    e <- w*1
    w <<- w1
    return(e)
  }
  
  if (is.null(z)) {
    if (!is.ts(mdl$z)) {
      if (is.null(mdl$z.name)) stop("argment z required")
      else mdl$z <- eval(parse(text = mdl$z.name), envir)
    }
  } else {
    mdl$z <- z
    mdl$z.name <- deparse(substitute(z))
  }
  
  if (mdl$bc && min(mdl$z) <= 0) mdl$bc <- FALSE
  w <- diffC(mdl$z, mdl$nabla, mdl$bc)
  # Differencing operators for each component
  D <- lapply(1:mdl$k, function(x) {
    nabla.ucarima(mdl, x, FALSE)
  })

  maxp <- max(sapply(1:mdl$k, function(x) mdl$ucm[[x]]$p))

  # Extract variance parameters  
  sigmas <- sapply(1:mdl$k, function(x) {
    mdl$ucm[[x]]$sig2
  })
  s2nms <- names(sigmas)
  if (all(s2nms == "sig2") || is.null(s2nms)) {
    s2nms <- paste0("s2.", 1:mdl$k)
    names(sigmas) <- s2nms  
    for (i in 1:mdl$k)
      names(mdl$ucm[[i]]$sig2) <- paste0("s2.", i)
  }
  nms <- s2nms[!duplicated(s2nms)]
  E <- t(sapply(nms, function(x) s2nms == x))
  E <- t(t(E)*sigmas)
  E <- apply(E, 1, function(x) {
    if (max(x) != 0) x/max(x) else x
  })
  E <- t(E)
   
  sigmas <- sigmas[!duplicated(s2nms)]
  s2con <- nms[sigmas == max(sigmas)][1]
  s2con <- (nms == s2con)
  s2 <- sigmas[s2con]
  s2fix <- (sigmas == 0) | s2con
  sigmas <- sigmas/s2

  if (!is.null(mdl$xreg)) {
    kX <- ncol(mdl$xreg)
    X <- sapply(1:kX, function(i) diffC(mdl$xreg[, i], mdl$nabla, FALSE))
    X <- X[1:length(w), ]
    if (kX == 1) X <- as.matrix(X)
    mdl$a <- tryCatch({
      as.numeric(solve(t(X) %*% X, t(X) %*% w))
    }, error = function(e) {
      stop("Singular design matrix in regression")
    })    
    if (length(colnames(mdl$xreg)) == kX)
      names(mdl$a) <- colnames(mdl$xreg)
    else
      names(mdl$a) <- paste0("a", 1:kX)
    par <- mdl$a
  } else {
    kX <- 0
    par <- NULL
  }
  
  k1 <- length(mdl$param)
  if(k1 > 0) par <- c(par, unlist(mdl$param))
  
  k2 <- length(sigmas[!s2fix])
  if (k2 > 0) par <- c(par, log(sqrt(sigmas[!s2fix])))
  else if(kX + k1 == 0) stop("All parameters are fixed")
  
  ll0 <- -sqrt(.Machine$double.xmax)
  opt <- tryCatch({
    optim(par, .lluca, hessian = TRUE, method = method, ...)
  }, error = function(e) {
    stop("Optimization failed: ", conditionMessage(e))
  })  
  # Check convergence
  if (opt$convergence != 0) {
    warning("Optimization did not converge (code: ", opt$convergence, ")")
  }  
  
  # Update model with estimated parameters  
  if (kX + k1 > 0) mdl <- .update_ucarima(mdl, opt$par[1:(kX + k1)])
  # Update variance components  
  for (i in 1:mdl$k) {
    mdl$ucm[[i]]$sig2 <- E[nms == names(mdl$ucm[[i]]$sig2), i]*
      sigmas[nms == names(mdl$ucm[[i]]$sig2)]*s2
  }
  
  mdl$um <- as.um(mdl)
  if (!is.null(mdl$a))
    mdl$res <- residuals(mdl$um, mdl$z - as.matrix(mdl$xreg[1:length(mdl$z),]) %*% mdl$a)
  else
    mdl$res <- residuals(mdl$um, mdl$z)
  
  # Compute auxiliary statistics  
  mdl$aux <- list()
  mdl$aux$n <- length(w)
  mdl$aux$logLik <- -opt$value
  mdl$aux$sig2 <- mdl$um$sig2
  mdl$aux$aic <- 2.0*(opt$value + (k1+k2))/mdl$aux$n
  mdl$aux$sigmas <- cbind(value = sigmas*s2, ratio = sigmas)
  # Compute gradient and variance-covariance matrix  
  res <- .resuca(opt$par, ...)
  J <- tryCatch({
    numDeriv::jacobian(.resuca, opt$par, ...)
  }, error = function(e) {
    warning("Failed to compute Jacobian: ", conditionMessage(e))
    matrix(NA, nrow = length(res), ncol = length(opt$par))
  })  
  
  mdl$aux$g <- t(J) %*% res

  # Variance-covariance matrix with error handling
  mdl$aux$varb <- tryCatch({
    solve(t(J) %*% J) * (sum(res^2) / length(res))
  }, error = function(e) {
    warning("Failed to compute variance-covariance matrix: ", conditionMessage(e))
    matrix(NA, nrow = length(opt$par), ncol = length(opt$par))
  })    

  # Parameter summary with standard errors
  if (!anyNA(mdl$aux$varb)) {
    se <- sqrt(diag(mdl$aux$varb))
  } else {
    se <- rep(NA, length(opt$par))
  }    
  
  mdl$aux$par <- cbind(par = opt$par, se = se)
  mdl$optim <- opt
  
  return(mdl)
  
}

#' @rdname diagchk
#' @description
#' For objects of class \code{ucarima}, this method calls \code{diagchk.um}
#' internally to perform diagnostic checking.
#' @param mdl an object of class \code{ucarima}
#' @export
diagchk.ucarima <- function(mdl, ...) {
  diagchk.um(as.um(mdl), ...)
}

#' @rdname decomp
#' @export
decomp.ucarima <- function(mdl, ...) {
  stopifnot("mdl must be a ucarima object" = inherits(mdl, "ucarima"))
  stopifnot("model must contain time series data" = !is.null(mdl$z))
  
  if (is.null(mdl$um)) 
    mdl$um <- as.um(mdl)
  # Incorporate AR polynomials into components if present  
  if (!is.null(mdl$ar)) {
    mdl$ucm <- lapply(mdl$ucm, function(x) {
      modify.um(x, ar = mdl$ar, fit = FALSE)
    })
  }
  
  # Extract components using Wiener-Kolmogorov filter  
  uc1 <- sapply(1:mdl$k, function(x) {
    v <- wkfilter.um(mdl$um, mdl$ucm[[x]], z = mdl$z)    
    if (mdl$bc) exp(v) else v
  })
  uc1 <- cbind(mdl$z, uc1)
  colnames(uc1) <- c(mdl$z.name, names(mdl$ucm))
  rownames(uc1) <- NULL
  uc1 <- ts(uc1, end = end(mdl$z), frequency = frequency(mdl$z))
  
  return(uc1)
}

#' Residuals of fitted UCARIMA models
#'
#' \code{residuals.ucarima} generates the residuals of a fitted \code{\link{ucarima}}
#' object.
#'
#' @param object an object of class \code{\link{ucarima}}.
#' @param method character. Either "exact" or "conditional" residuals.
#' @param ... additional arguments.
#' @return A \code{ts} object containing the residuals of the SSM.
#' @details These residuals are calculated by first converting the
#'   \code{\link{ucarima}} object to a \code{\link{um}} object and then using the
#'   \code{\link{residuals.um}} function.
#' @export
residuals.ucarima <- function(object, method = c("exact", "cond"), ...) {
  stopifnot(inherits(object, "ucarima"))
  stopifnot(!is.null(object$z))
  method <- match.arg(method)
  um1 <- as.um(object)
  um1$bc <- FALSE
  if (object$bc) z <- log(object$z)
  else z <- object$z 
  N <- length(z)
  if (!is.null(object$xreg)) 
    z <- z - object$xreg[1:N, ] %*% object$a
  res <- residuals.um(um1, z, method = method, ...)
  res <- ts(res, end = end(object$z), frequency = frequency(object$z))
  return(res)
}





#' @rdname outliers
#' @export
outliers.ucarima <- function(mdl, types = c("AO", "LS", "TC"), dates = NULL, c = 3, 
                             calendar = FALSE, easter = FALSE, 
                             resid = c("exact", "cond"), n.ahead = 0, p.value = 1, 
                             ...) {
  um1 <- as.um(mdl)
  if (!is.null(mdl$xreg)) {
    n <- length(mdl$z)
    if (mdl$bc) z <- log(mdl$z) - mdl$xreg[1:n, ] %*% mdl$a
    else z <- mdl$z - mdl$xreg[1:n, ] %*% mdl$a
  } else {
    if (mdl$bc) z <- log(mdl$z)
    else z <- mdl$z
  }
  um1$bc <- FALSE
  um1$param <- NULL
  tfm1 <- tfm(z, noise = um1, fit = FALSE, new.name = FALSE, envir = NULL)
  tfm1 <- outliers.tfm(tfm1, z, types, dates, c, calendar, easter, resid, n.ahead, 
                       p.value, TRUE, NULL, ...)
  if (is.null(tfm1$xreg))
    return(mdl)
  if (!is.null(mdl$xreg))
    mdl$xreg <- cbind(mdl$xreg, tfm1$xreg)
  else
    mdl$xreg <- tfm1$xreg
  fit(mdl)
}


#' @rdname as.ssm
#' @export
as.ssm.ucarima <- function(object, ...) {
  b <- c()
  s2 <- 0
  for (i in 1:object$k) {
    if (object$ucm[[i]]$p + object$ucm[[i]]$d + object$ucm[[i]]$q >  0) {
      sf <- as.ssm(object$ucm[[i]])
      if (is.null(b)) {
        b <- sf$b
        C <- sf$C
        S <- as.matrix(sf$S[-1, -1])
      } else {
        b <- c(b, sf$b)
        A <- sf$C
        C <- rbind(cbind(C, matrix(0, nrow = nrow(C), ncol = ncol(A))),
              cbind(matrix(0, nrow = nrow(A), ncol = ncol(C)), A))      
        A <- as.matrix(sf$S[-1, -1])
        S <- rbind(cbind(S, matrix(0, nrow = nrow(S), ncol = ncol(A))),
                   cbind(matrix(0, nrow = nrow(A), ncol = ncol(S)), A))      
      }
      s2 <- sf$S[1, 1]
    } else s2 <- s2 + object$ucm[[i]]$sig2    
  }
  S <- rbind(matrix(0, nrow = 1, ncol = ncol(S) + 1),
             cbind(matrix(0, nrow = nrow(S), ncol = 1), S))
  S[1, 1] <- s2
  sf <- ssm(b = b, C = C, S = S)
  sf$y <- object$z
  sf$y.name <- object$z.name
  return(sf)
}  


#' @rdname phi
#' @param i integer. Omit this component model.
#' @param ar logical. If TRUE, the common AR polynomial is included.
#' @param lp logical indicating the type of return: \code{\link{lagpol}} object
#'   or numeric vector.
#' @export
phi.ucarima <- function(x, i = NULL, ar = TRUE, lp = TRUE, ...) {
  stopifnot(inherits(x, "ucarima"))
  h <- 1:x$k
  if (!is.null(i)) h <- h[-abs(i)]
  pol <- 1
  for (j in h)
    pol <- polymultC(pol, x$ucm[[j]]$phi)
  if ( (ar == TRUE) && !is.null(x$ar))
    pol <- polymultC(pol, polyexpand(x$ar))
  if (lp)
    return(as.lagpol(pol))
  else
    return(pol)
}

#' @rdname nabla
#' @param i integer. Omit this component model.
#' @param lp logical indicating the type of return: \code{\link{lagpol}} object
#'   or numeric vector.
#'
#' @export
nabla.ucarima <- function(x, i = NULL, lp = TRUE, ...) {
  stopifnot(inherits(x, "ucarima"))
  j <- 1:x$k
  if (!is.null(i) && i != 0) j <- j[-abs(i)]
  pol <- polyexpand(lapply(x$ucm[j], function(x) x$nabla))
  if (lp) return(as.lagpol(pol))
  else return(pol)
}

  


#' Print ucarima models
#' @rdname print
#' @export
print.ucarima <- function(x, ...) {
  if (is.null(x$optim)) {
    nms <- names(x$ucm)
    for (i in 1:x$k) {
      if (is.um(x$ucm[[i]])) {
        cat(nms[i], "\n") 
        equation(x$ucm[[i]], FALSE)
      }
    }
  } else {
    print(x$aux$par)
    cat("\n");
    print(x$aux$sigmas)
    cat("\nlog likelihood: ", x$aux$logLik)
    cat("\nResidual standard error: ", sqrt(x$aux$sig2))
    cat("\naic:", x$aux$aic)
  }
  return(invisible(NULL))
}

#' Equation of ucarima model
#' @rdname equation
#' @export
equation.ucarima <- function(x, ...) {
  nms <- names(x$ucm)
  for (i in 1:x$k) {
    if (is.um(x$ucm[[i]]))
      cat(nms[i], "\n", equation(x$ucm[[i]], FALSE, ...), "\n")
  }
  return(invisible(NULL))
}


#' @export
coef.ucarima <- function(object, ...) {
  (object$param)
}

#' @rdname autocov
#' @param ma logical; if true, autocovariances are computed for the MA model. By
#'   default, ma = FALSE and autocovariances are computed for the ARMA model.
#' @export
autocov.ucarima <- function(mdl, ma = FALSE, ...) {
  um1 <- as.um(mdl)
  if (ma) {
    um1$ar <- NULL
    um1$phi <- 1
    um1$p <- 0
  }
  autocov.um(um1, ...)
}

#' @rdname autocorr
#' @export
autocorr.ucarima <- function(x, ...) {
  um1 <- as.um(x)
  autocorr.um(um1, ...)
}



# Internal functions for ucarima objects

#' @keywords internal
#' @noRd
.ucarima2um <- function(x, ...) {
  mdl <- x$ucm[[1]]
  for (i in 2:x$k) {
    mdl <- mdl + x$ucm[[i]]
  }
  if (!is.null(x$ar))
    mdl <- modify.um(mdl, ar = x$ar, fit = FALSE)
  mdl$z <- x$z.name
  mdl$bc <- x$bc
  return(mdl)
}

#' Update an object \code{ucarima}
#'
#' \code{.update_ucarima} updates an object of class \code{ucarima} with a new
#' vector of parameters
#'
#' @param x an object of class \code{ucarima}.
#' @param b a numeric vector.
#'
#' @return An object of class \code{um}.
#' @keywords internal
#' @noRd
.update_ucarima <- function(x, b, ...) {
  if (!is.null(x$a)) x$a <- b[1:length(x$a)]
  if (is.null(names(x$param))) return(x)
  x$param[] <- b[names(x$param)] 
  x$is.adm <- TRUE
  if (!is.null(x$ar)) {
    x$ar <- lapply(x$ar, function(pol) {
      pol <- .update.lagpol(pol, b)
      ok <- admreg.lagpol(pol)
      if (!ok) x$is.adm <<- FALSE
      pol
    })
  }
  
  for (i in 1:x$k) {
    x$ucm[[i]] <- .update_um(x$ucm[[i]], b)
    if (!x$ucm[[i]]$is.adm) x$is.adm <- FALSE
  }  
  
  return(x)
  
}

