## tfarima/R/stsm.R
## Jose L Gallego (UC)

#' Structural time series models
#'
#' \code{stsm} creates an S3 object representing a time-invariant structural
#' time series model:
#'
#' y(t) = b'x(t) + u(t) (observation equation), 
#' x(t+j) = Cx(t+j-1) + v(t) (state equation),
#' j = 0 for contemporaneous form or 1 for future form.
#'
#' @param y an object of class \code{ts}.
#' @param b vector of constants.
#' @param C matrix of constants.
#' @param S covariance matrix of the error vector (u_t, v_t).
#' @param xreg matrix of regressors.
#' @param bc logical. If TRUE logs are taken.
#' @param cform logical. If TRUE station equation is given in contemporaneous
#'   form, otherwise it is written in future form.
#'
#' @return An object of class \code{stsm}.
#'
#' @references
#'
#' Durbin, J. and Koopman, S.J. (2012) Time Series Analysis by State Space
#' Methods, 2nd ed., Oxford University Press, Oxford.
#'
#' Harvey, A.C. (1989) Forecasting, Structural Time Series Models and the Kalman
#' Filter. Cambridge University Press, Cambridge.  
#'
#' @examples
#' # Local level model
#' b <- 1
#' C <- as.matrix(1)
#' stsm1 <- stsm(Nile, b, C, S = diag(c(irr = 1, lvl = 0.5)) )
#' stsm1
#' @export
#' 
stsm <- function(y, b, C, S, xreg = NULL, bc = FALSE, cform = TRUE) {
  if (!is.matrix(C)) {
    if (length(C) == 1) C <- as.matrix(C)
    else stop("C argument must be a matrix")
  }
  stopifnot(nrow(C) == ncol(C))
  if (!is.matrix(S)) {
    stop("S argument must be a matrix")
  }
  stopifnot(nrow(S) == ncol(S))
  if (is.matrix(b)) 
    b <- as.vector(b)
  stopifnot(length(b) == nrow(C) || (length(b)+1) == nrow(S))
  
  if (missing(y)) y <- NULL
  if (!is.null(y)) {
    y.name <- deparse(substitute(y))
    if (!is.null(xreg)) {
      if (is.matrix(xreg))
        stopifnot(nrow(xreg) == length(y))
      else if (is.vector(xreg)||is.ts(xreg)) {
        stopifnot(length(xreg) == length(y))
        xreg <- matrix(xreg, ncol = 1)
      } else stop("wrong xreg argument")
    }
  } else y.name <- NULL
  if (!is.null(xreg)&&!is.null(y))
    a <- solve(t(xreg) %*% xreg, t(xreg) %*% y)
  else
    a <- NULL
  
  mdl <- list(y = y, y.name = y.name, bc = bc, b = b, C = C, S = S,
              xreg = xreg, a = a, cform = cform)
  class(mdl) <- "stsm"
  mdl
}

#' Alternative form for STS model
#'
#' \code{altform} converts a STS model given in contemporaneous form into future
#' form, and vice versa.
#'
#' @param mdl an object of class \code{stsm}.
#' @param ... other arguments.
#' @return An object of class \code{stsm}.
#'
#' @examples
#' # Local level model
#' b <- 1
#' C <- as.matrix(1)
#' stsm1 <- stsm(Nile, b, C, S = diag(c(irr = 1, lvl = 0.5)) )
#' stsm1
#' stsm2 <- altform(stsm1)
#' 
#' @export
altform <- function (mdl, ...) { UseMethod("altform") }

#' @rdname altform
#' @export
altform.stsm <- function(mdl, ...) {
  if (mdl$cform) {
    s11 <- mdl$S[1, 1] + sum((mdl$b %*% mdl$S[-1, -1])*mdl$b) + 2*sum(mdl$b*mdl$S[1, -1])
    mdl$S[-1, 1] <- mdl$S[1, -1] <- mdl$S[1, -1] + (mdl$b %*% mdl$S[-1, -1])
    mdl$S[1, 1] <- s11
    mdl$b <- as.numeric(mdl$b %*% mdl$C)
    mdl$cform <- FALSE
  } else {
    b <- mdl$b %*% solve(mdl$C)
    s11 <- mdl$S[1, 1] + sum((b %*% mdl$S[-1,-1])*b) - 2*sum(b*mdl$S[1, -1])
    mdl$S[-1, 1] <- mdl$S[1, -1] <- mdl$S[1, -1] - (b %*% mdl$S[-1, -1])
    mdl$S[1, 1] <- s11
    mdl$b <- as.numeric(b)
    
    mdl$cform <- TRUE
  }
  mdl$S[abs(mdl$S) < .Machine$double.eps] <- 0
  mdl
}

#' Non stationary AR polynomial of the reduced form
#' 
#' \code{nabla} returns the non stationary AR polynomialof the reduced form of
#'  an object of the \code{stsm} class.
#'
#' @param mdl an object of class \code{stsm}.
#' @param tol tolerance to check if a root is close to one. 
#' 
#' @return 
#' A numeric vector \code{c(1, a1, ..., ad)}
#' 
#' @examples
#' stsm1 <- stsm(b = 1, C = 1, S = diag(c(0.8, 0.04)))
#' nabla(stsm1)
#' @export
nabla.stsm <- function(mdl, tol = 1e-4) {
  lf <- .LeverrierFaddeev(mdl$b, mdl$C)
  r <- polyroot(lf$p) 
  nabla <- any(abs(abs(r)-1) <= tol)
  if (nabla) {
    nabla <- roots2lagpol(r[abs(abs(r)-1) <= tol])$Pol
  } else {
    nabla <- 1
  }
  return(as.lagpol(nabla, coef.name = "i"))
}

#' @rdname autocov
#' @param arma logical. If TRUE, the autocovariances for the stationary ARMA
#'   model of the reduced form are computed. Otherwise, the autocovariances are
#'   only computed for the MA part.
#' @param varphi logical. If TRUE, the varphi polynomial of the reduced form is
#'   also returned.
#' @param tol tolerance to check if a root is close to one.    
#' @examples
#' # Local level model
#' stsm1 <- stsm(b = 1, C = 1, S = diag(c(irr = 0.8, lvl = 0.04))
#' autocov(stsm1)
#'
#' @export
autocov.stsm <- function(mdl, lag.max = NULL, arma = TRUE, varphi = FALSE,
                         tol = 1e-4, ...) {
  lf <- .LeverrierFaddeev(mdl$b, mdl$C)
  d <- lf$p
  k <- length(d)
  if (is.null(lag.max)) lag.max <- k
  else lag.max <- abs(lag.max) + 1
  A <- lf$A
  if (mdl$cform)
    A <- cbind(t(A), rep(0, k-1))
  else
    A <- cbind(rep(0, k-1), t(A))
  A <- rbind(d, A)
  r <- polyroot(d)
  ar <- any(abs(abs(r)-1) > tol)
  if (ar && arma) {
    g <- rep(0, k)
    phi <- roots2lagpol(r[abs(abs(r)-1) > tol])$Pol
    psi <- A
    p <- length(phi)
    if (k > 1) {
      for (i in 2:k) {
        for (j in 2:min(i,p))
          psi[, i] <- psi[, i] - phi[j]*psi[, i-j+1]
      }
    }
    for (i in 1:k) {
      sum <- 0
      for (j in i:k) 
        sum <- sum + sum((A[, j] %*% mdl$S) * psi[, j-i+1])
      g[i] <- g[i] + sum
    }
    B <- toeplitz(phi)
    if (p > 1) {
      B[upper.tri(B)] <- 0
      B1 <- sapply(1:p, function(x) c(0, phi[-(1:x)], rep(0,x-1)))
      B <- B + t(B1)
    }
    x <- solve(B, g[1:p])
    if (lag.max > p && lag.max > k+1) {
      x <- c(x, rep(0, lag.max-p))
      g <- c(g, rep(0, lag.max-k-1))
      for(i in (p+1):lag.max) {
        x[i] <- g[i]
        for (j in 2:p)
          x[i] <- x[i] - phi[j]*x[i-j+1]
      }
    }
    g <- x
  } else {
    g <- rep(0, k)
    for (i in 1:k) {
      sum <- 0
      for (j in i:k) 
        sum <- sum + sum((A[, j] %*% mdl$S) * A[, j-i+1])
      g[i] <- sum
    }
    if (length(g) < lag.max) 
      g <- c(g, rep(0, lag.max - length(g)))
    else
      g <- g[1:lag.max]
  } 
  if (varphi)
    return(list(g = g, varphi = lf$p))
  else
    return(g)
}

#' Estimation of a STS model
#'
#' \code{fit} is used to estimate a stsm model.
#'
#' @param mdl an object of class \code{\link{stsm}}.
#' @param updmod  function to update the parameters of the model.
#' @param par argument passed to the \code{updmod} function.
#' @param fixed vector of logical values indicating which parameters are fixed
#'   (TRUE) or estimated (FALSE).
#' @param show.iter logical value to show or hide the estimates at the different
#'   iterations.
#' @param tol tolerance to check if a root is close to one.
#' @param ... other arguments.
#'
#' @return An object of class "stsm" with estimated parameters.
#'
#' @examples
#' # Local level model
#' b <- 1
#' C <- as.matrix(1)
#' stsm1 <- stsm(Nile, b, C, s2v = c(lvl = 0.5), s2u = c(irr = 1), fit = FALSE)
#' stsm1 <- fit(stsm1, method = "L-BFGS-B")
#' @export
fit.stsm <- function(mdl, updmdl, par, fixed = NULL, show.iter = FALSE, 
                     tol = 1e-4, ...) {
  .loglikrf <- function(par, ...) {
    if (bfixed) {
      parfix[!fixed] <<- par[1:k]
      lst <- updmdl(parfix, ...)
    } else lst <- updmdl(par[1:k], ...)
    if (bb) mdl$b <<- lst$b
    if (bC) mdl$C <<- lst$C
    if (bS) mdl$S <<- lst$S
    if (!is.null(mdl$xreg)) {
      ll <- llrfC(w-X%*%par[-(1:k)], nabla, mdl$b, mdl$C, mdl$S, s2, mdl$cform)
    } else {
      ll <- llrfC(w, nabla, mdl$b, mdl$C, mdl$S, s2, mdl$cform)
    }
    if (show.iter)
      print(c(ll = ll, s2 = s2, par))
    -ll
  }
  
  .res <- function(par, ...) {
    if (bfixed) {
      parfix[!fixed] <<- par[1:k]
      lst <- updmdl(parfix, ...)
    } else lst <- updmdl(par[1:k], ...)
    if (bb) mdl$b <<- lst$b
    if (bC) mdl$C <<- lst$C
    if (bS) mdl$S <<- lst$S
    if (!is.null(mdl$xreg)) {
      e <- resrfC(w-X%*%par[-(1:k)], nabla, mdl$b, mdl$C, mdl$S, s2, mdl$cform)
    } else {
      e <- resrfC(w, nabla, mdl$b, mdl$C, mdl$S, s2, mdl$cform)
    }
    e
  }
  
  .sigmas <- function(par, ...) {
    j <- length(par)
    for (i in 1:j) {
      p1 <- par[i]
      par[i] <- NA
      lst <- updmdl(par[1:j], ...)
      if (any(is.na(lst$S)))
        par[i] <- diag(mdl$S)[is.na(diag(lst$S))][1]
      else par[i] <- p1
    }
    par
  }
  
  if(is.null(mdl$y)) stop("missing time series")
  # Check update function
  lst <- updmdl(par, ...)
  stopifnot(is.list(lst))
  nm <- names(lst)
  if (!all(nm %in% c("b", "C", "S") )) {
    stop(paste0("invalid names: ", nm[!(nm %in% c("b", "C", "S"))]))
  }
  bb <- any(nm == "b")
  bC <- any(nm == "C")
  bS <- any(nm == "S")
  if (bb) {
    stopifnot(length(lst$b) == length(mdl$b))
    mdl$b <- lst$b
  }
  if (bC) {
    stopifnot(all(dim(lst$C) == dim(mdl$C)))
    mdl$C <- lst$C
  }
  if (bS) {
    stopifnot(all(dim(lst$S) == dim(mdl$S)))
    mdl$S <- lst$S
  }

  lf <- .LeverrierFaddeev(mdl$b, mdl$C)
  if (bC) {
    r <- polyroot(lf$p) 
    nabla <- any(abs(abs(r)-1) <= tol)
    if (nabla) {
      nabla <- roots2lagpol(r[abs(abs(r)-1) <= tol])$Pol
      phi <- polydivC(lf$p, nabla, FALSE)
    } else {
      nabla <- 1
      phi <- lf$p
    }
  } else {
    nabla <- lf$p
    phi <- 1
  }
  w <- diffC(mdl$y, nabla, mdl$bc)
  k <- length(par)
  if (is.null(fixed))
    fixed <- rep(FALSE, k)
  else {
    if (length(fixed) > k)
      stop("wrong dimension for fixed argument")
    else if (length(fixed) < k) {
      if (all(names(fixed) %in% names(par))) {
        f <- rep(FALSE, k)
        names(f) <- names(par)
        f[names(fixed)] <- fixed 
        fixed <- f
      } else stop("missing fixed parameters")
    }
  } 
  parfix <- par
  if (any(fixed)) {
    par <- par[!fixed]
    k <- length(par)
    bfixed <- TRUE
  } else bfixed <- FALSE

  if (!is.null(mdl$xreg)) {
    X <- sapply(1:ncol(mdl$xreg), function(col) diffC(mdl$xreg[,col], nabla, FALSE))
    a <- as.numeric(solve(t(X)%*%X, t(X) %*% w))
    names(a) <- paste0("a", 1:length(a))
    par <- c(par, a)
  }
  s2 <- c(1,0)
  if (k == 0) {
    ll0 <- .loglikrf(par, ...)
    mdl$S <- mdl$S*s2[1]
    return(mdl)
  }
  opt <- optim(par, .loglikrf, hessian = TRUE, ...)
  S <- mdl$S*s2[1]
  aux <- list()
  aux$n <- n
  res <- .res(opt$par, ...)
  J <- numDeriv::jacobian(.res, opt$par, ...)
  g <- t(J) %*% res
  varb <- solve(t(J) %*% J)*(sum(res^2)/length(res))
  mdl$S <- S
  par <- .sigmas(parfix, ...)
  aux$sigmas <- cbind(value = par[par != parfix],
                      ratio = par[par != parfix]/max(par[par != parfix]))
  if (length(par) >  length(aux$sigmas)) {
    par <- .sigmas(opt$par, ...)
    aux$par <- cbind(par = par[!(par != opt$par)],
                      se = sqrt(diag(varb)[!(par != opt$par)]))
  } else aux$par <- NULL
  a <- c(1, mdl$b)
  aux$s2 <- sum((a %*% S)*a)
  aux$g <- g
  aux$varb <- varb
  mdl$optim <- opt
  mdl$aux <- aux
  return(mdl)
}

#' Log-likelihood of an STS model
#' 
#' \code{logLik} computes the exact or conditional log-likelihood of object of 
#' the class \code{stsm}.   
#' 
#' @param mdl an object of class \code{stsm}.
#' @param y an object of class \code{ts}.
#' @param method exact or conditional.
#' @param tol tolerance to check if a root is close to one. 
#' @param ... additional parameters.   
#' 
#' @return 
#' The exact or conditional log-likelihood.
#' @export
logLik.stsm <-function(mdl, y = NULL, method = c("exact", "cond"), ...) {
  method <- match.arg(method)
  nabla <- nabla.stsm(mdl, ...)
  if (!is.null(y)) 
    w <- diffC(y, nabla$Pol, mdl$bc)
  else if(!is.null(mdl$y))
    w <- diffC(mdl$y, nabla$Pol, mdl$bc)
  else
    stop("missing time series");
  s2 <- 1
  if (!is.null(mdl$xreg)) {
    X <- sapply(1:ncol(mdl$xreg), function(col)
      diffC(mdl$xreg[,col], nabla, FALSE))
    ll <- llrfC(w - X %*% par[-(1:k)], nabla$Pol, mdl$b, mdl$C, mdl$S, s2)
  } else {
    ll <- llrfC(w, nabla$Pol, mdl$b, mdl$C, mdl$S, s2)
  }
  return(ll)
}

#' @export
AIC.stsm <- function(object, k = 2, ...) {
  ll <- -object$optim$value
  b <- object$optim$par
  aic <- (-2.0*ll+k*length(b))/object$aux$n
}

#' Initial conditions for Kalman filter
#'
#' \code{init} provides starting values to initialise the Kalman filter
#'
#' @param mdl an object of class \code{stsm}.
#' @param ... additional arguments.
#'
#' @return An list with the initial state and its covariance matrix.
#'
#' @export
init <- function (mdl, ...) { UseMethod("init") }

#' @rdname init
#' @export
init.stsm <- function(mdl, ...) {
  d <- length(mdl$b)
  y <- mdl$y[1:d]
  if (mdl$bc) y <- log(y)
  iC <- solve(mdl$C)
  biC <- mdl$b
  v <- c(1, rep(0, d))
  X <- mdl$b
  for (i in 2:d) {
    biC <- biC%*%iC
    X <- rbind(biC, X)
    v <- c(v, c(0, -biC))
  }
  V <- v
  for (i in 2:d) {
    v <- c(rep(0, d+1), v[1:((d+1)*(d-1))])
    V <- rbind(V, v)
  }
  I <- diag(1, d, d)
  S1 <- kronecker(I, mdl$S)
  V <- V %*% S1 %*% t(V)
  iX <- solve(X)
  V <- iX %*% V %*% t(iX)
  xd = iX %*% y
  list(xd = xd, Pd = V, X = X)
}

#' Initialization Kalman filter
#'
#' \code{ikf} computes the starting values x0 and P0 by generalized least
#' squares using the first n observations.
#' @param mdl an object of class \code{stsm}.
#' @param n number of observations used in the estimation.
#' @param ... additional arguments.
#'
#' @return An list with the initial state and its covariance matrix.
#'
#' @export
ikf <- function (mdl, ...) { UseMethod("ikf") }

#' @rdname ikf
#' @param y optional time series if it differes from model series.
#' @param n integer, number of observations used to estimate the inical
#'   conditions. By default, n is the length of y.
#' @export
ikf.stsm <- function(mdl, y = NULL, n = 0,  ...) {
  cform <- mdl$cform
  if (cform)
    mdl <- altform(mdl)
  d <- length(mdl$b)
  if (is.null(y))
    y <- mdl$y
  if (n <= d)
    n <- length(y)
  y <- y[1:n]
  if (mdl$bc) y <- log(y)
  if (!is.null(mdl$xreg))
    y <- y - mdl$xreg[1:n,] %*% mdl$a 
  x <- mdl$b
  X <- x
  for (i in 2:n) {
    x <- x %*% mdl$C
    X <- rbind(X, x)
  }
  x1 <- matrix(0, nrow = d, ncol = 1)
  P1 <- matrix(0, nrow = d, ncol = d)
  v <- matrix(0, nrow = n, ncol = 1)
  s2v <- matrix(0, nrow = n, ncol = 1)
  if(!kf0C(y, mdl$b, mdl$C, mdl$S, x1, P1, v, s2v))
    stop("Kalman filter algorithm failed")
  y <- v/sqrt(s2v)
  X <- sapply(1:d, function(j) { 
    kf0C(X[,j], mdl$b, mdl$C, mdl$S, x1, P1, v, s2v)
    v/sqrt(s2v)
  })
  iXX <- solve(t(X) %*% X)
  Xy <- t(X) %*% y
  b <- iXX %*% Xy
  s2 <- sum((y - X%*%b)^2)/(n-d)
  vb <- s2*iXX
  if (cform)
    list(x0 = b, P0 = vb)
  else
    list(x1 = b, P1 = vb)
}

#' Kalman filter for STS models
#'
#' \code{kf} computes the innovations and the conditional states with the Kalman
#' filter algorithm.
#'
#' @param mdl an object of class \code{stsm}.
#' @param ... additional arguments.
#'
#' @return An list with the innovations, the conditional states and their
#'   covariance matrices.
#'
#' @export
kf <- function (mdl, ...) { UseMethod("kf") }

#' @rdname kf
#' @param y time series to be filtered when it differs from the model series.
#' @param x1 initial state vector.
#' @param P1 covariance matrix of x1.
#' @param filtered logical. If TRUE, the filtered states x_{t|t} and their covariances
#'   matrices P_{t|t} are returned. Otherwise, x_{t|t-1} and P_{t|t-1} are 
#' @export
kf.stsm <- function(mdl, y = NULL, x1 = NULL, P1 = NULL, filtered = FALSE, ...) {
  cform <- mdl$cform
  if (cform) 
    mdl <- altform(mdl)
  d <- length(mdl$b)
  if (is.null(x1)||is.null(P1)) {
    ic <- ikf(mdl, y)
    if (is.null(x1)) 
      x1 <- ic[[1]]
    if (is.null(P1)) 
      P1 <- ic[[2]]
  }
  if (is.null(y))
    y <- mdl$y
  n <- length(y)
  if (mdl$bc) y <- log(y)
  if (!is.null(mdl$xreg))
    y <- y - mdl$xreg %*% mdl$a 
  X <- matrix(0, nrow = n, ncol = d)
  PX <- matrix(0, nrow = d*n, ncol = d)
  v <- matrix(0, nrow = n, ncol = 1)
  s2v <- matrix(0, nrow = n, ncol = 1)
  kfC(y, mdl$b, mdl$C, mdl$S, x1, P1, v, s2v, X, PX, cform, filtered, FALSE)
  logdet <- mean(log(s2v))
  s2 <- sum(v^2/s2v)/n
  ll <- -0.5*n*(1 + log(2*pi) + logdet + log(s2))
  v <- ts(v, end = end(mdl$y), frequency = frequency(mdl$y))
  X <- ts(X, end = end(mdl$y), frequency = frequency(mdl$y))
  list(ll = ll, se = sqrt(s2), X = X, PX = PX, u = v, s2u = s2v)
}

#' Kalman smoother for STS models
#'
#' \code{ks} computes smoothed states and their covariance matrices.
#'
#' @param mdl an object of class \code{stsm}.
#' @param ... additional arguments.
#'
#' @return An list with the smoothed states and their covariance matrices.
#'
#' @export
ks <- function (mdl, ...) { UseMethod("ks") }

#' @rdname ks
#' @param x1 initial state vector.
#' @param P1 covariance matrix of x1.
#' @export
ks.stsm <- function(mdl, x1 = NULL, P1 = NULL, ...) {
  cform <- mdl$cform
  if (cform) 
    mdl <- altform(mdl)
  if (is.null(x1)||is.null(P1)) {
    ic <- ikf(mdl)
    if (is.null(x1)) 
      x1 <- ic[[1]]
    if (is.null(P1)) 
      P1 <- ic[[2]]
  }
  d <- length(mdl$b)
  y <- mdl$y
  n <- length(y)
  if (mdl$bc) y <- log(y)
  if (!is.null(mdl$xreg))
    y <- y - mdl$xreg %*% mdl$a 
  X <- matrix(0, nrow = n, ncol = d)
  PX <- matrix(0, nrow = d*n, ncol = d)
  ksC(y, mdl$b, mdl$C, mdl$S, x1, P1, X, PX, cform)
  X <- ts(X, end = end(mdl$y), frequency = frequency(mdl$y))
  list(X = X, PX = PX)
}

#' Basic Structural Time Series models
#'
#' \code{bsm} creates/estimates basic structural models for seasonal time
#' series.
#'
#' @param y an object of class \code{ts}, with frequency 4 or 12.
#' @param bc logical. If TRUE logs are taken.
#' @param seas character, type of seasonality (Harvey-Durbin (hd), Harvey-Todd
#'   (ht), Harrison-Steven (ht))
#' @param par real vector with the error variances of each unobserved
#'   component.
#' @param fixed logical vector to fix parameters. 
#' @param xreg matrix of regressors.
#' @param fit logical. If TRUE, model is fitted.
#' @param updmdl function to update the parameters of the BSM.
#' @param ... additional arguments.
#'
#' @return An object of class \code{stsm}.
#'
#' @references
#'
#' Durbin, J. and Koopman, S.J. (2012) Time Series Analysis by State Space
#' Methods, 2nd ed., Oxford University Press, Oxford.
#'
#' @examples
#'
#' bsm1 <- bsm(AirPassengers, bc = TRUE)
#'
#' @export
#' 
bsm <- function(y, bc = FALSE, seas = c("hd", "ht", "hs"), 
                par = c(irr = 0.75, lvl = 1, slp = .05, seas = 0.075), 
                fixed = c(lvl = TRUE),
                xreg = NULL, fit = TRUE, updmdl = NULL, ...) {
  y.name <- deparse(substitute(y))
  s <- frequency(y)
  stopifnot(s > 1)
  seas <- match.arg(seas, c("hd", "ht", "hs"))
  if (is.null(updmdl))
    par <- log(par)
  else
    S <- updmdl(par, ...)$S
  if (seas == "ht") {
    b <- c(1, 0, 1, rep(0, s - 2))
    C <- matrix(0, s + 1, s + 1)
    C[1, 1] <- 1
    C[1, 2] <- 1
    C[2, 2] <- 1
    C[3, 3:(s+1)] <- -1
    diag(C[4:(s+1), 3:s]) <- 1
    if (is.null(updmdl)) {
      updmdl <- .updHT
      S <- .updHT(par, s)$S
    }
  } else if (seas == "hs") {
    b <- c(1, 0, 1, rep(0, s - 1))
    C <- matrix(0, s + 2, s + 2)
    C[1, 2] <- 1
    C[2, 2] <- 1
    diag(C[3:(s+1), 4:(s+2)]) <- 1
    C[s+2, 3] <- 1
    if (is.null(updmdl)) {
      updmdl <- .updHS
      S <- .updHS(par, s)$S
    }
  } else {
    X <- matrix(0, s+2, s+1)
    X[, 1] <- rep(1, s+2)
    X[, 2] <- 1:(s+2)
    t <- 1:(s+2)
    C <- sapply(1:floor(s/2), function(f) {
      cos(2*pi*f*t/s)
    })
    X[, seq(3, s+2, 2)] <- C
    C <- sapply(1:floor((s-1)/2), function(f) {
      sin(2*pi*f*t/s)
    })
    X[, seq(4, s+1, 2)] <- C
    C <- solve(X[1:(s+1),], X[2:(s+2), ])
    C[abs(C) < .0001] <- 0
    b <- X[1, ] %*% solve(C)
    b[abs(b) < .0001] <- 0
    if (is.null(updmdl)) {
      updmdl <- .updHD
      S <- .updHD(par, s)$S
    }
  }
  mdl <- stsm(y, b, C, S, xreg, bc, FALSE)
  mdl$y.name <- y.name
  if (fit)
    mdl <- fit.stsm(mdl, updmdl, par = par, fixed = fixed, .s = s, ...)
  class(mdl) <- c("stsm", "bsm")
  return(mdl)
}

.updHT <- function(par, .s, ...) {
  S <- matrix(0, .s+2, .s+2)
  diag(S)[1:4] <- exp(par[1:4])
  list(S = S)
}

.updHD <- function(par, .s, ...) {
  S <- matrix(0, .s+2, .s+2)
  diag(S) <- exp(c(par[1:3], rep(par[4], .s-1)))
  S[.s+2, .s+2] <- S[.s+2, .s+2]/2
  list(S = S)
}

.updHS <- function(par, .s, ...) {
  k <- .s+3
  S <- matrix(0, k, k)
  diag(S)[1:3] <- exp(par[1:3]) 
  S[4:k, 4:k] <- exp(par[4])*(diag(.s, .s) - matrix(1, .s, .s))/(.s-1)
  list(S = S)
}

.LeverrierFaddeev <- function(b, C) {
  n <- nrow(C)
  m <- ncol(C)
  if (n != m)
    stop("Argument 'C' must be a square matrix.")
  p <- rep(1, n + 1)
  I <- diag(1, n)
  C1 <- I
  A <- matrix(0, n+1, n)
  A[1, ] <- b
  for (k in 2:(n+1)) {
    p[k] <- -1 * sum(diag(C %*% C1))/(k - 1)
    C1 <- C %*% C1 + p[k] * I
    A[k, ] <- b %*% C1
  }
  A <- A[-n-1, ]
  if (is.vector(A)) A <- as.matrix(A)
  list(p = p, A = A)
}

#' @export
print.stsm <- function(x, ...) {
  if (!is.null(x$optim)) {
    if (!is.null(x$aux$par))
      print(x$aux$par)
    print(x$aux$sigmas)
    cat("\nlog likelihood: ", -x$optim$value)
    cat("\nResidual standard error: ", sqrt(x$aux$s2))
    cat("\naic:", AIC(x))
  } else {
    print("b")
    print(x$b)
    print("C")
    print(x$C)
    print("S")
    print(x$S)
  }
  return(invisible(NULL))
}

#' Reduce form for STS model
#'
#' \code{rform} finds the reduce form for a STS model.
#'
#' @param mdl an object of class \code{stsm}.
#' @param ... other arguments.
#'
#' @return An object of class \code{um}.
#'
#' @examples
#'
#' b <- 1
#' C <- as.matrix(1)
#' stsm1 <- stsm(b = b, C = C, Sv = c(lvl = 1469.619), s2u = c(irr = 15103.061))
#' rf1 <- rform(stsm1)
#' nabla(rf1)
#' theta(rf1)

#' @export
rform <- function (mdl, ...) { UseMethod("rform") }

#' @rdname rform
#' @param tol tolerance to check if a root is close to one.    
#' @export
rform.stsm <- function(mdl, tol = 1e-4, ...) {
  if (any(diag(mdl$S) < 0)) 
    stop("inadmissible model: negative variances")
  lst <- autocov(mdl, arma = FALSE, varphi = TRUE, lag.max = length(mdl$b))
  ma <- as.vector(autocov2MA(lst$g))
  r <- polyroot(lst$varphi)
  ar <- any(abs(abs(r)-1) > tol)
  if (ar) {
    ar <- roots2lagpol(r[abs(abs(r)-1) > tol])$Pol
    ar <- as.lagpol(ar, coef.name = "ar")
  } else {
    ar <- NULL
  }  
  i <- any((abs(r)-1) <= tol)
  if (i) {
    i <- roots2lagpol(r[(abs(r)-1) <= tol])$Pol
    i <- as.lagpol(i, coef.name = "i")
  } else {
    i <- NULL
  }  
  
  um1 <- um(bc = mdl$bc, ar = ar, i = i, ma = as.lagpol(ma[-1], coef.name = "ma")
     , sig2 = ma[1])
  um1$z <- mdl$y.name
  um1
}
