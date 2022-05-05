## tfarima/R/stsm.R
## Jose L Gallego (UC)

#' Structural Time Series models
#'
#' \code{stsm} creates an S3 object representing a time-invariant structural
#' time series model.
#'
#' y_t = b'x_t + u_t (observation equation), 
#' x_t = Cx_t-1 + v_t (state equation).
#' 
#' @param y an object of class \code{ts}.
#' @param b vector of constants.
#' @param C matrix of constants.
#' @param s2u variance of the error u_t in the observation equation.
#' @param s2v variances of the vector error v_t in the state equation.
#' @param fSv function to create the covariance matrix of v_t.
#' @param xreg matrix of regressors.
#' @param bc logical. If TRUE logs are taken.
#' @param fit logical. If TRUE, model is fitted.
#' @param ... other arguments.
#'
#' @return An object of class \code{stsm}.
#'
#' @references
#'
#' Durbin, J. and Koopman, S.J. (2012) Time Series Analysis
# by State Space Methods, 2nd ed., Oxford University Press, Oxford.
#'
#' @examples
#' # Local level model
#' b <- 1
#' C <- as.matrix(1)
#' stsm1 <- stsm(Nile, b, C, s2v = c(lvl = 0.5), s2u = c(irr = 1))
#' stsm1
#' @export
#' 
stsm <- function(y, b, C, fSv, s2v, s2u = 1, xreg = NULL, bc = FALSE, fit = TRUE, ...) {
  dots <- list(...)
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
  if (missing(fSv)||is.null(fSv)) {
    stopifnot(length(s2v) == length(b))
    fSv <- .sDiag
    Sv <- fSv(s2v, length(b))
  } else Sv <- fSv(s2v, length(b), ...)
  if (!is.null(xreg))
    a <- solve(t(xreg) %*% xreg, t(xreg) %*% y)
  else
    a <- 0
  mdl <- list(y = y, y.name = y.name, bc = bc, b = b, C = C, Sv = Sv, s2u = s2u,
              fSv = fSv, s2v = s2v, xreg = xreg, a = NULL, dots = dots)
  class(mdl) <- "stsm"
  if (fit && !is.null(y))
    mdl <- fit.stsm(mdl)
  mdl
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
#' @param s2u variance of the error u_t.
#' @param s2v variances of the error vector v_t.
#' @param xreg matrix of regressors.
#' @param fSv function to create the covariance matrix of v_t.
#' @param ... other arguments.
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
                s2v = c(lvl = .2, slp = .05, seas = 0.075), s2u = .1,
                xreg = NULL, fSv = NULL,  ...) {
  s <- frequency(y)
  seas <- match.arg(seas, c("hd", "ht", "hs"))
  stopifnot(s != 4 || s!= 12)
  if (seas == "ht") {
    b <- c(1, 0, 1, rep(0, s - 2))
    C <- matrix(0, s + 1, s + 1)
    C[1, 1] <- 1
    C[1, 2] <- 1
    C[2, 2] <- 1
    C[3, 3:(s+1)] <- -1
    diag(C[4:(s+1), 3:s]) <- 1
    D <- diag(1, s+1, 3)
    if (is.null(fSv))
      fSv <- .sHT
  } else if (seas == "hs") {
    b <- c(1, 0, 1, rep(0, s - 1))
    C <- matrix(0, s + 2, s + 2)
    C[1, 2] <- 1
    C[2, 2] <- 1
    diag(C[3:(s+1), 4:(s+2)]) <- 1
    C[s+2, 3] <- 1
    D <- NULL
    if (is.null(fSv))
      fSv <- .sHS
  } else {
    t <- 1:(s+2)
    C1 <- cbind(rep(1, s+2), 1:(s+2), 
                cos(2*pi*1*t/s), sin(2*pi*1*t/s),
                if (s == 12) {
                  cbind(cos(2*pi*2*t/s), sin(2*pi*2*t/s),
                        cos(2*pi*3*t/s), sin(2*pi*3*t/s),
                        cos(2*pi*4*t/s), sin(2*pi*4*t/s),
                        cos(2*pi*5*t/s), sin(2*pi*5*t/s),
                        cos(2*pi*6*t/s))
                } else {
                  cos(2*pi*2*t/s)
                })
    C <- solve(C1[1:(s+1),], C1[2:(s+2), ])
    C[abs(C) < .0001] <- 0
    D <- NULL
    b <- C1[1, ] %*% solve(C)
    b[abs(b) < .0001] <- 0
    if (is.null(fSv))
      fSv <- .sHD
  }
  stsm(y, b, C, fSv, s2v, s2u, xreg, bc, TRUE, ...)
}

.sHT <- function(s2v, k, ...) {
  Sv <- diag(0, k, k)
  Sv[1,1] <- s2v[1]
  Sv[2,2] <- s2v[2]
  Sv[3, 3] <- s2v[3]
  Sv
}

.sHD <- function(s2v, k, ...) {
  Sv <- diag(0, k, k)
  Sv[1,1] <- s2v[1]
  Sv[2,2] <- s2v[2]
  if (length(s2v) == k) {
    diag(Sv)[3:k] <- s2v[3:k]
  } else if (length(s2v) > 3) {
      j <- 3
      for (i in 3:k) {
        Sv[i, i] <- s2v[floor(j)]
        j <- j + 0.5
      }
  } else {
    diag(Sv)[3:k] <- s2v[3]
    #Sv[k, k] <- 0.5*s2v[3]
  }
  Sv
}

.sHS <- function(s2v, k, ...) {
  Sv <- diag(0, k)
  Sv[1, 1] <- s2v[1]
  Sv[2, 2] <- s2v[2]
  s <- k - 2
  Sv[3:k, 3:k] <- s2v[3]*(diag(s, s) - matrix(1, s, s))/(s-1)
  Sv
}

.sDiag <- function(s2v, k, ...) {
  diag(s2v, k, k)
}

.sPsi <- function(s2v, k, ...) {
  s2v %*% t(s2v)
}

#' @rdname autocov
#' @examples
#' # Local level model
#' b <- 1
#' C <- as.matrix(1)
#' stsm1 <- stsm(b = b, C = C, s2v = c(lvl = 1469.619), s2u = c(irr = 15103.061))
#' autocov(stsm1)
#' 
#' @export
autocov.stsm <- function(mdl, ...) {
  lf <- .LeverrierFaddeev(mdl$b, mdl$C)
  d <- lf$p
  A <- lf$A
  k <- length(d)
  g <- rep(0, k)
  for (i in 1:k) {
    sum <- 0
    for (j in i:k) 
      sum <- sum + d[j]*d[j-i+1]
    g[i] <- sum
  }
  g = g*mdl$s2u
  k <- k - 1
  for (i in 1:k) {
    sum <- 0
    for (j in i:k) 
      sum <- sum + sum((A[j, ] %*% mdl$Sv) * A[j-i+1, ])
    g[i] <- g[i] + sum
  }
  g
}

#' Estimation of a STS model
#'
#' \code{fit} fits the stsm to the time series y.
#'
#' @param mdl an object of class \code{\link{stsm}}.
#' @param method  argument of the \code{optim} function.
#' @param show.iter logical value to show or hide the estimates at the different
#'  iterations.
#' @param ... other arguments.  
#'
#' @return An object of class "stsm" with the estimated variances.
#'
#' @examples
#' # Local level model
#' b <- 1
#' C <- as.matrix(1)
#' stsm1 <- stsm(Nile, b, C, s2v = c(lvl = 0.5), s2u = c(irr = 1), fit = FALSE)
#' stsm1 <- fit(stsm1, method = "L-BFGS-B")
#' @export
fit.stsm <- function(mdl, method = "BFGS", show.iter = FALSE, ...) {
  
  .s2i <- function(par) {
    if (!is.null(mdl$xreg))
      s2star <- llrfC(w - X %*% par[-(1:k)], lf[[1]], lf[[2]], mdl$Sv, mdl$s2u, TRUE)
    else
      s2star <- llrfC(w, lf[[1]], lf[[2]], mdl$Sv, mdl$s2u, TRUE)
    s2 <<- s2*s2star
    pos <- c(s2i, s2pos)
    s2i <<- pos[s2[pos] == max(s2)][1]
    s2 <<- s2/s2[s2i]
    s2pos <<- pos[pos != s2i]
    par[1:k] <- log(s2[s2pos])
    par
  }
  
  .upmdl <- function(par) {
    s2[s2pos] <<- exp(par[1:k])
    if (max(s2) > 1000) {
      return(FALSE)
    }
    mdl$Sv <<- do.call(mdl$fSv, c(list(s2[-s2k], k = lb), mdl$dots)) 
    mdl$s2u <<- s2[s2k]
    TRUE
  }
  
  .loglikrf <- function(par) {
    if(!.upmdl(par))
      return(ll0)
    if (!is.null(mdl$xreg))
      ll <- -llrfC(w - X %*% par[-(1:k)], lf[[1]], lf[[2]], mdl$Sv, mdl$s2u, FALSE)
    else
      ll <- -llrfC(w, lf[[1]], lf[[2]], mdl$Sv, mdl$s2u, FALSE)
    if (show.iter)
      print(c(-ll, s2))
    ll
  }
  
  lb <- length(mdl$b)
  if (is.null(names(mdl$s2u)))
    s2 <- c(mdl$s2v, irr = mdl$s2u)
  else
    s2 <- c(mdl$s2v, mdl$s2u)
  
  s2k <- length(s2)
  s2i <- (1:s2k)[s2 == max(s2)][1]
  s2 <- s2/s2[s2i]
  s2pos <- (1:s2k)[s2 > .Machine$double.eps][-s2i]
  par <- log(s2[s2pos])
  k <- length(par)
  
  lf <- .LeverrierFaddeev(mdl$b, mdl$C)
  w <- diffC(mdl$y, lf[[1]], mdl$bc)
  if (!is.null(mdl$xreg)) {
    X <- sapply(1:ncol(mdl$xreg), function(col) diffC(mdl$xreg[,col], lf[[1]], FALSE))
    a <- solve(t(X)%*%X, t(X) %*% w)
    names(a) <- paste0("a", 1:length(a))
    par <- c(par, a)
  }
  ll0 <- .loglikrf(par)
  opt <- optim(par, .loglikrf, hessian = TRUE, method = method)
  if (max(s2) > 1) {
    ll0 <- opt$value
    par <- .s2i(opt$par)
    opt <- optim(par, .loglikrf, hessian = TRUE, method = method)
  }
  if (!is.null(mdl$xreg)) {
    s2star <- llrfC(w - X %*% par[-(1:k)], lf[[1]], lf[[2]], mdl$Sv, mdl$s2u, TRUE)
    mdl$a <- par[-(1:k)]
  } else
    s2star <- llrfC(w, lf[[1]], lf[[2]], mdl$Sv, mdl$s2u, TRUE)
  s2[s2pos] <- exp(opt$par[1:k])
  s2 <- s2*s2star
  mdl$Sv <- do.call(mdl$fSv, c(list(s2[-s2k], k = lb), mdl$dots)) 
  mdl$s2u <- s2[s2k]
  mdl$s2v <- s2[-s2k]
  mdl$optim <- opt
  mdl$s2i <- s2i
  mdl
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

#' Estimation of a STS model by the method of moments
#'
#' \code{fit2autocov} fits a STS model to a vector of theoretical
#' autocovariances.
#'
#' @param mdl an object of class \code{stsm}.
#' @param g a vector of theoretical autocovariances (gamma[k], k= 0, ..., K).
#' @param method optimation method.
#' @param show.iter logical. If TRUE, estimates at each iteration are printed.
#' @param ... other arguments.
#'
#' @return An object of class \code{stsm}.
#'
#' @export
fit2autocov <- function (mdl, ...) { UseMethod("fit2autocov") }

#' @rdname fit2autocov
#' @examples
#' um1 <- um(Nile, i = 1, ma = 1)
#' g <- autocov(um1, lag.max = 1)
#' # Local level model
#' b <- 1
#' C <- as.matrix(1)
#' stsm1 <- stsm(Nile, b, C, s2v = c(lvl = 0.5), s2u = c(irr = 1), fit = FALSE)
#' stsm2 <- fit2autocov(stsm1, g)
#' stsm2
#' @export
#' 
fit2autocov.stsm <- function(mdl, g, method = "BFGS", show.iter = FALSE, ...) {
  .f <- function(par, ...) {
    mdl$Sv <<- do.call(mdl$fSv, c(list(s2v = exp(par), k = lb), mdl$dots)) 
    g0 <- autocov.stsm(mdl)
    f <- sqrt(sum((g0-g)^2))
    if (show.iter)
      print(c(f = f, par))
    f
  }
  k <- length(mdl$b) + 1
  if (k != length(g))
    stop(paste0("The number of autocovariances should be ", k))
  mdl$s2u <- unname(abs(g[k]))
  par <- log(mdl$s2v)
  lb <- length(mdl$b)
  opt <- optim(par, .f, method = method)
  nms <- c("irr", names(mdl$s2v))
  mdl$s2v <- unname(exp(opt$par))
  mdl$Sv <- do.call(mdl$fSv, c(list(s2v = mdl$s2v, k = length(mdl$b)), mdl$dots)) 
  mdl$opt <- opt
  lf <- .LeverrierFaddeev(mdl$b, mdl$C)
  w <- diffC(mdl$y, lf[[1]], mdl$bc)
  mdl$opt$logLik <- llrfC(w, lf[[1]], lf[[2]], mdl$Sv, mdl$s2u, FALSE)
  mdl$optim <- NULL
  names(mdl$s2v) <- nms[-1]
  names(mdl$s2u) <- nms[1]
  mdl
}

#' @export
print.stsm <- function(x, ...) {
  if (!is.null(x$optim)) {
    s2 <- c(x$s2v, x$s2u)
    print(cbind(Var = s2, Ratio = s2/s2[x$s2i]))
    cat("\nlog likelihood: ", -x$optim$value)
  } else if (!is.null(x$opt)) {
      s2 <- c(x$s2v, x$s2u)
      print(cbind(Var = s2, Ratio = s2/max(s2)))
      cat("\nlog likelihood: ", x$opt$logLik)
      cat("\nValue obj. func.: ", x$opt$value)
  } else {
    print("b")
    print(x$b)
    print("C")
    print(x$C)
    print("s2u")
    print(x$s2u)
    print("Sv")
    print(x$Sv)
  }
}

#' Reduce form for STS model
#'
#' \code{rform} finds the reduce form for a STS model.
#' 
#' @param mdl an object of class \code{stsm}.
#' @param ... other arguments.
#' 
#' @return 
#' An object of class \code{um}.
#' 
#' @examples 
#' 
#' b <- 1
#' C <- as.matrix(1)
#' stsm1 <- stsm(b = b, C = C, s2v = c(lvl = 1469.619), s2u = c(irr = 15103.061))
#' rf1 <- rform(stsm1)
#' nabla(rf1)
#' theta(rf1)

#' @export
rform <- function (mdl, ...) { UseMethod("rform") }

#' @rdname rform
#' @export
rform.stsm <- function(mdl, ...) {
  lf <- .LeverrierFaddeev(mdl$b, mdl$C)
  d <- lf$p
  A <- lf$A
  k <- length(d)
  g <- rep(0, k)
  for (i in 1:k) {
    sum <- 0
    for (j in i:k) 
      sum <- sum + d[j]*d[j-i+1]
    g[i] <- sum
  }
  g = g*mdl$s2u
  k <- k - 1
  for (i in 1:k) {
    sum <- 0
    for (j in i:k) 
      sum <- sum + sum((A[j, ] %*% mdl$Sv) * A[j-i+1, ])
    g[i] <- g[i] + sum
  }
  ma <- as.vector(acovtomaC(g))
  um(bc = mdl$bc, i = as.lagpol(d), ma = as.lagpol(c(1, -ma[-1])), sig2 = ma[1])
}
