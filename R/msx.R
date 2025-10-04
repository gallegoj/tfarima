## tfarima/R/msx.R
## Jose L Gallego (UC)

#' Minimum signal extraction
#'
#' \code{msx} extracts a signal from a time series.
#'
#' @param mdl an object of class \code{\link{um}} or \code{\link{tfm}}.
#' @param ... additional arguments.
#'
#' @return An object of class \code{msx}.
#'
#' @export
msx <- function (mdl, ...) { UseMethod("msx") }

#' @rdname msx
#' @param ar autoregressive lag polynomial for the signal.
#' @param i integrated lag polynomial for the signal.
#' @param single logical. \code{TRUE} for single signal and \code{FALSE} for
#'   multiple signals.
#' @param canonical logical value to set or not the canonical requirement.
#' @param cwfact  method to compute the Cramer-Wold factorization: roots of the
#'   AGCF, Laurie/Wilson algorithms.
#' @param pfrac  method to compute the partial fraction decomposition of the
#'   spectrum: extended Euclidean algorithm (gcd) or solving linear system (solve).
#' @param tol tolerance for null or unit values.
#' @param envir the environment.
#'
#' @export
msx.um <- function(mdl, ar = NULL, i = NULL, single = FALSE, canonical = FALSE, 
                   cwfact = c("roots", "wilson", "best"), 
                   pfrac = c("gcd", "solve"),
                   tol = 1e-5, envir = parent.frame(), ...) {
  
  cwfact <- match.arg(cwfact)
  pfrac <- match.arg(pfrac)
  if (is.null(ar) && is.null(i))
    stop("error: missing ARI polynomial for signal")
  
  .fw <- function(w, num, den) {
    x <- 2*cos(2*pi*w)
    num <- sum(num*(x^(0:(length(num)-1))))
    den <- abs(sum(den*(x^(0:(length(den)-1)))))
    return(num/den)
  }

  # ARI lag polynomial for signal
  if (!is.null(mdl$ar)) {
    if (!is.null(ar)) {
      ar <- .lagpol0(ar, "ar", envir = envir)
      stopifnot(is.lagpol.list(ar))
      phi1 <- polyexpand(ar)
    } else phi1 <- NULL
    if (length(phi1) > mdl$p + 1)
      stop("error: wrong order for input AR polynomial.")
    if (length(phi1) == mdl$p + 1) {
      if (any(abs(phi1 - mdl$phi) > tol))
        stop("error: incompatible input AR polynomial.")
      phi2 <- NULL
    } else if (length(phi1) > 1) {
      r <- polydivC(mdl$phi, phi1, TRUE, tol)
      if (any(abs(r) > tol))
        stop("error: incompatible input AR polynomial.")
      phi2 <- as.vector(polydivC(mdl$phi, phi1, FALSE, tol))
    } else phi2 <- polyexpand(mdl$ar)
    if (!single) {
      phi <- lapply(ar, function(x) x$Pol)
      if (!is.null(phi2))
        phi <- c(phi, list(phi2))
    }
  } else {
    if (!is.null(ar))
      stop("error: incompatible input AR polynomial.")
    phi <- NULL; phi1 <- NULL; phi2 <- NULL
  }

  if (!is.null(mdl$i)) {
    if (!is.null(i)) {
      i <- .lagpol0(i, "i", envir = envir)
      stopifnot(is.lagpol.list(i))
      nabla1 <- polyexpand(i)
    } else nabla1 <- NULL
    if (length(nabla1) > mdl$d + 1)
      stop("error: wrong order for input I polynomial.")
    if (length(nabla1) == mdl$d + 1) {
      if (any(abs(nabla1 - mdl$nabla) > tol))
        stop("error: incompatible input I polynomial.")
      nabla2 <- NULL
    } else {
      if (!is.null(nabla1)) {
        r <- polydivC(mdl$nabla, nabla1, TRUE, tol)
        if (any(abs(r) > tol))
          stop("error: incompatible input I polynomial.")
        nabla2 <- as.vector(polydivC(mdl$nabla, nabla1, FALSE, tol))
      } else nabla2 <- polyexpand(mdl$i)
    }
    if (!single) {
      nabla <- lapply(i, function(x) x$Pol)
      if (!is.null(nabla2))
        nabla <- c(nabla, list(nabla2))
    }
  } else {
    if (!is.null(i))
      stop("error: incompatible input I polynomial.")
    nabla <- NULL; nabla1 <- NULL; nabla2 <- NULL
  }
  
  # Partial fraction:
  # u(x)/v(x) = r(x)/v(x) + q(x)
  # r(x) = u1(x)/v1(x) + ... + uk(x)/vk(x)
  u <- wold.pol(tacovC(1, mdl$theta, 1, mdl$q))
  v <- wold.pol(tacovC(1, polymultC(mdl$phi, mdl$nabla), 1, mdl$p + mdl$d))
  q <- as.vector(polydivC(u, v, FALSE, tol))
  r <- as.vector(polydivC(u, v, TRUE, tol))
  if (single) {
   vlist <- list(polyexpand(c(phi1, nabla1)))
   if (!is.null(phi2) || !is.null(nabla2))
     vlist <- c(vlist, list(polyexpand(list(phi2, nabla2))))
  } else {
    if (!is.null(phi) && !is.null(nabla)) vlist <- c(phi, nabla)
    else if (!is.null(phi)) vlist <- phi
    else vlist <- nabla
  }
  vlist <- lapply(vlist, function(x) {
    wold.pol(tacovC(1, x, 1, length(x) - 1))
  })
  k <- length(vlist)
  if (pfrac == "gcd")
    pf <- .pf_sum(r, vlist, tol = tol)
  else
    pf <- .pf_solve(r, vlist, tol)
  fw <- sapply(1:k, function(x) {
    num <- pf$ulist[[x]]; den <- vlist[[x]]
    w <- seq(0, 0.5, 0.0001)
    fw <- sapply(w, function(w) .fw(w, num, den))
    fwm <- min(fw)
    wm <- w[fw == fwm]
    c(wm[1], fwm) # if multiple minima, return the first
  })
  if (canonical) {
    pf$ulist <- lapply(1:k, function(x) {
      zeroes <- rep(0, length(vlist[[x]]) - length(pf$ulist[[x]]))
      c(pf$ulist[[x]], zeroes) - fw[2, x]*vlist[[x]]
    })
  }
  malist <- lapply(pf$ulist, function(x) cwfact(wold.pol(x, "p"), method = cwfact))
  cwfe <- sapply(malist, function(x) x[[4]])
  malist <- lapply(malist, function(x) x$th)

  is.admissible <- TRUE
  if (single) {
    models <- vector("list", k)
    j <- 1
    op <- malist[[j]]
    if (fw[2, j] < 0) {
      is.admissible <- FALSE
      warning(paste0("Spectrum for signal", j, " takes negative values\n"))
    }
    models[[1]] <- um(ar = ar, i = i, ma = as.lagpol(op/op[1], 1, "theta"),
                      sig2 = sign(op[1])*op[1]^2 * mdl$sig2, warn = FALSE)
    if (k == 2) {
      j <- 2
      op <- malist[[j]]
      if (fw[2, j] < 0) {
        is.admissible <- FALSE
        warning(paste0("Spectrum for signal", j, " takes negative values\n"))
      }
      if (is.null(ar)) ar <- mdl$ar
      else ar <- as.lagpol(phi2, 1, "phi")
      if (is.null(i)) i <- mdl$i
      else i <- as.lagpol(nabla2)
      models[[2]] <- um(ar = ar, i = i, ma = as.lagpol(op/op[1], 1, "theta"), 
                        sig2 = sign(op[1])*op[1]^2 * mdl$sig2, warn = FALSE)
    }
  } else {
    kar <- length(ar) # AR operators provided by the user
    kar1 <- length(phi) 
    ki <- length(i)
    models <- vector("list", k)
    for (j in 1:k) {
      op <- malist[[j]]
      if (fw[2, j] < 0) {
        is.admissible <- FALSE
        warning(paste0("Spectrum for signal", j, " takes negative values\n"))
      }
      if (j <= kar) {
        models[[j]] <- um(ar = ar[[j]], ma = as.lagpol(op/op[1], 1, "theta"),
                          sig2 = sign(op[1])*op[1]^2 * mdl$sig2, warn = FALSE)
      } else if (j <= kar1) {
        models[[j]] <- um(ar = as.lagpol(phi[[j]], 1, "phi"), 
                          ma = as.lagpol(op/op[1], 1, "theta"),
                          sig2 = sign(op[1])*op[1]^2 * mdl$sig2, warn = FALSE)
      } else {
        if (j - kar1 <= ki) i1 <- i[[j - kar1]]
        else i1 <- as.lagpol(nabla[[j - kar1]])
        models[[j]] <- um(i = i1, ma = as.lagpol(op/op[1], 1, "theta"),
                          sig2 = sign(op[1])*op[1]^2 * mdl$sig2, warn = FALSE)
      }
    }    
  }
  names(models) <- paste0("signal", 1:length(models))
  qs <- q
  if (canonical) {
    for (j in seq_along(malist)) 
      qs[1] <- qs[1] + fw[2, j]
  }
  if (length(qs) > 1)
    irreg <- um(ma = as.lagpol(qs[-1]), sig2 = qs[1]*mdl$sig2, warn = FALSE)
  else
    irreg <- um(sig2 = qs[1]*mdl$sig2, warn = FALSE)
  models[["noise"]] <- irreg
  if (irreg$sig2 < 0) {
    is.admissible <- FALSE
    warning(paste0("Model for noise is not admissible\n"))
  }
  pf1 <- lapply(1:k, function(x) {
    lst <- list(pf$ulist[[x]], vlist[[x]], fw[1, x], fw[2, x])
    names(lst) <- paste0(c("u", "v", "w", "fw"), x)
    return(lst)
    })
  names(pf1) <- paste0("pf", 1:k)
  pf1 <- c(frac = list(list(u  = u, v = v, q = q, r = r)), pf1)
  
  # pf <- list(u = u, v = v, q = q, r = r, f = f, 
  #            eps = elist, cwfe = cwfe)
  x <- c(list(um = mdl), models, pf = list(pf1), is.admissible = is.admissible)
  class(x) <- "msx"
  if (!is.admissible)
    warning("Non-admissible decomposition.")
  return(x)
}


#' Reconstruct original ARIMA model from msx decomposition
#'
#' Reconstructs the original ARIMA model by summing all component models
#' from an msx decomposition. Used to verify that the decomposition is correct.
#'
#' @param msx An object of class "msx" containing the decomposition of an 
#'   ARIMA model into unobservable component ARIMA models.
#' 
#' @return A "um" object representing the reconstructed original ARIMA model.
#' 
#' @examples
#' airl <- um(i = "(1-B)(1-B12)")
#' msx1 <- msx(airl, i = 2)
#' add_msx(msx1) # returns the airl model
#' 
#' @seealso \code{\link{add_um}}
#' @export
add_msx <- function(msx) {
  k <- length(msx) - 2
  stopifnot(k > 2)
  um1 <- msx[[2]] 
  for (i in 3:k)
    um1 <- add_um(um1, msx[[i]])
  return(um1)
}


#' @rdname print
#' @export
print.msx <- function(x, ...) {
  k <- length(x) - 2
  nms <- names(x)
  for (i in 1:k) {
    if (is.um(x[[i]])) {
      a <- letters[i]
      if (i == 1) z <- "Z"
      else if (i == k) z <- "N"
      else if (k == 3 && i == 2) z <- "S"
      else z = paste0("S", i-1)
      cat(nms[i], "\n")
      equation(x[[i]], FALSE, z = z, a = a)
      cat("\n")
    }
  }
  return(invisible(NULL))
}

#' @rdname wkfilter
#' @export
wkfilter.msx <- function(object, z = NULL, tol = 1e-5, envir = parent.frame(), ...) {
  k <- length(object) - 2
  stopifnot(k > 2)
  if (!is.ts(z)) z <- z.um(object[[1]], z, envir)
  start <- start(z)
  s <- frequency(z)
  X <- sapply(2:k, function(i) {
    wkfilter(object[[1]], object[[i]], z, "series", tol, NULL)
  })
  X <- ts(cbind(z, X), start = start, frequency = s)
  colnames(X) <- c("Series", paste0("Signal", 1:(k-2)), "Noise")
  return(X)
}




#' Extended Euclidean algorithm for polynomials
#'
#' \code{gcd_poly} computes the greatest common divisor of two polynomials with
#' real coefficients, and returns the Bézout coefficients \( u_1(x) \) and \(
#' u_2(x) \) such that: \deqn{v1(x) * u1(x) + v2(x) * u2(x) = gcd(v1, v2)}
#'
#' @param v1 numeric vector with the coefficients of v1(x).
#' @param v2 numeric vector with the coefficients of v2(x).
#' @param tol tolerance for zero coefficients.
#'
#' @return A list with the polynomials u1(x) and u2(x), and the gcd.
#'
#' @examples
#' v1 <- c(1, 1)   # 1 + x
#' v2 <- c(1, -1)  # 1 - x
#' res <- gcd_poly(v1, v2)
#'
#' @noRd
#' @keywords internal
.gcd_poly <- function(v1, v2, tol = 1e-10) {
  .clean_pol <- function(p, tol = 1e-10) {
    p[abs(p) < tol] <- 0
    if (any(p != 0)) {
      d <- max(which(p != 0))
      return(p[1:d])
    } else {
      return(0)
    }
  }
  
  .sub <- function(p1, p2, tol = 1e-10) {
    p1 <- .clean_pol(p1, tol)
    p2 <- .clean_pol(p2, tol)
    n <- max(length(p1), length(p2))
    p1 <- c(p1, rep(0, n - length(p1)))
    p2 <- c(p2, rep(0, n - length(p2)))
    return(.clean_pol(p1 - p2, tol))
  }
  
  r0 <- .clean_pol(v1, tol)
  r1 <- .clean_pol(v2, tol)
  s0 <- c(1); s1 <- c(0)
  t0 <- c(0); t1 <- c(1)
  while (any(abs(r1) > tol)) {
    q <- polydivC(r0, r1, FALSE, tol)
    r2 <- polydivC(r0, r1, TRUE, tol)
    s2 <- .sub(s0, polymultC(q, s1), tol)
    t2 <- .sub(t0, polymultC(q, t1), tol)
    r0 <- r1; r1 <- r2
    s0 <- s1; s1 <- s2
    t0 <- t1; t1 <- t2
  }
  
  gcd <- r0
  if (length(gcd) == 1 && abs(gcd) < tol) {
    warning("GCD is 0: v1 and v2 are not coprime.")
  } else if (abs(gcd[length(gcd)] - 1) > tol) {
    factor <- gcd[length(gcd)]
    s0 <- s0 / factor
    t0 <- t0 / factor
    gcd <- gcd / factor
  }
  return(list(u1 = as.vector(.clean_pol(s0, tol)),
              u2 = as.vector(.clean_pol(t0, tol)),
              gcd = as.vector(.clean_pol(gcd, tol))))
}

#' Partial fraction decomposition
#'
#' \code{.pf_sum} finds the partial fraction decomposition: \deqn{u(x)/v(x) =
#' u1(x)/v1(x) + ... + vn(x)*un(x)}.
#'
#' @param u numeric vector with the coefficients of u(x).
#' @param vlist list of n numeric vectors with the coefficients of \code{v1(x),
#'   ..., vn(x)}.
#' @param tol tolerance for zero.
#'
#' @return A list of numeric vectors \code{u1(x), ..., un(x)}.
#' @examples
#' v1 <- c(1, -1)     # 1 - x
#' v2 <- c(1, 0, 1)   # 1 + x^2
#' v3 <- c(1, 1)      # 1 + x
#' res <- .pf_sum(u = 1, list(v1, v2, v3))
#'
#' @seealso \code{\link{gcd_poly}} for the two-polynomial version.
#' @noRd
#' @keywords internal
.pf_sum <- function(u, vlist, tol = 1e-10) {
  stopifnot(u[1] > 0)
  n <- length(vlist)
  stopifnot(n > 0)
  # if (n == 1) 
  #   return(vlist)
  if (is.lagpol.list(vlist))
    vlist <- lapply(vlist, function(x) x$Pol)
  ulist <- list()
  J <- 1:n
  for (i in J) {
    p <- 1
    for (j in J[-i])
      p <- polymultC(p, vlist[[j]])
    res <- .gcd_poly(vlist[[i]], p, tol)
    ulist[[i]] <- res$u2
    if (abs(res$gcd - 1) > tol) 
      warning(paste0("fraction: ", i, " factor is not coprime"))
  }

  qlist <- list()
  if (length(u) > 1) {
    for (i in seq_along(ulist)) {
      pol <- polymultC(ulist[[i]], u)
      ulist[[i]] <- as.vector(polydivC(pol, vlist[[i]], TRUE, tol))
      qlist[[i]] <- as.vector(polydivC(pol, vlist[[i]], FALSE, tol))
    }
  }  else {
    ulist <- lapply(ulist, function(x) x*u)
  }
  return(list(ulist = ulist, qlist = qlist))
}

#' Partial fraction decomposition
#'
#' \code{.pf_solve} calculates the partial fraction decomposition of a 
#' rational function given its numerator, denominator, and factors.
#'
#' @param u A numeric vector representing the coefficients of the numerator
#'   polynomial in increasing powers of \(x\).
#' @param vlist A list of numeric vectors representing the coefficients of the
#'   factors used in the decomposition.
#' @param tol Numeric tolerance for determining numerical equality. Defaults 
#' to `1e-6`.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{ulist}{A list of numeric vectors \code{u1(x), ..., un(x)}.}
#'   \item{A}{The coefficient matrix used to solve for the partial fractions.}
#'   \item{b}{The right-hand side vector of the linear system used to solve for 
#'   partial fractions.}
#'   \item{x}{The solution vector, representing the coefficients of the 
#'   partial fractions.}
#' }
#'
#' @examples
#' # Partial fraction decomposition of (1 - 0.8B)(1 - 0.8B12)/(1-B)(1-B12)
#' # Factors: (1 - B)^2 and (1 + B + ... + B^11)
#' um1 <- um(i = "(1-B)(1-B12)", ma = "(1 - 0.8B)(1 - 0.8B12)")
#' pf <- .partialFractions(um1$ma, um1$i, list(c(1, -2, 1), rep(1, 12)))
#' print(pf)
#'
#' @noRd
#' @keywords internal
.pf_solve <- function(u, vlist, tol = 1e-10) {
  k <- length(vlist)
  stopifnot(k > 0)
  l <- 1
  for (i in 1:k) 
    l <- l + length(vlist[[i]]) - 1 
  A <- lapply(1:k, function(x) {
    lf <- length(vlist[[x]])
    f1 <- polyexpand(vlist[-x], names = FALSE)
    f1 <- c(f1, rep(0, l - length(f1) - 1))
    sapply(0:(lf-2), function(x) c(rep(0, x), f1[1:(l-x-1)]))
  })
  A <- do.call(cbind, A)
  b <- c(u, rep(0, l - length(u) - 1))
  x <- try(solve(A, b), silent = TRUE)
  
  if (inherits(x, "try-error")) {
    x <- try(ginv(A) %*% b, silent = TRUE)
    if (inherits(x, "try-error")) {
      stop("Partial fraction decomposition was not possible.")
    }
  }    
  
  l <- sapply(1:k, function(x) {
    length(vlist[[x]]) - 1
  })
  l <- c(0, cumsum(l))
  ulist <- lapply(1:k, function(j) {
    x[(l[j]+1):l[j+1]]
  })
  list(ulist = ulist, A = A, b = b, x = x)
}

#' Wold polynomial
#' 
#' Transforming a palindromic polymonial into a Wold polynomial/ 
#' Computing the Cramer-Wold factorization
#'
#' \code{wold.pol} can be used with three purposes:
#' 
#' (1) to transform a self-reciprocal or palindromic polynomial
#'  a_0 + a_1(B+F) + a_2(B^2+F^2) + ... + a_p(B^p+F^p) 
#'  into a Wold polynomial
#'  b_0 + b_1(B+F) + b_2(B+F)^2 + ... + b_p(B+F)^p;
#'  
#'  (2) to revert the previous transformation to obtain the palindromic
#'  polynominal from a Wold polynomial and
#'  
#'  (3) to compute the Cramer-Wold factorization: b(B+F) = c(B)c(F). 
#'
#' @param x numeric vector, coefficients of a palindromic or a Wold polynomial.
#' @param type character indicating the type of polynomial: (1) Wold polynomial, 
#' (2) Palindromic polynomial and (3) Cramer-Wold factor.
#' @param tol tolerance to check if an autocovariance is zero.
#' @return Numeric vector.
#'
#' @examples
#' wold.pol(c(6, -4, 1))

#' @export
wold.pol <- function(x, type = c("wold", "palindromic", "cramer-wold"), 
                     tol = 1e-5) {
  type <- tolower(type)
  type <- match.arg(type)
  if (startsWith("cramer-wold", type)) {
    if (x[1] == 0) 
      return(0)
    l <- length(x)
    while(abs(x[l]) < tol) {
      x <- x[-l]
      l <- l - 1
    }
    if (l == 1)
      return(x)
    if (l > 2)
      x <- wold.pol(x)
    r <- polyroot(x)
    n <- length(r)
    r <- r[order(abs(Arg(r)))]
    im <- FALSE
    r <- sapply(r, function(x) {
      if (abs(x) < tol^2) {
        im <<- !im
        if (im) return(c(1i, -1i))
        else return(c(-1i, 1i))
      }
      if (abs(x-2) < tol^2) return(c(1.0, 1.0))
      if (abs(x+2) < tol^2) return(c(-1.0, -1.0))
      d <- sqrt(x^2 - 4)
      r1 <- (x + d)/2
      r2 <- (x - d)/2
      if (Mod(r1) >= Mod(r2)) return(c(r1, r2))
      else return(c(r2, r1))
    })
    if (n > 1) r <- r[1, ]
    else r <- r[1]
    ma <- 1
    for (r1 in r)
      ma <- c(ma, 0) - c(0, ma/r1)
    ma <- Re(ma)
    ma <- c(x[l]/tail(ma, 1), ma)
    if (!is.finite(ma[1]))
      return(wold.pol(x[-l], "c"))
    else return(ma)
  }
  x <- as.numeric(x)
  n <- length(x)
  if (n < 3) return(x)
  if (startsWith("palindromic", type)) {
    b <- rep(0, n)
    b[1:2] <- x[1:2]
    for (i in 3:n) {
      for (j in 0:(i-1)) {
        if (i-2*j<1) break
        b[i-2*j] <- b[i-2*j] + x[i]*choose(i-1,j)
      }
    }
    return(b)
  }
  p0 <- rep(0, n)
  p1 <- rep(0, n)
  p2 <- rep(0, n)
  b <- rep(0, n)
  b[1:2] <- x[1:2]
  p0[1] <- 2
  p1[2] <- 1
  for (i in 3:n) {
    p2[1:i] <- c(0, p1[1:(i-1)]) - p0[1:i]
    b[1:i] <- b[1:i] + x[i]*p2[1:i]
    p0 <- p1
    p1 <- p2
  }
  return(b)
}

#' Cramer-Wold Factorization
#'
#' \code{cwfact} performs the Cramer-Wold factorization of the generating
#' autocovariance function of a pure moving average (MA) process, expressed as:
#' \deqn{g(x) = \theta(x)\theta(x^{-1})} where \deqn{g(x) = g_0 + g_1(x +
#' x^{-1}) + \dots + g_q(x^q + x^{-q})} and \deqn{\theta(x) = \theta_0 +
#' \theta_1 x + \dots + \theta_q x^q}
#'
#' The factorization can be computed by finding the roots of the polynomial
#' \eqn{g(x)}, or using the iterative Wilson (1969) algorithm as implemented by
#' Laurie (1981).
#'
#' @param g A numeric vector with the autocovariance coefficients \code{c(g0,
#'   g1, ..., gq)}.
#' @param th Optional numeric vector with initial values for the MA coefficients
#'   \eqn{\theta(x)}.
#' @param method A character string specifying the factorization method to use.
#'   Options are \code{"roots"} (default) and \code{"wilson"}.
#' @param tol A numeric tolerance for convergence (only used for \code{method =
#'   "wilson"}). Default is \code{1e-8}.
#' @param iter.max Maximum number of iterations for the Wilson method. Default
#'   is \code{100}.
#'
#' @return A numeric vector containing the moving average coefficients
#'   \code{c(theta_0, ..., theta_q)}.
#'
#' @details The implementation for \code{method = "laurie"} is a custom R
#' adaptation of Algorithm AS 175 from Laurie (1981).
#'
#' @examples
#' g <- autocov(um(ma = "1 - 0.8B"), lag.max = 1)
#' cwfact(g, method = "roots")
#' cwfact(g, method = "wilson")
#'
#' @references
#'
#' Wilson, G. T. (1969). Factorization of the covariance generating
#' function of a pure moving average process. \emph{SIAM Journal on Numerical
#' Analysis}, 6(1), 1–7.
#'
#' Laurie, D. P. (1981). Cramer-Wold Factorization. \emph{Journal of the Royal
#' Statistical Society Series C: Applied Statistics}, 31(1), 86–93.
#'
#' @export
cwfact <- function(g, th = NULL, method = c("roots", "wilson", "best"), 
                   tol = 1e-8, iter.max = 500) {
  l <- length(g)
  while(abs(g[l]) < tol) {
    g <- g[-l]
    l <- l - 1
  }
  stopifnot(l > 0)
  if (g[1] < 0) 
    method <- "roots"
  
  if (length(g) == 1) {
    return(list(th = sqrt(g), converged = NULL, iter = NULL, 
                res_norm = 0, method = NULL))
  }
  if (g[1] == 0) 
    return(c(0, 0))
  method <- match.arg(method)
  if (method == "roots"||method == "best") {
    cw1 <- .cwfact_roots(g, tol)
    cw2 <- .cwfact_roots_pal(g, tol)
    if (cw1$res_norm > cw2$res_norm) 
      cw1 <- cw2
    if (method == "roots")
      return(cw1)
  }  
  if (method == "wilson"||method == "best") {
    cw2 <- .Laurie(g, th, tol, iter.max)
    cw3 <- .Wilson(g, th, tol, iter.max)
    if (cw2$res_norm > cw3$res_norm) 
      cw2 <- cw3
    if (method == "wilson")
      return(cw2)
  } 
  
  if (cw1$res_norm < cw2$res_norm) return(cw1)
  else return(cw2)

}

#' @keywords internal
#' @noRd
.cwfact_roots <- function(g, tol = 1e-8) {
  l <- length(g)
  if (l == 2) 
    r <- -g[1]/g[2]
  else 
    r <- polyroot(wold.pol(g))
  r <- as.complex(r)
  n <- length(r)
  r <- r[order(abs(Arg(r)))]
  im <- FALSE
  r <- sapply(r, function(x) {
    if (abs(x) < tol) {
      im <<- !im
      if (im) return(c(1i, -1i))
      else return(c(-1i, 1i))
    }
    if (abs(x-2) < tol) return(c(1.0, 1.0))
    if (abs(x+2) < tol) return(c(-1.0, -1.0))
    d <- sqrt(x^2 - 4)
    r1 <- (x + d)/2
    r2 <- (x - d)/2
    if (Mod(r1) <= Mod(r2)) return(c(r1, r2))
    else return(c(r2, r1))
  })
  if (n > 1) r <- r[1, ]
  else r <- r[1]
  th <- 1
  for (r1 in r)
    th <- c(th, 0) - c(0, r1*th)
  th <- Re(th)
  s2 <- g[l]/th[l]
  th <- th*(sqrt(abs(s2)))*sign(s2)
  if (!is.finite(th[1]))
    return(.cwfact_roots(g[-l], tol = tol))
  gs <- sapply(1:l, function(x) sum(th[x:l]*th[1:(l+1-x)]))
  f <- gs - g
  return(list(th = th, converged = NULL, iter = NULL, 
              res_norm = sqrt(sum(f^2)), method = ".cwfact_roots"))
}

#' @keywords internal
#' @noRd
.cwfact_roots_pal <- function(g, tol = 1e-8) {
  l <- length(g)
  r <- polyroot(c(rev(g), g[-1]))
  r <- r[order(Mod(r))]
  r <- r[1:(l-1)]
  r <- r[order(abs(Arg(r)))]
  th <- 1
  for (r1 in r)
    th <- c(th, 0) - c(0, r1*th)
  th <- Re(th)
  s2 <- g[l]/th[l]
  th <- th*(sqrt(abs(s2)))*sign(s2)
  if (!is.finite(th[1]))
    return(.cwfact_roots_pal(g[-l], tol = tol))
  gs <- sapply(1:l, function(x) sum(th[x:l]*th[1:(l+1-x)]))
  f <- gs - g
  return(list(th = th, converged = NULL, iter = NULL, 
              res_norm = sqrt(sum(f^2)), method = ".cwfact_roots_pal"))
}


#' Cramer-Wold factorization based on the Bauer algorithm
#'
#' \code{.Bauer} is an R adaptation of the \code{BAUER} subroutine from
#' Algorithm AS 175.
#' 
#' @inheritParams cwfact
#'
#' @return A list with the MA coefficients and an error flag.
#'
#' @seealso \code{\link{cwfact}}.
#'
#' @keywords internal
#' @noRd
.Bauer <- function(g, iter.max = 500, tol = 1e-8) {
  q <- length(g)
  g0 <- g[1]
  g <- g/g0
  th <- c(g[-1], 0)
  th0 <- th
  for (iter in seq_len(iter.max)) {
    if (g[1] <= abs(g[q]))
      break
    r <- -th[1] / g[1]
    if (abs(r) > 1) 
      break
    for (i in seq_len(q-1)) {
      g[i] <- g[i] + r * th[i]
      th[i] <- th[i+1] + r * g[i+1]
    }
    if (max(abs(th - th0)) < tol) break
    th0 <- th
  }
  if (g[1] <= 0) th <- th
  else th = sqrt(g0)*g/sqrt(g[1])
  gs <- sapply(1:q, function(x) sum(th[x:q]*th[1:(q+1-x)]))
  f <- gs - g
  return(list(th = th, converged = iter < iter.max, iter = iter, 
              res_norm = sqrt(sum(f^2)), method = ".Bauer"))
}

#' Laurie Recursion (AS 175)
#'
#' \code{.Laurie} performs one recursion step of the Wilson factorization
#' algorithm. It is an R adaptation of the \code{recurs} subroutine from
#' Algorithm AS 175 (Laurie, 1981).
#'
#' @inheritParams cwfact
#'
#' @return A list with the updated MA coefficients and an error flag.
#'
#' @seealso \code{\link{cwfact}}.
#' @keywords internal
#' @noRd
.Laurie <- function(g, th = NULL, tol = 1e-8, iter.max = 500) {
  if (is.null(th)) {
    b <- .Bauer(g, iter.max, tol)
    th <- b$th
    th[1] <- abs(th[1])
  }
  stopifnot(length(g) == length(th))
  q <- length(g)
  for (iter in seq_len(iter.max)) {
    th1 <- .LaurieIter(g, th) 
    if (th1[1] <= 0) 
      break
    th1 <- th1 + 0.5 * th
    if (abs(th1[q]/th1[1]) > 1)
      th1[q] <- 1/th1[q]
    if (max(abs(th1 - th)/th[1]) < tol) break
    th <- th1
  }
  gs <- sapply(1:q, function(x) sum(th[x:q]*th[1:(q+1-x)]))
  f <- gs - g
  return(list(th = th, converged = iter < iter.max, iter = iter, 
              res_norm = sqrt(sum(f^2)), method = ".Laurie"))
}

#' @keywords internal
#' @noRd
.LaurieIter <- function(g, th) {
  q <- length(th) - 1
  m <- q + 1
  g[1] <- 0.5 * g[1]
  for (i in seq_len(q)) {
    if (th[1] <= 0) 
      return(th)
    s <- g[m] / th[1]
    g[1:m] <- g[1:m] - s*th[m:1]
    g[m] <- s
    s <- th[m] / th[1]
    th[1:m] <- th[1:m] - s*th[m:1]
    th[m] <- s
    m <- m - 1
  }
  if (th[1] <= 0)
    return(th)

  th[1] <- g[1] / th[1]
  for (i in 2:(q+1)) {
    s <- th[i]
    th[i] <- 0
    th[1:i] <- th[1:i] - s*th[i:1]
    th[i] <- th[i] + g[i]
  }
  return(th)
}

#' Wilson algorithm to compute MA parameters from autocovariances
#'
#' \code{.Wilson} performs the factorization of the covariance generating
#' function of a pure moving average process.
#'
#' @inheritParams cwfact
#'
#' @return A list with the estimated coefficients, convergence status, iteration
#'   count, and residual norm..
#'
#' @references
#'
#' Wilson, G. T. (1969). Factorization of the covariance generating
#' function of a pure moving average process. \emph{SIAM Journal on Numerical
#' Analysis}, 6(1), 1–7.
#'
#' @seealso \code{\link{cwfact}}.
#' @keywords internal
#' @noRd
.Wilson <- function(g, th = NULL, tol = 1e-8, iter.max = 500) {
  q <- length(g)
  g0 <- g[1]
  g <- g/g0
  if (is.null(th))
    th <- c(sqrt(abs(g[1])), rep(0, q-1))
  stopifnot(length(th) == q)
  for (iter in seq_len(iter.max)) {
    gs <- sapply(1:q, function(x) sum(th[x:q]*th[1:(q+1-x)]))
    f <- gs - g
    if (all(abs(f) < tol)) 
      break
    T1 <- sapply(1:q, function(x) c(th[x:q], rep(0, x - 1)))
    T2 <- sapply(1:q, function(x) c(th[x:1], rep(0, q - x)))
    sol <- try(solve(T1 + T2, f), silent = TRUE)
    if (inherits(sol, "try-error")) {
      th <- th - solve(T1 + T2 + 0.001 * diag(q), f)
    } else {
      th <- th - sol
    }    
    th <- th / sqrt(sum(th^2))    
  }
  th <- th*sqrt(g0)
  return(list(th = th, converged = iter < iter.max, iter = iter, 
              res_norm = sqrt(sum(f^2)), method = ".Wilson"))
}

#' @keywords internal
#' @noRd
.cwfact_check <- function(g, th) {
  g0 <- tacovC(1, th[-1], th[1], length(g)-1)
  max(abs(g - g0))
}

#' Refine roots using Newton-Raphson
#'
#'   \code{.refine_roots} improves the accuracy of the roots (real or complex)
#'   computed with the \code{polyroot} by using the Newton-Raphson method.
#'
#' @param roots Inexact roots (real or complex number).
#' @param coeffs Numeric vector of polynomial coefficients.
#' @param tol Convergence tolerance. Default is 1e-15.
#' @param max_iter Maximum number of iterations. Default is 50.
#'
#' @return The refined roots.
#'
#' @keywords internal
.refine_roots <- function(roots, coeffs, tol = 1e-15, max_iter = 50) {
  nr <- length(roots)
  n <- length(coeffs) - 1
  refined <- rep(FALSE, nr)
  px <- double(nr)
  for (i in 1:nr) {
    z <- roots[i]
    p <- coeffs[n+1]
    for(j in n:1)
      p <- p * z + coeffs[j]
    px[i] <- p 
    if (is.finite(p) && abs(p) > tol) {
      for(k in 1:max_iter) {
        p <- coeffs[n+1]
        dp <- 0
        for(j in n:1) {
          dp <- dp * z + p
          p <- p * z + coeffs[j]
        }
        if (!is.finite(p) || !is.finite(dp)) break
        if(abs(p) < tol || abs(dp) < tol) break
        z_new <- z - p / dp
        if (!is.finite(z_new)) break
        if(abs(z_new - z) < tol) break
        z <- z_new
      }
      if (is.finite(p) && abs(p) < abs(px[i])) {
        roots[i] <- z
        refined[i] <- TRUE
      } 
    }
  }    
  return(list(x = roots, px = px, refined = refined))
}

