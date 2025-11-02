## tfarima/R/lagpol.R
## Jose L Gallego (UC)

#' Lag polynomials
#'
#' \code{lagpol} creates a lag polynomial of the form: 
#' \deqn{(1 - coef_1 B^s - ... - coef_d B^{sd})^p}
#' This class of lag polynomials is defined by:
#' \itemize{
#' \item the base lag polynomial \eqn{1 - coef_1 B^s - ... - coef_d B^{sd}},
#' \item the exponent `p` of the base lag polynomial (default is `p = 1`),
#' \item the spacing parameter `s` in sparse lag polynomials 
#' (default is `s = 1`),
#' \item the vector of `d` coefficients `c(coef_1, ..., coef_d)`, which can be 
#'  mathematical expresions dependent on `k` parameters 
#'  `c(param_1, ..., param_k)`.
#'  }
#' 
#' @param param a vector/list of named parameters. These parameters can be used 
#' within the coefficient expressions.
#' @param coef an optional vector of mathematical expressions defining the 
#' coefficients of the lag polynomial. If \code{NULL}, the names in \code{param} 
#' are used.
#' @param s an integer specifying the lag spacing or seasonal period.
#' @param p an integer specifying the exponent applied to the base lag polynomial.
#' @param lags an optional vector of lags for sparse polynomials. If \code{NULL},
#'  lags are determined by \code{s}.
#'
#' @return \code{lagpol} An object of class `lagpol` with the following
#'   components:
#' \describe{   
#'   \item{\code{coef}}{vector of coefficients c(coef_1, ..., coef_p) provided to 
#'   create the lag polynomial.}
#'   \item{\code{pol}}{base lag polynomial vector: 
#'   \eqn{1 - \text{coef}_1 B^s - \dots - \text{coef}_d B^{sd}}.}
#'   \item{\code{Pol}}{lag polynomial raised to the power `p`. 
#'   If `p = 1`, this equals `pol`.}
#'   }
#'   
#' @examples
#' # Simple AR(1) lag polynomial: 1 - 0.8B
#' lagpol(param = c(phi = 0.8))
#'
#' # AR(2) lag polynomial with seasonal lag s = 4: 1 - 1.2B^4 + 0.6B^8
#' lagpol(param = c(phi1 = 1.2, phi2 = -0.6), s = 4)
#'
#' # Integration operator squared: (1 - B)^2 = 1 - 2B + B^2
#' lagpol(param = c(delta = 1), p = 2)
#'
#' # Lag polynomial using explicit coefficients
#' lagpol(coef = c("1"), p = 2)  # (1 - B)^2 = 1 - 2B + B^2
#'
#' # Custom coefficients defined by mathematical expressions
#' lagpol(param = c(theta = 0.8), coef = c("2*cos(pi/6)*sqrt(theta)", "-theta"))
#' 
#' @export
lagpol <- function(param = NULL, s = 1, p = 1, lags = NULL, coef = NULL) 
{
  if(is.null(param) && is.null(coef)) {
    stop("0 arguments passed to 'lagpol' which requires 'param' or 'coef")
  }

  if ( is.numeric(param) ) param <- as.list(param)
  else if ( is.null(param) ) param <- list()
  else if (!is.list(param)) stop("'param' argument must be a list")
  
  param <- param[!duplicated(names(param))]
  names <- names(param)
  if ((is.null(names) || any(names == "")) && is.null(coef)) 
    stop("All parameters in 'param' must be named.")
  
  if (is.null(coef)) coef <- names(param)
  coef <- lapply(coef, function(x) parse(text = x))
  cc <- sapply(coef, function(x) eval(x, envir = param))
  
  s <- as.integer(s)
  p <- as.integer(p)
  if (s < 1 || p < 1) stop("'s' and 'p' must be >= 1.")  

  ncoef <- length(cc)
  if (!is.null(lags)) {
    lags <- unique(sort(lags))
    if (length(lags) != ncoef || any(lags <= 0)) stop("Invalid 'lags' vector.")    
  } else {
    lags <- (1:ncoef) * s
  }
  
  nlags <- length(lags)
  pol <- c(1, rep(0, lags[nlags]))
  pol[lags+1] <- -cc
  
  if (p > 1) Pol <- as.numeric( polyraiseC(pol, p) )
  else Pol <- pol
  
  names(pol) <- paste0("B^", 0:(length(pol)-1))
  names(Pol) <- paste0("B^", 0:(length(Pol)-1))
  
  lp <- list(pol = pol, Pol = Pol, d = length(pol) - 1, s = s, p = p, 
              lags = lags, nlags = nlags, param = param, 
              k = length(param), coef = coef, ncoef = ncoef)
  class(lp) <- c("lagpol")
  
  #if (!admreg.lagpol(lp, FALSE)) warning("Roots inside the unit circle")
  
  return(lp)
  
}

# Exported methods of lagpol

#' Lag polynomial converter
#' 
#' \code{as.lagpol} converts a numeric vector \code{c(1, -a_1, \dots, -a_d)} 
#' into a lag polynomial \eqn{(1 - a_1 B - ... - a_p B^p)}.
#' 
#' @param pol the numeric vector to be converted into an object of class lag 
#' polynomial.
#' @param p the exponent of the lag polynomial, positive integer.
#' @param coef.name name prefix for coefficients, character.
#' 
#' @return An object of class \code{lagpol}.
#' 
#' @examples
#' as.lagpol(c(1, -0.8))
#' as.lagpol(c(1, 0, 0, 0, -0.8))
#' @export 
as.lagpol <- function(pol, p = 1, coef.name = "a") {
  if (is.lagpol(pol)) return(pol)
  if (!is.numeric(pol)) return(NULL)
  pol <- as.numeric(pol)
  pol[1] <- 1
  k <- length(pol) - 1
  if (k == 0) return(NULL)
  if (!is.numeric(p) || length(p) != 1 || p <= 0 || p != as.integer(p)) {
    stop("'p' must be a positive integer.")
  }
  if (!is.character(coef.name) || length(coef.name) != 1) {
    stop("'coef.name' must be a character string.")
  }  
  pol <- as.numeric(pol)
  pol <- pol[-1]
  if (all(pol == 0)) {
    stop("'pol' must have at least one non-zero coefficient.")
  }  
  lags <- which(pol != 0)
  pol <- -pol[lags]
  names(pol) <- paste0(coef.name, lags)
  lp <- lagpol(pol, lags = lags, p = p)
  return(lp)
}

#' Inverse of a lag polynomial
#'
#' \code{inv} inverts a lag polynomial until the indicated lag.
#'
#' @param lp an object of class \code{lagpol}.
#' @param lag.max largest order of the inverse lag polynomial.
#' @param ... additional arguments.
#' 
#' @return \code{inv} returns a numeric vector with the coefficients
#' of the inverse lag polynomial truncated at lag.max.  
#' 
#' @examples
#' inv(as.lagpol(c(1, 1.2, -0.8))) 
#' 
#' @export
inv <- function (lp, ...) { UseMethod("inv") }

#' @rdname inv
#' @export
inv.lagpol <- function(lp, lag.max = 10, ...) {
  if (!is.numeric(lag.max) || lag.max < 1 || lag.max != as.integer(lag.max)) {
    stop("'lag.max' must be a positive integer.")
  }  
  if (is.lagpol(lp)) p <- lp$Pol
  else if(is.numeric(lp) && length(lp) > 1) {
    p <- lp
  } else stop("'lp' must be an object of class 'lagpol'.")

  lp1 <- tryCatch({
    as.numeric(polyratioC(1, p, lag.max))
  }, error = function(e) {
    stop("Error in computing the inverse: ", e$message)
  })  
    
  lp1 <- as.numeric( polyratioC(1, p, lag.max) )
  names(lp1) <- paste0("B^", 0:(length(lp1)-1))
  return(lp1)
}


#' Roots of lag polynomials
#'
#' \code{roots} computes the roots of lag polynomials from polynomial objects and 
#' time series models that contain lag polynomials as components.
#'
#' @param x A model object containing lag polynomials ("um", "tfm") or 
#'   a lag polynomial object ("lagpol").
#' @param ... Additional arguments passed to methods.
#' 
#' @return Returns a summary table with the roots of each \code{lagpol}.
#' 
#' @export
roots <- function(x, ...) { UseMethod("roots") }

#' @rdname roots
#' @export
roots.default <- function(x, ...) {
  if (is.lagpol.list(x)) {
    lapply(x, roots, ...)
  } else if (is.numeric(x) && length(x) > 1) {
    roots(as.lagpol(x), ...)
  } else {
    warning("x must be a lagpol object")
  }
}

#' @rdname roots
#' @param table Logical. If TRUE returns detailed table, if FALSE complex vector.
#' @param tol Tolerance for identifying distinct roots.
#' @examples
#' roots(c(1, 1.2, -0.8))
#' @export
roots.lagpol <- function(x, table = TRUE, tol = 1e-5, ...) {
  m <- 1
  if (is.lagpol(x)) {
    if (x$nlags == 1 && x$pol[x$d+1] == -1) {
      f <- 2*pi/x$d
      r <- sapply(0:(x$d-1), function(k) complex(1, cos(f*k), sin(f*k)))
    } else r <- polyroot(x$pol)
    m <- x$p
  } else if(is.numeric(x) & length(x) > 1) {
    r <- polyroot(x)
  } else {
    stop("error: x must be a lag polynomial")
  }
  
  if (table) { # to compare if two roots are almost equal
    .eq <- function (x, y, tol) { 
      if (all(abs(x - y) < tol)) return(TRUE)
      else(FALSE)
    }
    f <- acos( Re(r)/Mod(r) )/(2.0*pi)
    t <- cbind(Re(r), Im(r), Mod(r), f, 1/f, m)
    nr <- nrow(t)
    indx <- rep(TRUE, nr)
    if (nr > 1) { # remove repeated roots
      t <- t[order(t[, 4]), ] # roots ordered by frequency
      for (i in 1:nr) {
        if (indx[i] & i < nr) {
          for (j in (i+1):nr) {
            if ( .eq(t[j, 4], t[i, 4], tol) ) {
              if( .eq(t[j, 1:2], t[i, 1:2], tol) ) {
                t[i, 6] <- t[i, 6] + t[j, 6]
                indx[j] <- FALSE
              }
            } else {
              break
            }
          }
        }
      }
      t <- t[indx, ]
    }
    
    if (!is.matrix(t))
      t <- matrix(t, nrow = length(t)/6, ncol = 6)
    
    colnames(t) <- 
      c("Real", "Imaginary", "Modulus", "Frequency", "Period", "Mult.")
    return(t)
  } 
  
  return(rep(r, m))
  
}



#' Lag polynomial from roots
#'
#' \code{roots2lagpol} creates a lag polynomial from its roots.
#'
#' @param x a vector of real and/or complex roots.
#' @param lp logical. If TRUE, a object of class \code{lagpol} is returned;
#'   otherwise a numeric vector is returned.
#' @param ... additional arguments.   
#'
#' @return A lag polynomial or a numeric vector.
#'
#' @examples
#' roots2lagpol(polyroot(c(1, -1)))
#' roots2lagpol(polyroot(c(1, -1, 1)))
#'
#' @export
roots2lagpol <- function(x, lp = TRUE, ...) {
  if (!is.numeric(x) && !is.complex(x)) {
    stop("'x' must be a numeric or complex vector of roots.")
  }
  if (any(x == 0)) {
    stop("Roots cannot be zero; division by zero is undefined.")
  }  
  x <- 1/x
  pol <- 1
  for (r in x)
    pol <- c(pol, 0) - c(0, pol*r)
  pol <- Re(pol)
  if (lp)
    as.lagpol(pol, ...)
  else
    pol
}

#' Unit circle
#'
#' \code{unitcircle} plots the inverse roots of a lag polynomial together the
#' unit circle.
#'
#' @param lp an object of class \code{lagpol}.
#' @param s integer, seasonal period.
#' @param ... additional arguments.
#'
#' @return \code{unitcircle} returns a NULL value.
#'
#' @examples
#' unitcircle(as.lagpol(c(1, rep(0, 11), -1)))
#'
#' @export
unitcircle <- function (lp, ...) { UseMethod("unitcircle") }

#' @rdname unitcircle
#' @export
unitcircle.default <- function(lp, ...) {
  if (is.lagpol.list(lp)) {
    unitcircle.lagpol(lp, ...)
  } else if (is.numeric(lp)) {
    unitcircle(as.lagpol(lp), ...)
  } else {
    warning("lp must be a lagpol object")
  }
}

#' @rdname unitcircle
#' @export
unitcircle.lagpol <- function(lp, s = 12, ...) {
  t <- seq(0, 2*pi, length.out = 500)
  y <- sin(t)
  x <- cos(t)
  plot(x, y, type = "l", ylim = c(-1, 1), xlim = c(-1, 1), 
       xlab = "Real", ylab = "Imaginary", ...)
  if (is.lagpol.list(lp)) {
    k <- length(lp)
    for (i in 1:k) {
      r <- polyroot(rev(lp[[i]]$pol))
      x <- Re(r)
      y <- Im(r)
      points(x, y)
    }
  } else {
    r <- polyroot(rev(lp$pol))
    x <- Re(r)
    y <- Im(r)
    points(x, y)
  }
  if (s == 12 || s == 4) {
    r <- polyroot(c(-1, rep(0, s-1), 1))
    x <- Re(r)
    y <- Im(r)
    for (i in 1:s) {
      segments(0, 0, x[i], y[i], col = "gray")
    }
  }
  return(invisible(NULL))
}

#' Lag polynomial factorization
#'
#' \code{factors} extracts the simplifying factors of a polynomial in the lag
#' operator by replacing, if needed, its approximate unit or real roots to exact
#' unit or real roots.
#'
#' @param lp an object of class \code{lagpol}.
#' @param tol numeric tolerance for identifying real and unit roots (default is 
#' \code{1e-5}).
#' @param expand logical value to indicate whether the expanded polynomial
#'   should be returned.
#' @param ... additional arguments.
#'
#' @return \code{factors} returns a list with the simplifying factors of the lag
#'   polynomial or the expanded polynomial.
#'
#' @examples
#' factors( as.lagpol(c(1, rep(0, 11), -1)) )
#'
#' @export
factors <- function (lp, ...) { UseMethod("factors") }

#' @rdname factors
#' @param full logical value. If TRUE, the lag polynomial is completely
#'   factored. Otherwise, it is factored separating positive real roots from the
#'   others.
#' @param tol tolerance for nonzero coefficients.
#' @param expand logical value to indicate whether or not the factored lag
#'   polynomial must be expanded.
#' @export
factors.lagpol <- function(lp, full = TRUE, tol = 1e-5, expand = FALSE, ...) {
  if (is.numeric(lp)) 
    lp <- as.lagpol(lp)
  if (is.null(lp)) {
    if (expand) return(1)
    else return(NULL)
  }
  tol <- abs(tol)
  stopifnot(is.lagpol(lp))
  # Calcular raíces del polinomio
  t <- tryCatch({
    polyrootsC(lp$pol)
  }, error = function(e) {
    stop("Error in computing polynomial roots: ", e$message)
  })  
  
  if (nrow(t) > 1)
    t <- t[order(t[, 5], decreasing = TRUE), ]
  k <- 1
  p <- list()
  for (i in seq_len(nrow(t))) {
    x <- t[i, ]
    if (x[2] > -tol) {
      f <- x[4]
      param <- x[3]
      if(abs(abs(param)-1) < tol) 
        param <- sign(param)*1
      coef.name <- paste0("c", k)
      names(param) <- coef.name
      if (abs(f - 0.5) < tol) {
        coef <- c(paste0("-abs(", coef.name, ")"))
      } else if (f < tol) {
        f <- 0
        coef <- c(paste0("abs(", coef.name, ")"))
      } else {
        coef <- c(paste0("2*cos(2*pi*", f, ")*sqrt(abs(", coef.name, "))"),
                  paste0("-abs(", coef.name, ")"))
        param <- param^2
      }
      p[[k]] <- lagpol(param = param, coef = coef, p = x[6] + lp$p - 1)
      k <- k + 1
    } 
  }
  
  if (expand) return(polyexpand(p))
  
  if (!full && length(p) > 1) {
    i <- t[, 4] < tol
    if (any(i) && any(!i)) {
      p <- c(p[1:sum(i)], list(as.lagpol(polyexpand(p[-(1:sum(i))]))))
    } else if(all(!i)) {
      if (nrow(t) == 5 || nrow(t) == 13)
        p <- c(p[1], list(as.lagpol(polyexpand(p[-1]))))
      else
        p <- as.lagpol(polyexpand(p))
    }
  }
  return(p)
}

#' @title Print Method for Lag Polynomial Objects
#' @description Prints objects of class \code{lagpol}.
#' @param x An object of class \code{lagpol}.
#' @param digits Integer. Number of significant digits to display. Default is 2.
#' @param width Integer. Maximum number of characters per line. Default is the console width.
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns the \code{lagpol} object.
#' @export
print.lagpol <- function(x, digits = 2, width = getOption("width"), ...) {
  if (is.lagpol(x)) {
    if (x$d == 0) {
      print(1)
      return(invisible(x))
    }
    if (x$p > 1) {
      p <- paste0("(", as.character.lagpol(x, digits, TRUE), ")^",
                 x$p, " = ")
      p <- paste0(p, as.character.lagpol(x, digits, FALSE))
    } else {
      p <- as.character.lagpol(x, digits, FALSE)
    }
  } else p <- as.character.lagpol(x, digits, FALSE)
  
  pc <- nchar(p)
  m2 <- 0
  while(m2 < pc) {
    m1 <- m2 + 1
    m2 <- min(pc, m2 + width)
    if(m2 < pc)
      while(substring(p, m2, m2) != " " && m2 > m1 + 1)
        m2 <- m2 - 1
    cat(substring(p, m1, m2), "\n")
  }
  invisible(x)  
}

#' Print non-normalized polynomial as a lag polynomial
#'
#' Prints a non-normalized polynomial as a \code{lagpol} object, preserving the 
#' original coefficients.
#'
#' @param pol Numeric vector with the coefficients of a non-normalized polynomial. 
#'        The first element can be any numeric value, not necessarily 1.
#' @param digits Integer. Number of significant digits to display. Default is 2.
#' @param width Integer. Maximum number of characters per line. Default is the console width.
#' @return Invisibly returns the input vector.
#' @export
printLagpol <- function(pol, digits = 2, width = getOption("width")){
  print.lagpol(pol, digits = digits, width = width)
  invisible(pol)  
}

#' Prints a list of \code{lagpol} objects.
#'
#' @param llp A list of \code{lagpol} objects.
#' @param digits Integer. Number of significant digits to display. Default is 2.
#' @param width Integer. Maximum number of characters per line. Default is the console width.
#' @return Invisibly returns the list of \code{lagpol} objects.
#' @export
printLagpolList <- function(llp, digits = 2, width = getOption("width")){
  k <- length(llp)
  for (i in 1:k) {
    lp <- llp[[i]]
    if (is.lagpol(lp)) {
      if (lp$d == 0) {
         p <- "1"
      }
      if (lp$p > 1) {
        p <- paste0("(", as.character.lagpol(lp, digits, TRUE), ")^",
                   lp$p, " = ")
        p <- paste0(p, as.character.lagpol(lp, digits, FALSE))
      } else {
        p <- as.character.lagpol(lp, digits, FALSE)
      }
    } else p <- as.character.lagpol(lp, digits, FALSE)
    if (i > 1) txt <- paste0(txt, "   [", i, "] ", p)
    else txt <- paste0("[1] ", p)
  }
  
  pc <- nchar(txt)
  m2 <- 0
  while(m2 < pc) {
    m1 <- m2 + 1
    m2 <- min(pc, m2 + width)
    if(m2 < pc)
      while(substring(txt, m2, m2) != " " && m2 > m1 + 1)
        m2 <- m2 - 1
    cat(substring(txt, m1, m2), "\n")
  }
  invisible(llp)  
}

# Internal functions

#' @title Convert lag polynomial to character
#' @description Converts a \code{lagpol} object into its character
#'   representation.
#' @param x A \code{lagpol} object or a numeric vector representing the lag
#'   polynomial.
#' @param digits Integer. Number of significant digits to display. Default is 2.
#' @param pol Logical. If \code{TRUE}, uses the \code{pol} slot; otherwise, uses
#'   the \code{Pol} slot.
#' @param eq Logical. If \code{TRUE}, formats the output for equations. Default
#'   is \code{FALSE}.
#' @param eps Numeric. Threshold below which coefficients are considered zero.
#'   Default is \code{5.96e-08}.
#' @param ... Additional arguments (currently unused).
#' @return A character string representing the lag polynomial.
#' @note
#' Based on as.character.polynomial() function from Bill Venables and Kurt 
#' Hornik and Martin Maechler (2019). polynom: A Collection of Functions to 
#' Implement a Class for Univariate Polynomial Manipulations. 
#' R package version 1.4-0.
#' @keywords internal
#' @noRd
as.character.lagpol <- function(x, digits = 2, pol = FALSE, eq = FALSE, 
                                eps = 5.96e-08, ...) {
  if (!is.lagpol(x) && !is.numeric(x)) {
    stop("'x' must be either a lagpol object or a numeric vector.")
  }  
  if (is.lagpol(x)) {
    if (pol) p <- x$pol else p <- x$Pol
  } else p <- x
  p <- signif(p, digits = digits)
  d <- length(p) - 1
  names(p) <- 0:d
  p <- p[abs(p) > eps]
  if (length(p) == 0) return("0")
  
  signs <- ifelse(p < 0, "- ", "+ ")
  if (signs[1] == "+ ") signs[1] <- ""
  
  np <- names(p)
  p <- as.character(abs(p))
  p[p == "1" & np != "0"] <- ""
  
  if (eq) pow <- paste0("B'^", np, "*'")
  else pow <- paste0("B^", np)
  pow[np == "0"] <- ""
  pow[np == "1"] <- "B"
  paste(signs, p, pow, sep = "", collapse = " ")
}

# Internal function to check if an object is of class 'lagpol'
is.lagpol <- function(x) {
  return(inherits(x, "lagpol"))
}

# Internal function to check if an object is a list of 'lagpol' objects.
is.lagpol.list <- function(lpl) {
  if(!base::is.list(lpl)) return(FALSE)
  all(vapply(lpl, is.lagpol, logical(1)))
}

#' Common factors of two lag polynomials
#'
#' \code{common.factors} finds the common factors or the greatest common divisor
#' (GCD) of two lag polynomials.
#'
#' @param lp1,lp2 Lag polynomials or numeric vectors to be converted into
#'   \code{lagpol}.
#' @param tol Numeric tolerance to determine if coefficients are considered
#'   zero. Default is \code{1e-5}.
#' @param expand Logical. If \code{TRUE}, returns the expanded product of the
#'   common factors. Default is \code{FALSE}.
#'
#' @return A list of common factors or the expanded polynomial if \code{expand =
#'   TRUE}.
#'
#' @examples
#' p1 <- as.lagpol(c(1, -1))
#' p2 <- as.lagpol(c(1, 0, 0, 0, -1))
#' common.factors(p1, p2)
#'
#' @keywords internal
#' @noRd
common.factors <- function(lp1, lp2, tol = 1e-5, expand = FALSE) {
  if (!is.lagpol(lp1)) lp1 <- as.lagpol(lp1) 
  if (!is.lagpol(lp2)) lp2 <- as.lagpol(lp2) 
  if (is.null(lp1)||is.null(lp2)) {
    if (expand) return(1)
    else return(list())
  }
  f1 <- factors(lp1, tol)
  f2 <- factors(lp2, tol)
  l1 <- length(f1)
  l2 <- length(f2)
  lst <- list()
  for (i in seq_along(f1)) {
    pol <- f1[[i]]$pol
    for (j in seq_along(f2)) {
      r <- polydivC(f2[[j]]$pol, pol, TRUE, tol)
      if (all(abs(r) < tol, na.rm = TRUE)) {
        pol <- as.lagpol(pol, min(f2[[j]]$p, f1[[i]]$p))
        lst <- append(lst, list(pol))
        break
      }
    }
  }
  
  if (expand) 
    return(polyexpand(lst))
  else 
    return(lst)
}

#' Check if two lag polynomials have common roots
#'
#' This function checks whether two polynomials share at least one common root.
#'
#' @param lp1,lp2 Lag polynomials or numeric vectors to be converted into
#'   \code{lagpol}.
#' @param all logical. If TRUE, checks if all roots of the firs lag polynomial 
#' are present in the second one.  If FALSE (default), checks if at least one 
#' root is shared. 
#' @param tol Numeric tolerance to determine if the difference is considered
#'   zero. Default is \code{1e-5}.
#' @return A logical value
#'
#' @examples
#' x1 <- c(1, -1)
#' x2 <- c(1, -2, 1)
#'
#' hasCommonRoot(x1, x2) # Expected output: TRUE
#'
#' @noRd
#' @keywords internal
.hasCommonRoot <- function(lp1, lp2, all = FALSE, tol = 1e-5) {
  if (is.lagpol(lp1))
    lp1 <- lp1$pol
  else if (!is.numeric(lp1))
    stop("'lp1' argument must be a lagpol object")
  if (is.lagpol(lp2))
    lp2 <- lp2$pol
  else if (!is.numeric(lp2))
    stop("'lp2' argument must be a lagpol object")
  if (all) {
    if (length(lp1) <= length(lp2)) {
      r <- polydivC(lp2, lp1, TRUE, tol)
      return(all(abs(r) < tol))
    }
  } else {
    r <- polydivC(lp2, lp1, TRUE, tol)
    return(all(abs(r) < tol))
  }  
  return(FALSE) 
}

.parcounter <- function(llp, coef.name) {
  if (is.null(llp)) return(0)
  param <- lapply(llp, function(x) x$param)
  param <- unlist(param)
  if (is.null(param)) return(0)
  param <- param[!duplicated(names(param))]
  nm <- names(param)
  nm <- nm[startsWith(nm, coef.name)]
  if (length(nm) > 0) {
    i <- as.numeric(gsub(coef.name, "", nm)) 
    i <- i[is.numeric(i)]
    if (length(i) > 0) return(max(i))
    else return(0)
  } else {
    return(0)
  }
} 

## Internals models used in the estimation of ARIMA and related models
## admreg.lagpol, .update.lagpol

#' Admissible region for lag polynomial
#'
#' \code{admreg.lagpol} checks if the parameters of a lag polynomial take 
#' admissible values.
#'
#' @param lp An object of class \code{lagpol}.
#' @param ar Logical. If \code{TRUE}, roots must be outside the unit circle 
#' (AR process);  if \code{FALSE}, roots can also be on the unit circle 
#' (MA process).
#'
#' @return \code{TRUE} if the parameters are admissible, \code{FALSE} otherwise.
#'
#' @note For AR polynomials, roots must lie outside the unit circle. For MA 
#' polynomials, roots can be on or outside the unit circle.
#'
#' @examples
#' p1 <- lagpol(param = c(phi = 1.0))
#' admreg.lagpol(p1, ar = TRUE)
#' @keywords internal
#' @noRd
admreg.lagpol <- function(lp, ar = TRUE) {
  
  if(!all(is.finite(lp$pol))) return(FALSE)
  
  if (ar) {
    if (!abs(lp$pol[lp$d+1]) < 1) return(FALSE)
    if (lp$k == 1) return(TRUE) # first order
    if (lp$nlags == 2 && lp$lags[2] == lp$lags[1] + lp$s) { # second order
      phi1 <- -lp$pol[lp$lags[1]+1]
      phi2 <- -lp$pol[lp$lags[2]+1]
      if (abs(phi2) < 1 && phi2+phi1 < 1 && phi2-phi1 < 1) return(TRUE)
      else return(FALSE)
    }
  } else {
    if (abs(lp$pol[lp$d+1]) > 1) return(FALSE)
    if (lp$k <= 1) return(TRUE) # first order
    if (lp$nlags == 2 && lp$lags[2] == lp$lags[1] + lp$s) { # second order
      phi1 <- -lp$pol[lp$lags[1]+1]
      phi2 <- -lp$pol[lp$lags[2]+1]
      if (abs(phi2) > 1 || phi2 + phi1 > 1 || phi2 - phi1 > 1) return(FALSE)
      else return(TRUE)
    }
  }
  
  # general polynomial
  if (ar) {
    if( all( Mod(polyroot(lp$pol)) > 1) ) return(TRUE)
    else return(FALSE)
  } else {
    if( all( Mod(polyroot(lp$pol)) >= 1) ) return(TRUE)
    else return(FALSE)
  }
}

#' Expand a list of lag polynomials
#'
#' \code{polyexpand} multiplies two or more lag polynomials and returns the 
#' resulting polynomial.
#'
#' @param ... Multiple \code{lagpol} objects or a list of \code{lagpol} objects.
#' @param names Logical. If TRUE, put names to coeficcients.
#' @return A numeric vector representing the expanded polynomial.
#'
#' @examples
#' p1 <- as.lagpol(c(1, -0.8))
#' p2 <- as.lagpol(c(1, 0, 0, 0, -0.8))
#' polyexpand(p1, p2)
#'
#' @keywords internal
#' @noRd
polyexpand <- function(..., names = TRUE) {
  pols <- list(...)
  n <- length(pols)
  # check if ... is a list
  if (n == 1 && is.list(pols[[1]])) {
    pols <- pols[[1]]
    pols <- pols[!sapply(pols, is.null)]
    n <- length(pols)
    if (n == 0) return(1)
  }
  
  if (n == 0 || is.null(pols[[1]])) {
    return(1)
  }  
  
  pol <- if (is.lagpol(pols[[1]])) pols[[1]]$Pol else pols[[1]]  
  if (n == 1) {
    return(pol)
  }
  
  for(i in 2:n){
    if(is.lagpol(pols[[i]])){
      pol <- polymultC(pol, pols[[i]]$Pol)
    } else {
      pol <- polymultC(pol, pols[[i]])
    }
  }
  pol <- as.numeric(pol)
  if (names)
    names(pol) <- paste0("[B", 0:(length(pol)-1), "]")
  return(pol)
}

#' Update lag polynomial parameters
#'
#' \code{.update.lagpol} updates the parameters of a lag polynomial object.
#'
#' @param lp An object of class \code{lagpol}.
#' @param param A named numeric vector with the new parameter values.
#'
#' @return An updated object of class \code{lagpol}.
#'
#' @note This hidden function is called by the \code{fit.um} function to 
#' estimate multiplicative ARIMA(p, d, q) models.
#'
#' @keywords internal
#' @noRd
.update.lagpol <- function(lp, param) {
  lp$param[] <- param[names(lp$param)]
  cc <- sapply(lp$coef, function(x) eval(x, envir = lp$param))
  lp$pol[lp$lags+1] <- -cc
  
  if (lp$p > 1) lp$Pol <- as.vector( polyraiseC(lp$pol, lp$p) )
  else lp$Pol <- lp$pol 
  
  return(lp)
  
}

## Helper function to create lag polynomials

#' Create lag polynomial objects
#'
#' `lagpol0` is a flexible constructor for \code{\link{lagpol}} objects.
#' It accepts multiple input formats: polynomial orders \eqn{(d, s, p)}, 
#' literal equations (e.g., \eqn{1 - 0.8B}), or seasonal specifications for 
#' special lag polynomials like \eqn{AR/I/MA(s-1)}.
#'
#' @param op Polynomial specification in one of these formats:
#'   \itemize{
#'     \item Numeric vector \eqn{c(d, s, p)} specifying polynomial degree, 
#'           lag step (default 1), and exponent (default 1)
#'     \item Character equation (e.g., \code{"1 - 0.8B"})
#'     \item Seasonal period or frequency (e.g., \code{"12"} or \code{"1/12"})
#'     \item List combining multiple specifications
#'   }
#' @param type Operator type: \code{"ar"} (autoregressive), \code{"ma"} (moving
#'   average), or \code{"i"} (integration).
#' @param envir Environment for argument evaluation. Defaults to parent frame.
#'
#' @return List of \code{\link{lagpol}} objects.
#'
#' @examples
#' # AR(1) polynomial
#' lagpol0(op = 1, type = "ar")
#'
#' # From literal equation
#' lagpol0(op = "1 - 0.8B", type = "ar")
#'
#' # Multiple polynomials at once
#' lagpol0(op = list(1, "1 - 0.5B"), type = "ma")
#' 
#' # Seasonal polynomial
#' lagpol0(op = "12", type = "ar")
#' 
#' # Custom orders with seasonal component
#' lagpol0(op = c(2, 12, 1), type = "ar")
#'
#' @seealso \code{\link{lagpol}}
#' @export
lagpol0 <- function(op, type, envir = parent.frame()) {
  if (is.lagpol(op)) return(list(op))
  if (is.lagpol.list(op)) return(op)
  if (!is.environment(envir)) envir <- parent.frame() 
  
  if (!type %in% c("ar", "ma", "i")) 
    stop("'type' must be 'ar', 'ma' or 'i'.")
  if (!is.environment(envir)) envir <- parent.frame() 
  if (!is.list(op)) {
    if (length(op) == 1 || is.character(op)) op <- as.list(op)
    else op <- list(op) # c(1, 12)
  } 
  k <- length(op)
  nms <- names(op)
  if (!is.null(nms))
    nms <- nms[nms != ""]
  if (length(nms) < k) {
    if (k == 1) nms <- type
    else nms <- paste0(type, 1:k) 
    names(op) <- nms
  }
  llp <- list()
  for (i in seq_len(k)) {
    if (is.lagpol(op[[i]]))
      lp <- op[i]
    else if (is.numeric(op[[i]]))
      lp <- .lagpol1(op[i], type, envir)
    else if(is.character(op[[i]])) {
      if (grepl("B", op[[i]]))
        lp <- .lagpol2(op[i], type, envir)
      else
        lp <- .lagpol3(op[i], type, envir)
    } else stop("Invalid 'lagpol' specification")
    llp <- c(llp, lp)
  }
  names(llp) <- paste0(type, seq_along(llp))
  return(llp)
}

## Internal helper functions to create lag polynomials
# .lagpol1(), .lagpol2(), .lagpol3() 


#' Create lag polynomials from order specifications
#'
#' Generates lag polynomials from orders \code{c(d, s, p)} where \code{d} is
#' polynomial degree, \code{s} is lag step (default 1), and \code{p} is
#' exponent (default 1).
#'
#' @param orders Numeric vector or list of \code{c(d, s, p)} specifications.
#' @inheritParams .lagpol0
#'
#' @return List of \code{\link{lagpol}} objects.
#'
#' @examples
#' .lagpol1(orders = 2, type = "ar")
#' .lagpol1(orders = c(3, 2, 2), type = "ma")
#' .lagpol1(orders = list(phi = 2, Phi = c(1, 12)), type = "ar")
#'
#' @seealso \code{\link{.lagpol0}}
#' @keywords internal
#' @noRd
.lagpol1 <- function(orders, type, envir = parent.frame()) {
  if (!type %in% c("ar", "ma", "i")) 
    stop("'type' must be 'ar', 'ma' or 'i'.")
  if (!is.environment(envir)) envir <- parent.frame() 
  
  if (is.numeric(orders)) {
    if (orders[1] == 0) return(NULL)
    nms <- names(orders)
    orders <- list(orders)
    if (length(nms) == 1) names(orders) <- nms
  }
  nms <- names(orders)
  orders <- lapply(orders, function(x) {
    if (!is.numeric(x)) 
      stop("Invalid order specification for 'lagpol' object.")
    x <- c(x, 1, 1)[1:3]
    x[-1] <- as.integer(x[-1]) 
    if (any(c(x[1] < 0, x[-1] < 1)))
      stop("Invalid order specification for 'lagpol' object.")
    x
  })

  k <- length(orders)
  if (!is.null(nms))
    nms <- nms[nms != ""]
  if (length(nms) < k) {
    nms <- paste0(type, 1:k) 
    names(orders) <- nms
  }

  k <- 1
  llp <- lapply(orders, function(x) {
    d <- x[1]; s <- x[2]; p <- x[3]
    if (d < 1) {
      # Restricted frequency polynomial      
      param <- switch(type, ar = 0.5, ma = 0.6, i  = 1)
      coef <- paste0("-abs(", nms[k], k, ")")
      if (d != 0.5) {
        coef <- c(paste0("-2*cos(2*pi*", d, ")*sqrt(abs(", nms[k], k, "))"), 
                  coef)
      }
      names(param) <- paste0(nms[k], 1)
      lagpol(param = param, s = s, p = p, coef = coef)
    } else {
      # Standard polynomial      
      if (type == "ar") param <- as.list(c(rep(0.01, d - 1), 0.1))
      else if (type == "ma") param <- as.list(c(rep(0.02, d - 1), 0.2))
      else { 
        param <- as.list(1)
        p <- p + d - 1
      }
      if (grepl("\\d$", nms[k]) && length(param) > 1) 
        names(param) <- paste0(nms[k], ".", 1:length(param))
      else if (length(param) > 1) 
        names(param) <- paste0(nms[k], 1:length(param))
      else
        names(param) <- nms[k]
      if (type == "i") lp <- lagpol(coef = param, s = s, p = p)
      else lp <- lagpol(param = param, s = s, p = p)
    }
    k <<- k + 1
    lp
  })
  llp
}

#' Create lag polynomials from textual equations
#'
#' Parses character strings defining lag polynomial equations (e.g.,
#' \code{"1 - 0.5*B - 0.3*B^2"}) into \code{\link{lagpol}} objects.
#'
#' @param str Character vector of lag polynomial equations.
#' @inheritParams .lagpol0
#'
#' @return List of \code{\link{lagpol}} objects.
#' 
#' @examples
#' .lagpol2(str = "1 - 0.8*B - 0.3*B^2", type = "ar")
#' .lagpol2(str = c("1 + 0.5*B", "1 + 0.4*B + 0.2*B^2"), type = "ma")
#'
#' @seealso \code{\link{.lagpol0}}
#' @keywords internal
#' @noRd
.lagpol2 <- function(str, type, envir = parent.frame()) {
  if (!is.environment(envir)) envir <- parent.frame() 
  if (is.list(str)) str <- unlist(str)
  if (!is.character(str)) 
    stop("'str' must be a character vector.")
  if (!type %in% c("ar", "ma", "i")) 
    stop("'type' must be 'ar', 'ma' or 'i'.")
  nm <- names(str)
  if (is.null(nm)) nm <- type
  str <- gsub(" ", "", str)
  # Handle multiple polynomials: (...)(...) format  
  if (grepl("\\)\\(", str)) {
    str <- gsub("\\)\\(", "\\)|\\(", str)
    str <- strsplit(str, split = "\\|")[[1]]
    if (grepl("\\d$", nm) && length(str) > 1) 
      names(str) <- paste0(nm, ".", 1:length(str))
    else if (length(str) > 1)
      names(str) <- paste0(nm, 1:length(str))
    else
      names(str) <- nm
    
    llp <- list()
    for (i in 1:length(str))
      llp[i] <- .lagpol2(str[i], type = type, envir = envir)
    names(llp) <- names(str)
    return(llp)
  }
  
  # Remove outer parentheses  
  if (startsWith(str, "(") && endsWith(str, ")"))
    str <- substr(str, start = 2, stop = nchar(str) - 1)
  else if (startsWith(str, "("))
    str <- substr(str, start = 2, stop = nchar(str))
  # Check if first term is equal to 1
  i <- regexpr("\\+", str)[1]
  j <- regexpr("\\-", str)[1]
  if ( (j < 0) || (i < j && i > 0) )
    j <- i
  stopifnot(j > 1)
  txt <- substr(str, start = 1, stop = j-1)
  i <- eval(parse(text = txt), envir)
  stopifnot(as.integer(i) == 1)
  str <- substr(str, start = j, stop = nchar(str))

  # Exponent of the lag polynomial
  p <- 1
  i <- max(regexpr("\\)[[:alnum:][:punct:]]", str)[1])
  if (i > 1) {
    txt <- substr(str, i+1, nchar(str))
    if (!grepl("B", txt)) {
      if (startsWith(txt, "^")) txt <- substr(txt, 2, nchar(txt))
      p <- try(eval(parse(text = txt), envir), silent = TRUE)
      if (inherits(p, "try-error")) 
        stop(paste0("Invalid expression in 'str': ", str))
      str <- substr(str, 1, i-1)
    }
  }
  
  # Vector of coefficients  
  str <- gsub("\\*B", "B", str)
  str <- gsub("B\\^", "B", str)
  a <- list(1)
  lags <- 0
  repeat {
    i <- regexpr("B", str, fixed = TRUE)[1]
    if (i < 1) break
    if (!(startsWith(str, "+") || startsWith(str, "-")))
      stop("Invalid expression for lag polynomial")
    txt <- substr(str, start = 1, stop = i-1)
    if (nchar(txt) > 1) {
      b <- try(eval(parse(text = txt), envir), silent = TRUE)
      if (inherits(b, "try-error")) 
        stop("Invalid expression for lag polynomial")
    } else if (txt == "-")
      b <- -1
    else
      b <- 1
    a <- c(a, list(-b))
    if ((i+1) <= nchar(str)){
      str <- substr(str, start = i+1, stop = nchar(str))
      i <- regexpr("\\+", str)[1]
      j <- regexpr("\\-", str)[1]
      if ( (j < 0) | (i < j & i > 0) ) j <- i
      if (j > 1) 
        txt <- substr(str, start = 1, stop = j-1)      
      else if (j==1)
        txt <- "1"
      else {
        txt <- substr(str, start = 1, stop = nchar(str))
        j <- nchar(str)
      }
      lags <- c(lags, as.integer(txt))
      str <- substr(str, start = j, stop = nchar(str))
    } else {
      lags <- c(lags, 1)
      break
    }
  }
  
  lags <- lags[-1]
  a[[1]] <- NULL
  s <- lags[1]
  if (s > 1) {
    if (any(lags%%s != 0)) s <- 1
  }
  if (grepl("\\d$", nm) && length(a) > 1) 
    names(a) <- paste0(nm, ".", 1:length(a))
  else if (length(a) > 1)
    names(a) <- paste0(nm, 1:length(a))
  else
    names(a) <- nm
  if (type == "i") 
    lp <- list(lagpol(coef = a, lags = lags, s = s, p = p))
  else 
    lp <- list(lagpol(param = a, lags = lags, s = s, p = p))
  names(lp) <- nm
  return(lp)
}

#' Create special seasonal lag polynomials
#'
#' Generates special lag polynomials: (1) AR/I/MA of order \eqn{s-1} for seasonal
#' period \eqn{s}, or (2) AR/I/MA(2) with specific frequency \eqn{k/s}.
#' Frequencies 0 and 0.5 are also supported.
#'
#' @param str Character vector specifying seasonal period or frequency
#'   (e.g., \code{"12"}, \code{"(0:6)/12"}, \code{"0/12"}).
#' @inheritParams .lagpol0
#'
#' @return List of \code{\link{lagpol}} objects.
#'
#' @examples
#' .lagpol3(str = "12", type = "ar")
#' .lagpol3(str = "(0:6)/12", type = "ma")
#' .lagpol3(str = "0/12", type = "i")
#'
#' @seealso \code{\link{.lagpol0}}, \code{\link{lagpol}}
#' @keywords internal
#' @noRd
.lagpol3 <- function(str, type, envir = parent.frame()) {
  if (is.list(str)) str <- unlist(str)
  if (!is.character(str)) 
    stop("'str' must be a character vector.")
  if (!type %in% c("ar", "ma", "i")) 
    stop("'type' must be 'ar', 'ma' or 'i'.")
  if (!is.environment(envir)) envir <- parent.frame() 
  param <- switch(type, ar = 0.3, ma = 0.6, i  = 1)
  nm <- names(str)
  if (is.null(nm)) nm <- type
  f <- lapply(str, function(x) {
    val <- try(eval(parse(text = x), envir), silent = TRUE)
    if (inherits(val, "try-error")) 
      stop(paste("Invalid expression in 'str':", x))
    return(val)    
  })
  f <- unlist(f)

  # Determine seasonal period  
  if (length(str) == 1) {
    k <- gregexpr("\\)\\/", str)[[1]] # several frequencies, (0:6)/12
    if (length(k) == 1 && k[1] > 0) {
      txt <- substr(str, k[1]+2, nchar(str))
      s <- eval(parse(text = txt), envir)
      if (s < 2) s <- 1
    } else {
      k <- gregexpr("\\/", str)[[1]] # single frequency 0/12
      if (length(k) == 1 && k[1] > 0) {
        txt <- substr(str, k[1]+1, nchar(str))
        s <- eval(parse(text = txt), envir)
        if (s < 2) s <- 1
      } else s <- 1
    }
  } else s <- 1
  counter <- 1
  if (grepl("\\d$", nm) && length(f) > 1) 
    nm <- paste0(nm, ".")
  lp <- lapply(f, function(x) {
    coef <- paste0(nm, counter)
    names(param) <- coef
    counter <<- counter + 1
    if (x > 1) {
      x <- as.integer(x)
      if (x > 1) {
        coef <- paste0("-abs(", coef, "^(", 1:(x-1), "/", x, "))")
      } else {
        x <- abs(x)
        coef <- paste0("-abs(", coef, "^(", 1:(x-1), "/", x-1, "))")
      }
      lagpol(param = param, coef = coef)
    } else if (x > 0 && x < 0.5) {
      if (s > 1) {
        coef1 <- paste0(2*cos(2*pi*x), "*", coef, "^(1/", s, ")")
        coef2 <- paste0("-", coef, "^(2/", s, ")",  sep = "")
      } else {
        coef1 <- paste0(2*cos(2*pi*x), "*sqrt(", coef, ")")
        coef2 <- paste0("-", coef)
      }
      lagpol(param = param, coef = c(coef1, coef2))
    } else if (x == 0.5) {
      if (s > 1) {
        coef1 <- paste("-abs(", coef, ")^(1/", s, ")", sep = "")
      } else  {
        coef1 <- paste("-abs(", coef, ")", sep = "")
      }
      lagpol(param = param, coef = coef1)
    } else if (x == 0) {
      if (s > 1) {
        coef1 <- paste("abs(", coef, ")^(1/", s, ")", sep = "")
      } else  {
        coef1 <- paste("abs(", coef, ")", sep = "")
      }
      lagpol(param = param, coef = coef1)
    } else stop(paste("Invalid frequency."))
  })
  names(lp) <- paste0(nm, 1:length(lp))
  return(lp)
}


# Utilities for Lag Polynomials

#' Evaluate the k-th derivative of a polynomial at point z
#'
#' @param pol Numeric vector of polynomial coefficients in ascending order
#'   (pol[1] = constant term, pol[2] = coefficient of z, etc.)
#' @param z Numeric value where the polynomial (or its derivative) is evaluated.
#' @param k Integer. Derivative order (0 = original polynomial).
#'
#' @return Numeric value of the k-th derivative of P(z).
#' @examples
#' pol <- c(1, 2, 3, 4)  # P(z) = 1 + 2z + 3z² + 4z³
#' polyderivEvalR(pol, 2, 0)  # 49
#' polyderivEvalR(pol, 2, 1)  # 62
#' @export
polyderivEvalR <- function(pol, z, k = 0) {
  if (is.lagpol(pol))
      pol <- pol$Pol
  if (k == 0L) {
    d <- 0
    for (i in seq_along(pol))
      d <- d * z + pol[length(pol) - i + 1L]
    return(d)
  }
  
  deriv <- pol
  for (r in seq_len(k)) {
    if (length(deriv) <= 1L) return(0)
    deriv <- deriv[-1L] * seq_along(deriv[-1L])
  }
  
  d <- 0
  for (i in seq_along(deriv))
    d <- d * z + deriv[length(deriv) - i + 1L]
  d
}

