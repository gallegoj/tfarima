## tfarima/R/lagpol.R
## Jose L Gallego (UC)

#' Lag polynomials
#'
#' \code{lagpol} creates a lag polynomial of the form \eqn{(1 - coef_1 B^s - ...
#' - coef_d B^sd)^p}. This class of lag polynomials is defined by a vector of d
#' coefficients c(coef_1, ..., coef_d), the powers s and p, and a vector of k
#' parameters c(param_1, ..., param_k). The vector c(coef_1, ..., coef_d) is
#' actually a vector of math expressions to compute the value of each
#' coefficient in terms of the parameters.
#'
#' @param param a vector/list of named parameters.
#' @param coef a vector of math expressions.
#' @param s the seasonal period, integer.
#' @param p the power of lag polynomial, integer.
#' @param lags a vector of lags for sparse polynomials.
#'
#' @return \code{lagpol} returns an object of class "lagpol" with the following
#'   components: \describe{ \item{coef}{Vector of coefficients c(coef_1, ...,
#'   coef_p) provided to create the lag polynomial.}
#'
#'   \item{pol}{Base lag polynomial, c(1, -coef_1, ..., -coef_d).}
#'
#'   \item{Pol}{Power lag polynomial when p > 1.} }
#'   
#' @examples
#' lagpol(param = c(phi = 0.8) )
#' lagpol(param = c(phi1 = 1.2, phi2 = -0.6), s = 4)
#' lagpol(param = c(delta = 1), p = 2)
#'
#' @export
lagpol <- function(param = NULL, s = 1, p = 1, lags = NULL, coef = NULL) 
{
  if ( is.numeric(param) ) param <- as.list(param)
  else if ( is.null(param) ) param <- list()
  else if (!is.list(param)) stop("'param' argument must be a list")
  param <- param[!duplicated(names(param))]
  names <- names(param)
  if (is.null(names) && is.null(coef)) stop("unnamed parameters")
  if (is.null(coef)) coef <- names(param)
  coef <- lapply(coef, function(x) parse(text = x))
  cc <- sapply(coef, function(x) eval(x, envir = param))
  s <- as.integer(s)
  p <- as.integer(p)
  
  ncoef <- length(cc)
  if (!is.null(lags)) {
    lags <- unique(sort(lags))
    if ( length(lags) != ncoef) stop("Invalid vector of lags")
  } else {
    lags <- (1:ncoef) * s
  }
  nlags <- length(lags)
  pol <- c(1, rep(0, lags[nlags]))
  pol[lags+1] <- -cc
  
  if (p > 1) Pol <- as.vector( polyraiseC(pol, p) )
  else Pol <- pol
  
  names(pol) <- paste("[B", 0:(length(pol)-1), "]", sep = "")
  names(Pol) <- paste("[B", 0:(length(Pol)-1), "]", sep = "")
  lp <- list(pol = pol, Pol = Pol, d = length(pol) - 1, s = s, p = p, 
              lags = lags, nlags = nlags, param = param, 
              k = length(param), coef = coef, ncoef = ncoef)
  class(lp) <- c("lagpol")
  
  if (!admreg.lagpol(lp, FALSE)) warning("Roots inside the unit circle")
  
  return(lp)
  
}

# Exported methods of lagpol

#' Lag polynomial
#' 
#' \code{as.lagpol} converts a numeric vector c(1, -a_1, \dots, -a_d) into 
#' a lag polynomial \eqn{(1 - a_1 B - ... - a_p B^p)}.
#' 
#' @param pol a numeric vector.
#' @param p integer power.
#' 
#' @return An object of class \code{lagpol}.
#' 
#' @examples
#' as.lagpol(c(1, -0.8))
#' as.lagpol(c(1, 0, 0, 0, -0.8))
#' @export 
as.lagpol <- function(pol, p = 1) {
  
  if (is.lagpol(pol)) return(pol)
  pol <- as.numeric(pol)
  stopifnot(pol[1] == 1)
  k <- length(pol) - 1
  if (k == 0) return(lagpol(NULL))
  pol[1] <- 0
  lags <- pol != 0
  pol <- -pol[lags]
  lags <- (0:k)[lags]
  names(pol) <- paste("a", lags, sep = "")
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
  if (is.lagpol(lp)) p <- lp$Pol
  else if(is.numeric(lp) & length(lp) > 1) p <- lp
  else stop("error: lp must be a lag polynomial")
  
  lp <- as.numeric( polyratioC(1, p, lag.max) )
  names(lp) <- paste("[B", 0:(length(lp)-1), "]", sep = "")
  return(lp)
}


#' Roots of a lag polynomial
#'
#' \code{roots.lagpol} computes the roots of a lag polynomial.
#'
#' @param x an object of class \code{lagpol}.
#' @param table logical. If TRUE, it returns a  five columns table showing the
#'   real and imaginary parts, the modulus, the frequency and the period of each
#'   root.
#' @param ... additional arguments.   
#'
#' @return A vector or a table.
#'
#' @examples
#' roots(c(1, 1.2, -0.8))
#'
#' @export
roots.lagpol <- function(x, table = TRUE, ...) {
  
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
    eq <- function (x, y, tol = 1e-5) { 
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
            if ( eq(t[j, 4], t[i, 4]) ) {
              if( eq(t[j, 1:2], t[i, 1:2]) ) {
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
    
    colnames(t) <- 
      c("Real", "Imaginary", "Modulus", "Frequency", "Period", "Mult.")
    return(t)
  } 
  
  return(rep(r, m))
  
}

#' @rdname roots.lagpol
#' @export
roots.default <- function(x, ...) {
  if (is.lagpol.list(x)) {
    lapply(x, roots)
  } else if (is.numeric(x)) {
    roots(as.lagpol(x))
  } else {
    warning("x must be a lagpol object")
  }
}

# Based on as.character.polynomial() function from
# Bill Venables and Kurt Hornik and Martin Maechler (2019). polynom: A
# Collection of Functions to Implement a Class for Univariate Polynomial
# Manipulations. R package version 1.4-0.
# #' @export
as.character.lagpol <- function(lp, digits = 2, pol = FALSE, eq = FALSE, ...) {
  
  if (is.lagpol(lp)) {
    if (pol) p <- lp$pol
    else p <- lp$Pol
  } else p <- lp
  p <- signif(p, digits = digits)
  d <- length(p) - 1
  names(p) <- 0:d
  p <- p[p != 0]
  
  signs <- ifelse(p < 0, "- ", "+ ")
  if (signs[1] == "+ ") signs[1] <- ""
  
  np <- names(p)
  p <- as.character(abs(p))
  p[p == "1" & np != "0"] <- ""
  
  if (eq) pow <- paste("B'^", np, "*'", sep = "")
  else pow <- paste("B^", np, sep = "")
  pow[np == "0"] <- ""
  pow[np == "1"] <- "B"
  #stars <- rep.int("*", length(p))
  #stars[p == "" | pow == ""] <- ""
  #p <- paste(signs, p, stars, pow, sep = "", collapse = " ")
  paste(signs, p, pow, sep = "", collapse = " ")
  
}

#' @export
print.lagpol <- function(x, digits = 2, ...){
  if (is.lagpol(x)) {
    if (x$d == 0) {
      print(1)
      return(invisible(x))
    }
    if (x$p > 1) {
      p <- paste("(", as.character.lagpol(x, digits, TRUE), ")^",
                 x$p, " = ", sep = "")
      p <- paste(p, as.character.lagpol(x, digits, FALSE), sep = "")
    } else {
      p <- as.character.lagpol(x, digits, FALSE)
    }
  } else p <- as.character.lagpol(x, digits, FALSE)
  
  pc <- nchar(p)
  ow <- max(35, getOption("width"))
  m2 <- 0
  while(m2 < pc) {
    m1 <- m2 + 1
    m2 <- min(pc, m2 + ow)
    if(m2 < pc)
      while(substring(p, m2, m2) != " " && m2 > m1 + 1)
        m2 <- m2 - 1
    cat(substring(p, m1, m2), "\n")
  }
  invisible(x)  
}

#' Print numeric vector as a lagpol object
#' 
#' @param pol numeric vectors with the coefficients of a normalized polynomial.
#' @param digits number of decimals.
#' 
#' @export
printLagpol <- function(pol, digits = 2){
  if (is.lagpol(pol)) {
    print(pol)
  } else {
    w0 <- pol[1]
    pol[1] <- 1
    pol <- as.lagpol(pol)
    pol$pol[1] <- w0
    pol$Pol[1] <- w0
    print.lagpol(pol, digits = digits)
  }
}

#' Print a list of lagpol objects
#' 
#' @param llp a list of \code{lagpol} objects.
#' @param digits number of decimals.
#' 
#' @export
printLagpolList <- function(llp, digits = 2){
  
  k <- length(llp)
  for (i in 1:k) {
    lp <- llp[[i]]
    if (is.lagpol(lp)) {
      if (lp$d == 0) {
        print(1)
        return(invisible(lp))
      }
      if (lp$p > 1) {
        p <- paste("(", as.character.lagpol(lp, digits, TRUE), ")^",
                   lp$p, " = ", sep = "")
        p <- paste(p, as.character.lagpol(lp, digits, FALSE), sep = "")
      } else {
        p <- as.character.lagpol(lp, digits, FALSE)
      }
    } else p <- as.character.lagpol(lp, digits, FALSE)
    if (i > 1) txt <- paste(txt, "   [", i, "] ", p, sep = "")
    else txt <- paste("[1] ", p, sep = "")
  }
  
  pc <- nchar(txt)
  ow <- max(35, getOption("width"))
  m2 <- 0
  while(m2 < pc) {
    m1 <- m2 + 1
    m2 <- min(pc, m2 + ow)
    if(m2 < pc)
      while(substring(txt, m2, m2) != " " && m2 > m1 + 1)
        m2 <- m2 - 1
    cat(substring(txt, m1, m2), "\n")
  }
  invisible(llp)  
}


# Non-exported and hidden functions

# Admissible region
#
# \code{admreg.lagpol} checks if the parameters of a lag polynomial take
# admissible values.
#
# @param lp an object of class \code{lagpol}.
# @param ar logical. If TRUE, roots must be outside unit circle; if FALSE, roots
#   can also be on unit circle.
#
# @return \code{TRUE} for admissible values and \code{FALSE} for not admissible
#   values.
#
# @note While roots of AR polynomials must lie inside of unit circle, those of
#   MA polynomials can lie on the unit circle.
#
# @examples
# p1 <- lagpol(param = c(phi = 1.0))
# pol <- admreg(p1, out = FALSE)
# 
admreg.lagpol <- function(lp, ar = TRUE) {
  
  if(!all(is.finite(lp$pol))) return(FALSE)
  
  if (ar) {
    if (lp$k == 1) {
      if (abs(lp$pol[lp$d+1]) < 1) return(TRUE)
      else return(FALSE)
    }
    phi1 <- -lp$pol[lp$lags[1]+1]
    if (lp$nlags == 1) { # first order
      if (abs(phi1) < 1.0) return(TRUE)
      else return(FALSE)
    }
    if (lp$nlags == 2 & lp$lags[2] == lp$lags[1] + lp$s) { # second order
      phi2 <- -lp$pol[lp$lags[2]+1]
      if (abs(phi2) < 1 & phi2+phi1 < 1 & phi2-phi1 < 1) return(TRUE)
      else return(FALSE)
    }
  } else {
    if (lp$k <= 1) {
      if (abs(lp$pol[lp$d+1]) > 1) return(FALSE)
      else return(TRUE)
    }
    phi1 <- -lp$pol[lp$lags[1]+1]
    if (lp$nlags == 1) { # first order
      if (abs(phi1) > 1.0) return(FALSE)
      else return(TRUE)
    }
    if (lp$nlags == 2 & lp$lags[2] == lp$lags[1] + lp$s) { # second order
      phi2 <- -lp$pol[lp$lags[2]+1]
      if (abs(phi2) > 1 | phi2 + phi1 > 1 | phi2 - phi1 > 1) return(FALSE)
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

is.lagpol <- function(lp) {
  return(inherits(lp, "lagpol"))
}


is.lagpol.list <- function(lpl) {
  if(!base::is.list(lpl)) return(FALSE)
  all( sapply(lpl, is.lagpol) )
}


# Lag polynomial factorization
#
# \code{polyfactors} extracts the simplifying factors of a polynomial
#  in the lag operator.
#
# @param pol an object of class \code{lagpol}.
# 
# @return \code{polyfactors} returns a list with the simplifying factors
# of the lag polynomial.  
# 
# @examples
# polyfactors( polyexpand( c(1, 1.2, -0.8), c(1, -2, 1) ) ) 
# 
factorize.lagpol <- function(lp) {
  t <- polyrootsC(lp$pol)
  f <- -1
  k <- 1
  p <- apply(t, 1, function(x) {
    if (f != x[4]) {
      f <<- x[4]
      param <- x[3]
      coef.name <- paste("f", k, sep = "")
      k <<- k + 1
      if (k == 0.5) {
        coef <- c(paste("-abs(", coef.name, ")", sep = ""))
      } else if (k < 0.0000001) {
        coef <- c(paste("abs(", coef.name, ")", sep = ""))
      } else {
        coef <- c(paste("-2*cos(2*pi*", f, ")*sqrt(abs(", coef.name, "))", 
                        sep = ""), paste("-abs(", coef.name, ")", sep = ""))
      }
      return(lagpol(param = param, coef = coef, p = x[6] + lp$p - 1))
    }
    else
      return(NULL)
  })
  
  p = p[!sapply(p, is.null)]
  
  return(p)
  
}


# Expand a list of lag polynomials
# 
# \code{polyexpand} multiplies two or more lag polynomials. 
# 
# @param ... list of lag polynomials.
#
# @ return A numeric vector
#
# @examples
# p1 <- as.lagpol(c(1, -0.8))
# p2 <- as.lagpol(c(1, 0, 0, 0, -0.8))
# pol <- polyexpand(p1, p2)
#
polyexpand <- function(...) {
  pols <- list(...)
  n <- length(pols)
  # check if ... is a list
  if (n==1) {
    if (inherits(pols[[1]], "list")) {
      pols <- pols[[1]]
      n <- length(pols)
      if (n == 0) return(1)
    } else if(is.null(pols[[1]]))
      return(1)
  }
  
  pol <- pols[[1]]
  if (is.lagpol(pol)){
    pol <- pol$Pol
  }
  
  if (n==1) {
    return(pol)
  }
  
  for(i in 2:n){
    if(is.lagpol(pols[[i]])){
      pol <- polymultC(pol, pols[[i]]$Pol)
    }
    else{
      pol <- polymultC(pol, pols[[i]])
    }
  }
  pol <- as.numeric(pol)
  names(pol) <- paste("[B", 0:(length(pol)-1), "]", sep = "")
  
  return(pol)
  
}


.parcounter <- function(llp, coef.name) {
  if (is.null(llp)) return(0)
  param <- lapply(llp, function(x) x$param)
  param <- unlist(param)
  if (is.null(param)) return(0)
  param <- param[!duplicated(names(param))]
  nm <- names(param)
  nm <- nm[startsWith(nm, coef.name)]
  if (!is.null(nm)) {
    max(as.numeric(gsub(coef.name, "", nm)))
  } else {
    0
  }
} 


# Update lag polynomials
# 
# \code{.update.lagpol} is a hidden function to update the values
# of the parameters of a lag polynomial.
# 
# @param lp an object of class \code{lagpol}..
# @param param new vector of parameters.
#
# @return an object of class \code{lagpol}.
#
# @note This hidden function is called by the \code{fit.um} function to 
# estimate multiplicative ARIMA(p, d, q) models.
.update.lagpol <- function(lp, param) {
  lp$param[] <- param[names(lp$param)]
  cc <- sapply(lp$coef, function(x) eval(x, envir = lp$param))
  lp$pol[lp$lags+1] <- -cc
  
  if (lp$p > 1) lp$Pol <- as.vector( polyraiseC(lp$pol, lp$p) )
  else lp$Pol <- lp$pol 
  
  return(lp)
  
}







## Internal helper functions to create lag polynomials
# .lagpol0(), .lagpol1(), .lagpol2(), .lagpol3() 

# Main helper function to create lag polynomials
#
# \code{.lagpol0} is an internal helper function to simplify the creation of lag
# polynomials providing for each lagpol:
# (1) the orders (d, s, p), 
# (2) the literal equation, e.g., 1 - 0.8B, 
# (3) the seasonal period s or frequency k/s for special lag polynomials 
# AR/I/MA(s-1) or AR/I/MA(2) with restricted frequency.
#
# @param op info about the lag polynomials, integers or character.
# @param type the type of operators: "ar", "i" or "ma".
# @param coef.name prefix name for the parameters.
# @param counter index for parameters.
#
# @return A list of lag polynomials.
#
# @note This hidden helper function is called by the um function to create
#   multiplicative ARIMA(p, d, q) models.
#
# @seealso \link{lagpol}
.lagpol0 <- function(op, type, coef.name, counter = 1) {
  
  if (is.character(op)) {
    op <- gsub(" ", "", op)
    if (regexpr("B", op)[1] > 1) {
      if (!startsWith(op, "(")) op <- paste("(", op, ")", sep = "")
      op <- .lagpol2(op, type, coef.name, counter)
    }
    else{ 
      op <- .lagpol3(op, type, coef.name, counter)
    }
  }
  
  if (is.lagpol(op)) return(list(op))
  if (is.lagpol.list(op)) return(op)
  
  if (is.numeric(op) || all(sapply(op, is.numeric))) {
    op <- .lagpol1(op, type, coef.name, counter)    
    if (is.null(op)) return(NULL)
    if (is.lagpol(op)) op <- list(op)
    return(op)
  }
  
  if (is.list(op)) {
    l <- lapply(op, function(x) {
      lp <- .lagpol0(x, type, coef.name, counter)
      if (type != "i") {
        if (!is.list(lp)) lp <- list(lp)
        counter <<- .parcounter(lp, coef.name) + 1
      }
      lp
    })
    
    if (!is.lagpol(l[[1]])) l <- unlist(l, recursive = FALSE)
    return(l)
  }
  stop("invalid lag operators")
}


# Helper function to create lag polynomials
#
# \code{.lagpol1} creates a list of lag polynomials from a list of orders c(d,
# s, p), where (d, s, p) are the orders of a lag polynomial in B^s of degree d
# raised to the power d, \eqn{(1 - coef_1 B^s - ... - coef_d B^sd)^p}. By
# default, s and p are equal to 1. The values and names of the parameters are
# assigned by the own function.
#
# @param orders a list containing the orders of one or several polynomials.
# @inheritParams .lagpol0
#
# @return A list of lag polynomials.
#
# @seealso \link{.lagpol0}
#
# 
.lagpol1 <- function(orders, type, coef.name = "", counter = 1) {
  
  if (is.numeric(orders)) {
    if (orders[1] == 0) return(NULL)
    orders <- list(orders)
  }
  if (coef.name == "") coef.name <- "a"
  if (counter < 0) counter <- abs(counter)
  j <- counter
  llp <- lapply(orders, function(x) {
    l <- length(x)
    d <- x[1]
    if (l == 1) {
      s <- 1
      p <- 1
    } else if (l == 2) {
      s <- as.integer(x[2])
      p <- 1
    } else if (l == 3) {
      s <- as.integer(x[2])
      p <- as.integer(x[3])
    } else stop("wrong orders for lagpol")
    if (d < 0||s < 1 || p < 1) stop("wrong orders for lagpol")
    if (d < 1) {
      if (type == "ar") param <- 0.5
      else if (type == "ma") param <- 0.6
      else if (type == "i") param <- 1
      else stop(paste("invalid ", type, " operator"))
      coef <- paste("-abs(", coef.name, j, ")", sep = "" )
      if (d != 0.5) {
        coef <- c( paste("-2*cos(2*pi*", d, ")*sqrt(abs(", coef.name, j, "))",
                         sep = "" ), coef )
      }
      names(param) <- paste(coef.name, j, sep = "" )
      j <<- j + 1
      lagpol(param = param, s = s, p = p, coef = coef)
    } else {
      if (type == "ar") param <- as.list(c(rep(0.01, d-1), 0.1))
      else if (type == "ma") param <- as.list(c(rep(0.02, d-1), 0.2))
      else if (type == "i") { 
        param <- as.list(1)
        p <- p + d - 1
      } else stop(paste("invalid ", type, " operator"))
      names(param) <- paste(coef.name, j:(j+length(param)-1), sep = "" )
      j <<- j + length(param)
      if (type == "i") lagpol(coef = param, s = s, p = p)
      else lagpol(param = param, s = s, p = p)
    }
  })
    
  return(llp)
  
}

# Helper function to create special lag polynomials
#
# \code{.lagpol2} creates a list of lag polynomials from their textual
# equations.
#
# @param str equations of of the lag polynomials, character.
# @inheritParams .lagpol0
#
# @return A list of lag polynomials.
#
# @seealso \link{.lagpol0}
#
.lagpol2 <- function(str, type, coef.name = "", counter = 1) {
  
  if (!startsWith(str, "1")) {
    if (coef.name == "") coef.name <- "a"
    l <- gregexpr("\\(", str)[[1]]
    r <- gregexpr("\\)", str)[[1]]
    n <- length(l)
    stopifnot(length(r) == n)
    txt <- list()
    start <- l[1]+1
    if (n > 1) { 
      for (i in 2:n) {
        if (l[i] > r[i-1]) {
          txt1 <- substr(str, start, l[i]-1)
          if (endsWith(txt1, ")")) txt1 <- substr(txt1, 1, nchar(txt1)-1)
          txt <- c(txt, txt1)
          start <- l[i]+1
        }
      }
    }
    n <- nchar(str)
    if (endsWith(str, ")")) n <- n-1
    txt <- c(txt, substr(str, start = start, stop = n))
    ll <- lapply(txt, .lagpol2)
    if (counter < 0) counter <- abs(counter)
    j <- counter
    ll <- lapply(ll, function(x) {
      param <- x[[1]]
      names(param) <- paste(coef.name, j:(j+length(param)-1), sep = "" )
      j <<- j + length(param)
      if (type == "i") 
        lagpol(coef = param, lags = x[[2]], s = x[[3]], p = x[[4]])
      else 
        lagpol(param = param, lags = x[[2]], s = x[[3]], p = x[[4]])
    })
    return(ll)
  }
  
  # Check if first term is equal to 1
  i <- regexpr("\\+", str)[1]
  j <- regexpr("\\-", str)[1]
  if ( (j < 0) | (i < j & i > 0) )
    j <- i
  stopifnot(j > 1)
  txt <- substr(str, start = 1, stop = j-1)
  i <- eval(parse(text = txt), .GlobalEnv)
  stopifnot(as.integer(i) == 1)
  str <- substr(str, start = j, stop = nchar(str))
  
  # Power of the lag polynomial
  i <- regexpr("\\)[1-9]", str)[1]
  if (i > 1) {
    txt <- substr(str, i+1, nchar(str))
    p <- eval(parse(text = txt), .GlobalEnv)
    str <- substr(str, 1, i-1)
  } else {
    i <- regexpr("\\)\\^[1-9]", str)[1]
    if (i > 1) {
      txt <- substr(str, i+2, nchar(str))
      p <- eval(parse(text = txt), .GlobalEnv)
      str <- substr(str, 1, i-1)
    } else {
      p <- 1
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
    if (!(startsWith(str, "+") | startsWith(str, "-")))
      stop("Invalid lag polynomial")
    txt <- substr(str, start = 1, stop = i-1)
    if (nchar(txt)>1)
      b <- eval(parse(text = txt), .GlobalEnv)
    else if (txt == "-")
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
    }
    else {
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
  
  return(list(a, lags, s, p))
  
}

# Helper function to create special lag polynomials
#
# \code{.lagpol3} creates a list of special lag polynomials: (1) AR/I/MA(s-1)
# polynomials or (2) AR/I/MA(2) polynomial with a specific frequency k/s, where
# s is the seasonal period.
#
# @param str seasonal perior or frequency, character.
# @inheritParams .lagpol0
#
# @return A list of lag polynomials.
#
# @seealso \link{.lagpol0}
#
.lagpol3 <- function(str, type, coef.name = "", counter = 1) {
  
  if (type == "ar") param <- 0.3
  else if (type == "ma") param <- 0.6
  else if (type == "i") param <- 1
  else stop(paste("invalid ", type, " operator"))
  
  if (coef.name == "") coef.name <- "a"
  f <- lapply(str, function(x) {
    eval(parse(text = x), .GlobalEnv)
  })
  f <- unlist(f)
  
  if (length(str) == 1) {
    k <- gregexpr("\\)\\/", str)[[1]] # several frequencies, (0:6)/12
    if (length(k) == 1 && k[1] > 0) {
      txt <- substr(str, k[1]+2, nchar(str))
      s <- eval(parse(text = txt), .GlobalEnv)
      if (s < 2) s <- 1
    } else {
      k <- gregexpr("\\/", str)[[1]] # single frequency 0/12
      if (length(k) == 1 && k[1] > 0) {
        txt <- substr(str, k[1]+1, nchar(str))
        s <- eval(parse(text = txt), .GlobalEnv)
        if (s < 2) s <- 1
      } else s <- 1
    }
  } else s <- 1
  
  lp <- lapply(f, function(x) {
    coef <- paste(coef.name, counter, sep = "" )
    names(param) <- coef
    counter <<- counter + 1
    if (x > 1) {
      x <- as.integer(x)
      if (x > 1) {
        coef <- paste("-abs(", coef, "^(", 1:(x-1), "/", x, "))", sep = "")
      } else {
        x <- abs(x)
        coef <- paste("-abs(", coef, "^(", 1:(x-1), "/", x-1, "))", sep = "")
      }
      lagpol(param = param, coef = coef)
    } else if (x > 0 && x < 0.5) {
      if (s > 1) {
        coef1 <- paste(2*cos(2*pi*x), "*", coef, "^(1/", s, ")", sep = "")
        coef2 <- paste("-", coef, "^(2/", s, ")",  sep = "")
      } else {
        coef1 <- paste(2*cos(2*pi*x), "*sqrt(", coef, ")", sep = "")
        coef2 <- paste("-", coef, sep = "")
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
    } else stop("invalid operator")
  })
  
  return(lp)
  
}
