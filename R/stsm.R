## tfarima/R/stsm.R
## Jose L Gallego (UC)

# S3 classes to handle Unobserved Components Time Series Models (UCM)
# uc: unobserved component (UC), building block for UCM
# ucm: unobserved component model  
# ssm: state space representation for UCM

#' Unobservable components (UC) for structural time series models
#'
#' \code{uc0} constructs unobservable components by specifying their state-space
#' representation matrices. This is the helper function used internally by
#' \code{\link{uc}} to create predefined components like trends, seasonals, and
#' cycles.

#' @param name character string specifying the component name (e.g., "trend",
#'   "seasonal")
#' @param param named vector or list containing parameter values. Names must
#'   match those used in symbolic expressions within matrices C and S
#' @param b numeric vector defining the component's contribution to the
#'   observation equation (Z matrix in state-space notation)
#' @param C numeric matrix or matrix with character expressions defining the
#'   transition matrix (T matrix). If NULL, defaults to identity matrix
#' @param S numeric matrix or matrix with character expressions defining the
#'   state error covariance matrix (Q matrix). If NULL, assumes no stochastic
#'   elements

#' @param name character, name of the UC.
#' @param param a vector or list with the parameters of the UC.
#' @param b vector of the UC in the observation equation.
#' @param C matrix of the UC in the transition matrix.
#' @param S covariance matrix of the error vector driving the UC.
#'
#' @return An S3 object of class \code{\link{uc}}.
#'
#' @note Matrices b, C and S can include symbolic expressions in terms of the 
#' parameters of the model. These parameters are defined in the \code{param} 
#' argument with their corresponding names. See the \code{\link{ssm}} function
#' for more information about b, C and S.
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
#'
#' # Local linear trend component
#' param <- c(s2_lvl = 0.05, s2_slp = 0.025)
#' b <- c(1, 0)
#' C <- matrix(c(1, 1, 0, 1), 2, 2, byrow = TRUE)
#' S <- matrix(c("s2_lvl", "0", "0", "s2_slp"), 2, 2)
#' trend <- uc0("trend", param, b, C, S)
#' 
#' # Cycle component
#' param <- c(phi = 0.8, per = 5, s2c = 0.01)
#' b <- c(1, 0)
#' C <- matrix(c("phi*cos(2*pi/per)", "phi*sin(2*pi/per)", 
#' "-phi*sin(2*pi/per)", "phi*cos(2*pi/per)"), 2, 2, byrow = TRUE)
#' S <- matrix(c("s2c*(1 - phi^2)", "0", "0", "s2c*(1 - phi^2)"))
#' cycle <- uc0("cycle", param, b, C, S)
#' 
#' @export
#' 
uc0 <- function(name, param, b, C = NULL, S = NULL) {
  if (missing(name) || !is.character(name) || length(name) != 1) {
    stop("'name' must be a single character string")
  }  
  if (missing(param)) {
    stop("argument 'param' is required")
  }  
  nms <- names(param)
  if (is.null(nms))
    stop("param must be a named vector")
  param <- as.list(param)
  
  if (is.null(S)) stop("argument S required")
  # irregular component
  if (is.null(b) && is.null(C)) {
    .b <- NULL
    .C <- NULL
    .S <- sapply(S, function(x) parse(text = x))
    lS <- length(.S)
    stopifnot(lS == 1)
    S <- sapply(.S, function(x) unname(eval(x, envir = param)))
    m <- 0    
  } else {
    m <- length(b)
    stopifnot(nrow(C) == m && ncol(C) == m)
    if (is.character(b)) {
      .b <- sapply(b, function(x) parse(text = x))
      b <- sapply(.b, function(x) unname(eval(x, envir = param)))
    } else .b <- NULL
    if (is.character(C)) {
      .C <- sapply(C, function(x) parse(text = x))
      C <- sapply(.C, function(x) unname(eval(x, envir = param)))
      C <- matrix(C, m, m)
    } else .C <- NULL
    
    lS <- 0
    if (is.character(S)) {
      .S <- sapply(S, function(x) parse(text = x))
      S <- sapply(.S, function(x) unname(eval(x, envir = param)))
      lS <- length(S)
      if (lS <= m) S <- diag(c(S, rep(0, m - lS)), m, m)
      else if (lS == m*m) S <- matrix(S, m, m)
      else stop("wrong S matrix")
    } else .S <- NULL
  }
  stopifnot(all(diag(S) >= 0) )
  param <- param[!duplicated(names(param))]
  comp <- list(name = name, param = param, m = m, b = b, C = C, S = S,
               .b = .b, .C = .C, .S = .S, lS = lS)
  class(comp) <- "uc"
  return(comp)
}

#' Unobservable components
#'
#' \code{uc} creates an S3 object representing a predefined UC (trend, seasonal,
#'  cycle, autoregressive and irregular component).
#'
#' @param type character. The type of component: trend, seasonal, cycle, AR or
#'   irregular.
#' @param s2 vector with the variances of the errors driving the UC.
#' @param s integer. Seasonal period for seasonal components.
#' @param form character. The form of the seasonal component (trigonometric or
#'   dummy seasonality).
#' @param per numeric. The period of the cycle component.
#' @param phi numeric. The damping factor of the cycle component.
#' @param counter integer. Counter for numbering phi cycle parameters.
#'
#' @return An object of class \code{\link{uc}}.
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
#' # Trend component
#' trend <- uc("trend", s2 = c(s2_lvl = 0.5, s2_slp = 0.025))
#' # Seasonal component
#' seas <- uc("seasonal", s2 = c(s2_seas = 0.5), s = 12, form = "trig")
#' # Cycle component
#' cycle1 <- uc("cycle", s2 = c(s2_c1 = 0.5), per = 5, phi = 0.8)
#' cycle2 <- uc("cycle", s2 = c(s2_c2 = 0.5), per = 10, phi = 0.8, counter = 2)
#'
#' @export
#' 
uc <- function(type, s2 = NULL, s = 12, form = c("trig", "dummy"), 
                  per = 5, phi = 0.9, counter = 1) {
  if (grepl("+", type, fixed = TRUE)[1])
    type <- unlist(strsplit(type, "+", fixed = TRUE))
  if (length(type) > 1) {
    lst <- lapply(type, function(x) uc(x, s2, s, form, per, phi, counter))
    names(lst) <- sapply(lst, function(x) x$name)
    return(lst)
  }
  type <- trimws(type)
  param <- NULL
  if (startsWith("trend", type) || startsWith(type, "trend")) {
    if (nchar(type) > 5) r <- as.integer(substr(type, 6, nchar(type)))
    else if (length(s2) > 0) r <- length(s2) - 1
    else r <- 1
    stopifnot(r >= 0)
    type <- "trend"
    if (!is.null(s2)) {
      if (length(s2) > r + 1)
        stop("too many variances")
    } else {
      s2 <- 0.5^(1:(r+1))
    }
    if (is.null(names(s2))) {
      l <- length(s2)
      if (l == 1) names(s2) <- "s2_lvl"
      else if (l == 2) names(s2) <- c("s2_lvl", "s2_slp")
      else names(s2) <- paste0("s2_trend", 0:(l-1))
    }
    b <- c(1, rep(0, r))
    C <- t(sapply(0:r, function(i) choose(0:r, i)))
    S <- names(s2)
    param <- s2
  } else if (startsWith("seasonal", type)) {
    stopifnot(s > 1)
    form <- match.arg(form)
    C <- matrix(0, s-1, s-1)
    S <- matrix(0, s-1, s-1)
    if (is.null(s2)) 
      s2 <- c(s2_seas = 0.05)
    if (form == "trig") {
      b <- rep(c(1, 0), floor(s/2))[1:(s-1)]
      if (s > 2) {
        for (i in 1:floor((s-1)/2)) {
          j <- i*2 - 1
          C[j+1, j+1] <- C[j, j] <- cos(2*pi*i/s)
          C[j, j+1] <- sin(2*pi*i/s)
          C[j+1, j] <- -C[j, j+1]
        }
      }
      if (s %% 2 == 0) C[s-1, s-1] <- -1
      l <- length(s2)
      if (l == 1) {
        if (is.null(names(s2))) names(s2) <- "s2_seas"
        S <- rep(names(s2), s-1)
        if (s %% 2 == 0 && s > 2) S[s-1] <- paste0(0.5, "*", names(s2))
      } else if (l == floor(s/2)) {
        if (is.null(names(s2))) 
          names(s2) <- paste0("s2_seas", 1:l)         
        S <- rep(names(s2), each = 2)[1:(s-1)]         
      } else if (l == (s-1)) {
        if (is.null(names(s2)))
          names(s2) <- paste0("s2_seas", 1:(s-1))
        S <- names(s2)         
      } else {
        stop("wrong number of sigmas")
      }
    } else if (form == "dummy") {
      stopifnot(length(s2) == 1)
      if (is.null(names(s2)))
        names(s2) <- "s2_seas"
      b <- rep(0, s-1)
      b[1] <- 1
      C[1, ] <- -1
      if (s > 2) {
        for (i in 2:(s-1))
          C[i, i-1] <- 1
      }
      S <- names(s2)
    } else stop("unknown seasonal component")
    param <- s2
  } else if (startsWith("ar", type) || startsWith(type, "ar")) {
    if (nchar(type) > 2) r <- as.integer(substr(type, 3, nchar(type)))
    else r <- length(phi)
    if (is.null(s2)){
      s2 <- 0.5
      names(s2) <- paste0("s2_ar", counter)
    } 
    stopifnot(length(s2) == 1)
    if (length(phi) < r)
      phi <- c(rep(0.01, r - 1), 0.1)
    if (length(names(phi)) != r)
      names(phi) <- paste0("phi", counter:(counter+r-1))
    b <- c(1, rep(0, r-1))
    C <- diag("0", r, r)
    C[1, ] <- names(phi)
    if (r > 1) {
      for (i in 2:r)
        C[i, i-1] <- "1"
    }
    S <- names(s2)
    param <- c(phi, s2)
  } else if (startsWith("sar", type) || startsWith(type, "sar")) {
    stopifnot(s > 1)
    if (nchar(type) > 3) r <- as.integer(substr(type, 4, nchar(type)))
    else r <- length(phi)
    if (is.null(s2)){
      s2 <- 0.5
      names(s2) <- paste0("s2_sar", counter)
    } 
    stopifnot(length(s2) == 1)
    if (length(phi) < r)
      phi <- c(rep(0.01, r - 1), 0.1)
    if (length(names(phi)) != r)
      names(phi) <- paste0("Phi", counter:(counter+r-1))
    b <- c(1, rep(0, r*s-1))
    C <- diag("0", r*s, r*s)
    if (length(phi) < r)
      phi <- c(0.9, rep(0.1, r-1))
    C[1, (1:r)*s] <- names(phi)
    if (r*s > 1) {
      for (i in 2:(r*s))
        C[i, i-1] <- "1"
    }
    S <- names(s2)
    param <- c(phi, s2)
  } else if (startsWith("cycle", type) || startsWith(type, "cycle")) {
    if (nchar(type) > 5) 
      per <- as.integer(substr(type, 6, nchar(type)))
    stopifnot(per > 2 && per < Inf)
    counter <- as.integer(per)
    if (is.null(s2)) {
      s2 <- c(0.05, 0.05)
      names(s2) <- rep(paste0("s2_cycle", counter), 2)
    }
    stopifnot(length(s2) == 1 || length(s2) == 2)
    if (is.null(names(s2))) 
      names(s2) <- paste0("s2_cycle", counter, ".", 1:length(s2))
    b <- c(1, 0)
    if (phi == 1) {
      C <- matrix(0, 2, 2)
      C[2, 2] <- C[1, 1] <- cos(2*pi/per)
      C[1, 2] <- sin(2*pi/per)
      C[2, 1] <- -C[1, 2] 
      param <- s2
    } else if (phi < 1 && phi > 0) {
      if (is.null(names(phi))) names(phi) <- paste0("phi_cycle", counter)
      if (is.null(names(per))) names(per) <- paste0("per_cycle", counter)
      C <- c(paste0("cos(2*pi/abs(", names(per), "))"),
             paste0("sin(-2*pi/abs(", names(per), "))"),
             paste0("sin(2*pi/abs(", names(per), "))"),
             paste0("cos(2*pi/abs(", names(per), "))"))
      C <- matrix( paste0("abs(", names(phi), ") * ", C), 2, 2)
      param <- c(phi, per, s2)
    } else stop("nonadmissible value for phi")
    S <- names(s2)
  } else if (startsWith("irregular", type)) {
    type <- "irregular"
    if (is.null(s2)) 
      s2 <- c(s2_irreg = 0.5)
    stopifnot(length(s2) == 1)
    if (is.null(names(s2))) 
      names(s2) <- "s2_irreg"
    param <- s2
    b <- NULL
    C <- NULL
    S <- names(s2)
  } else stop("unknown component")
  
  return(uc0(type, param, b, C, S))
  
}

#' Unobserved Components Time Series Models
#'
#' \code{ucm} creates an S3 object representing an UC model:
#'
#' z(t) = b'*x(t) + T(t - j) + S(t - j) + C(t - j) + AR(t - j) + I(t),
#'
#' where z(t) is a time series; x(t) is a set of regressors; T(t - j), S(t - j),
#' C(t - j), AR(t - j) and I(t) are the trend, seasonal, cycle, autoregressive
#' and irregular unobserved components; and i indicates whether the model is 
#' written in lagged form (j = 1) or contemporaneous form (j = 0). 
#' See \code{\link{uc}} and \code{\link{ssm}} for more details.
#'
#' @param z an object of class \code{ts}.
#' @param bc logical. If TRUE, logs are taken.
#' @param uc list of objects of class \code{\link{uc}} or a character with
#'   the name of a predefined model: "llm" local linear model, "lltm" local
#'   linear trend model, "bsm" basic structural model, "bsmd" basic structural
#'   model with dummy seasonality.
#' @param xreg optional design matrix of regressors.
#' @param cform logical. If TRUE, observation equation is given in 
#' contemporaneous form (j = 0); otherwise it is written in lagged form (j = 1).
#' @param fit logical. If TRUE, model is fitted.
#' @param s integer, seasonal period. Optional argument to create a UC model
#'   without providing a time series.
#' @param ... additional parameters for the \code{\link{fit.ssm}} function.   
#'
#' @return An object of class \link{ucm} and class \code{\link{ssm}}.
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
#' ucm1 <- ucm(Nile, uc = "llm")
#' ucm1
#' 
#' @export
#' 
ucm <- function(z, bc = FALSE, uc = NULL, xreg = NULL, cform = TRUE, 
                 fit = TRUE, s = 12, ...) {
  if (missing(z)) {
    z <- NULL
    z.name <- NULL
  } else {
    stopifnot(is.ts(z))
    z.name <- deparse(substitute(z))
  }
  
  if (is.null(uc) || is.character(uc)) {
    if (!is.null(z)) s <- frequency(z)
    if (s < 0) stop("seasonal period 's' must be greater than zero")
    if (is.null(uc)) {
      if (s > 1) uc <- "bsm"
      else ucm <- uc <- "lltm"
    }
    uc <- tolower(uc)
    uc <- switch(uc,
                 "bsm" = uc(c("trend", "seas", "irreg"), s = s),
                 "bsmd" = uc(c("trend","seas","irreg"), s=s, form="dummy"),
                 "llm" = uc(c("trend0", "irreg"), s = s),
                 "lltm" = uc(c("trend", "irreg"), s = s),
                 stop("unknown UC model")
    )
  } else {
    stopifnot(all(sapply(uc, is.uc)))
    i <- sapply(uc, function(x) x$name == "irregular")
    if (sum(i) > 1) stop("only an irregular component can be specified")
    if (!any(i)) {
      uc <- c(uc, list(uc("irregular", s2 = 0)))
    } else {
      if (!i[length(i)]) {
        uc <- c(uc[!i], list(uc[i]))
      }
    }
  }
  k <- length(uc)
  param <- list()
  m <- 0
  for (i in 1:k) {
    param <- c(param, uc[[i]]$param)
    m <- m + uc[[i]]$m
  }
  param <- param[!duplicated(names(param))]
  
  b <- rep(0, m)
  C <- matrix(0, m, m)
  S <- matrix(0, m+1, m+1)
  j1 <- 1
  if (k > 1) {
    for (i in 1:(k - 1)) {
      j2 <- j1 + uc[[i]]$m
      b[j1:(j2-1)] <- uc[[i]]$b
      C[j1:(j2-1), j1:(j2-1)] <- uc[[i]]$C
      S[(j1+1):j2, (j1+1):j2] <- uc[[i]]$S
      j1 <- j2
    }
  }
  S[1, 1] <- uc[[k]]$S
  
  ucm1 <- list(z = z, z.name = z.name, bc = bc, b = b, C = C, S = S, uc = uc,
                k = k, m = m, xreg = xreg, a = NULL, param = param, cform = cform)
  class(ucm1) <- c("ucm", "ssm")
  
  if (fit && !is.null(z)) 
    ucm1 <- fit.ssm(ucm1, method = "BFGS", ...)
  
  return(ucm1)
}

#' Time Invariant State Space Model
#'
#' \code{ssm} creates an S3 object representing a time-invariant state space
#' model:
#'
#' z(t) = b'x(t-j) + u(t) (observation equation), 
#' x(t) = Cx(t-1) + v(t) (state equation),
#' j = 0 for contemporaneous form or j = 1 for lagged form. 
#' Note: the lagged form (j=1) is equivalent to the future form 
#' x(t+1) = Cx(t) + v(t+1).
#'
#' @param z an object of class \code{ts}.
#' @param b vector of coefficients (m-dimensional).
#' @param C matrix of coefficients (m x m).
#' @param S covariance matrix of the error vector [u(t), v(t)], (m+1 x m+1).
#' @param xreg design matrix of regressors.
#' @param bc logical. If TRUE logs are taken.
#' @param cform logical. If TRUE state equation is given in contemporaneous
#'   form, otherwise it is written in lagged form.
#' @param tol tolerance to check if a value is zero or one.
#'
#' @return An object of class \code{ssm} containing:
#' \item{z}{the input time series}
#' \item{b}{observation coefficients}
#' \item{C}{state transition matrix}
#' \item{S}{error covariance matrix}
#' \item{xreg}{regressor matrix (if provided)}
#' \item{a}{regression coefficients for xreg (if computed)}
#' \item{z.name}{name of the input series}
#' \item{bc}{Box-Cox transformation indicator}
#' \item{m}{number of state variables}
#' \item{cform}{form indicator (contemporaneous vs lagged)}
#' \item{is.adm}{admissibility flag}
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
#' ssm1 <- ssm(Nile, b, C, S = diag(c(irr = 1, lvl = 0.5)) )
#' ssm1
#' 
#' @export
#' 
ssm <- function(z, b, C, S, xreg = NULL, bc = FALSE, cform = TRUE, tol = 1e-5) {
  
  if (!is.matrix(C)) {
    if (length(C) == 1) C <- as.matrix(C)
    else stop("C argument must be a matrix")
  }
  stopifnot("C matrix must be square" = nrow(C) == ncol(C))
  stopifnot("S argument must be a matrix" = is.matrix(S))
  stopifnot("S matrix must be square" = nrow(S) == ncol(S))
  if (nrow(S) == nrow(C)) {
    S <- rbind(rep(0, ncol(C)+1), cbind(rep(0, nrow(C)), S))
  }
  
  if (is.matrix(b)) b <- as.vector(b)
  stopifnot(length(b) == nrow(C) || (length(b)+1) == nrow(S))
  
  if (missing(z)) z <- NULL
  if (!is.null(z)) {
    z.name <- deparse(substitute(z))
    if (!is.null(xreg)) {
      if (is.matrix(xreg))
        stopifnot(nrow(xreg) == length(z))
      else if (is.vector(xreg)||is.ts(xreg)) {
        stopifnot(length(xreg) == length(z))
        xreg <- matrix(xreg, ncol = 1)
      } else stop("wrong xreg argument")
    }
  } else z.name <- NULL

  if (!is.null(xreg) && !is.null(z)) {
    a <- tryCatch({
      as.numeric(solve(t(xreg) %*% xreg, t(xreg) %*% z))
    }, error = function(e) {
      stop("Singular design matrix in regression")
    })    
  } else a <- NULL
  
  if(any(abs(S - t(S)) > tol))
    stop("non-symmetric S matrix")
  if(min(eigen(S, symmetric = TRUE, only.values = TRUE)$values) < -tol)
    stop("negative definite S matrix")
  if(any(abs(eigen(C)$values) > 1 + tol)) 
    stop("non-stationary C matrix")
  
  mdl <- list(z = z, b = b, C = C, S = S, xreg = xreg, a = a, z.name = z.name, 
              bc = bc, m = length(b), cform = cform, is.adm = TRUE)
  
  class(mdl) <- "ssm"
  return(mdl)
}

#' @rdname fit
#' @param z an object of class \code{\link{ts}}, univariate time series.
#' @param updateSSM  user function to update the parameters of the SS model. 
#' The function must take a model object and a parameter vector as inputs and
#' return an updated model object.
#' @param param a numeric vector of named parameters passed to the 
#' \code{updateSSM} function.
#' @param show.iter logical. If TRUE, displays parameter estimates at each iteration.
#' @param tol numeric. Tolerance to check if a root is close to one.
#' @param method a character string specifying the optimization method. See the 
#' \code{\link{optim}} function documentation for details. Common options include 
#' "BFGS" and "Nelder-Mead".
#' @param ... additional arguments for the \code{\link{optim}} function.
#'
#' @return An object of class "ssm" with the estimated parameters.
#'
#' @examples
#' # Predefined local level model
#' ucm1 <- ucm(Nile, uc = "llm", fit = FALSE)
#' ucm1 <- fit(ucm1)
#' ucm1
#' 
#' # User defined local level model
#' ssm1 <- ssm(Nile, b = 1, C = 1, S = diag(c(1, 0.5)) )
#' param <- c(irr = var(Nile), lvl = var(diff(Nile)))
#' updateSSM <- function(mdl, param) {
#' mdl$S[1,1] <- param[1]
#' mdl$S[2,2] <- param[2]
#' mdl
#' }
#' fit(ssm1, updateSSM = updateSSM, param = param)
#' 
#' @export
fit.ssm <- function(mdl, z = NULL, updateSSM, param, show.iter = FALSE, 
                    tol = 1e-4, method = "BFGS", ...) {

  fit_env <- new.env(parent = emptyenv())
  fit_env$res <- 0
  fit_env$scale <- 1
  fit_env$iter <- 0
  
  # log-likelihood for SS model
  .logLik <- function(par, type = 0, ...) {
    param1 <- param
    param1[free] <- par
    param1[is_variance & free] <- exp(param1[is_variance & free])
    mdl1 <- updateSSM(mdl, param1)
    if (mdl1$is.adm) {
      if (p > 0) uca <- .ssm2ucarima(mdl1, tol = tol)
      if (!is.null(mdl1$xreg)) {
        if (type == 0)
          ll <- llucaC(w - X %*% param1[1:kX], uca$phi, uca$MA, mdl1$S,
                       fit_env$scale, type)
        else {
          fit_env$res <- w - X %*% param1[1:kX] # w is changed by llucaC function
          ll <- llucaC(fit_env$res, uca$phi, uca$MA, mdl1$S, fit_env$scale, type)
        }
      } else {
        ll <- llucaC(w, uca$phi, uca$MA, mdl1$S, fit_env$scale, type)
      }
      if (show.iter) {
        if ( fit_env$iter == 0) {
          cat(format("iter", width = 4, justify = "right"))
          cat(format("logLik", width = 10, justify = "right"))
          cat(format("scale", width = 10, justify = "right"))
          for (nm in names(param1[is_variance]))
            cat(format(nm, width = 10, justify = "right"))
          cat("\n")
        }
        if (fit_env$iter %% 10 == 0) {
          cat(format(fit_env$iter, width = 4, justify = "right"))
          cat(format(round(ll, 6), width = 10, justify = "right"))
          cat(format(round(fit_env$scale, 6), width = 10, justify = "right"))
          for (par1 in param1[is_variance])
            cat(format(round(par1, 6), width = 10, justify = "right"))
          cat("\n")
        }
      }
      fit_env$iter <- fit_env$iter + 1
    } else ll <- fit_env$ll0
    return(-ll)
  }
  
  # Function to calculate two types of residuals
  .resuca <- function(par, type = 1, ...) {
    ll <- .logLik(par, type, ...)
    return(fit_env$res)
  }
  
  # Indentify variance parameters
  .varianceIndicator <- function(par) {
    npar <- length(par)
    is_variance <- rep(FALSE, npar)
    S <- diag(mdl$S)
    par1 <- par
    for (i in 1:npar) {
      par1[i] <- par1[i] + 100 # Check if changing a parameter affects S
      mdl1 <- updateSSM(mdl, par1[1:npar])
      if (any( (diag(mdl1$S) - S) > 10 )) 
        is_variance[i] <- TRUE
      par1[i] <- par[i]
    }
    is_variance
  }
  
  # Check if z is provided
  if (is.null(z)) {
    if(!is.ts(mdl$z))
      stop("argument z required")
  } else if (is.ts(z)) {
    mdl$z <- z
    mdl$z.name <- deparse(substitute(z))
  } else stop("z must be a time series")
  if (any(is.na(mdl$z))) stop("z has NA values")
  if (mdl$bc && min(mdl$z) <= 0) {
    warning("z has negative values: bc argument set to FALSE")
    mdl$bc <- FALSE
  }
  
  # Check updateSSM function and param argument 
  if (inherits(mdl, "ucm")) {
    if (missing(updateSSM)) updateSSM <- .update_ucm
    if (missing(param)) param <- mdl$param
  } 
  stopifnot("updateSSM function is required" = !missing(updateSSM),
            "param argument is requiered" = !missing(param))
  if (is.list(param)) param <- unlist(param)
  stopifnot("updateSSM must be a function" = is.function(updateSSM),
            "param must be numeric" = is.numeric(param))
  free <- rep(TRUE, length(param)) # free parameters
  # Fix zero variances and set scale
  is_variance <- .varianceIndicator(param)
  
  if (!any(is_variance)) 
    stop("No error variance to estimate.")
  for (i in 1:3) {
    fit_env$scale <- max(param[is_variance])
    param[is_variance] <- param[is_variance] / fit_env$scale
    param[is_variance & param < 0] <- 0
    free[is_variance] <- (param[is_variance] > .Machine$double.eps)
    idx <- which(is_variance & free & param[is_variance] == 1)
    if (length(idx) > 0) free[idx[1]] <- FALSE
    else stop("No free error variance to estimate.")
    if (i == 1) {
      uca <- .ssm2ucarima(mdl, tol = tol)
      p <- length(uca$phi) - 1
      w <- diffC(mdl$z, uca$nabla, mdl$bc)
      n <- length(w)
      if (!is.null(mdl$xreg)) {
        kX <- ncol(mdl$xreg)
        X <- sapply(1:kX, function(i) diffC(mdl$xreg[, i], uca$nabla, FALSE))
        X <- X[1:n, ]
        if (kX == 1) X <- as.matrix(X)
        a <- tryCatch({
          as.numeric(solve(t(X) %*% X, t(X) %*% w))
        }, error = function(e) {
          stop("Singular design matrix in regression")
        })    
        if (length(colnames(mdl$xreg)) == kX)
          names(a) <- colnames(mdl$xreg)
        else
          names(a) <- paste0("a", 1:kX)
        param <- c(a, param)
        free <- c(rep(TRUE, kX), free)
        is_variance <- c(rep(FALSE, kX), is_variance)
      }
    }
    mdl <- updateSSM(mdl, param)
    par <- param
    par[is_variance & free] <- log(par[is_variance & free])
    par <- par[free]
    mdl$optim <- optim(par, .logLik, hessian = TRUE, method = method, ...)
    param[free] <- mdl$optim$par
    param[is_variance & free] <- exp(param[is_variance & free])
    if (any(param[is_variance & free] > 1.01) && i < 3) {
      param[is_variance] <- param[is_variance]*fit_env$scale
    } else break
  }

  mdl$fit <- list()
  mdl$fit$n <- length(w)
  mdl$fit$param <- param
  mdl$fit$is_variance <- is_variance
  mdl$fit$free <- free
  mdl$fit$scale <- fit_env$scale*1
  show.iter <- FALSE
  res <- .resuca(mdl$optim$par)
  J <- numDeriv::jacobian(.resuca, mdl$optim$par)
  mdl$fit$g <- t(J) %*% res
  mdl$fit$varb <- tryCatch ({ 
    MASS::ginv(t(J) %*% J)*(sum(res^2)/length(res))
  }, error = function(e) {
    message("Error inverting matrix: ", e$message)
    return(NULL)
  })
  param[is_variance] <- param[is_variance]*mdl$fit$scale
  mdl <- updateSSM(mdl, param)
  mdl$fit$param <- param
  mdl$fit$sigmas <- param[is_variance]
  mdl$fit$sigmas <- cbind(Variance = mdl$fit$sigmas, 
                        Ratio = mdl$fit$sigmas/mdl$fit$scale)
  rownames(mdl$fit$sigmas) <- names(param)[is_variance]
  
  if (!is.null(mdl$xreg))
    mdl$a <- param[1:kX]
  mdl$fit$aic <- (2.0*mdl$optim$value+2*length(mdl$optim$par))/mdl$fit$n
  um1 <- .ssm2um(mdl)
  mdl$fit$sig2 <- um1$sig2

  return(mdl)
}

#' Print method for unobserved components
#'
#' @param x an object of class \code{\link{uc}}.
#' @param ... currently unused.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @export
print.uc <- function(x, ...) {
  print(x$name)
  if (x$m > 0) {
    A <- cbind(b = x$b, C = x$C, S = x$S)
    colnames(A) <- c("b", "C", rep("", x$m-1), "S", rep("", x$m-1))
  } else A <- x$S
  print(A)
  invisible(NULL)
}

#' Print method for \code{\link{ssm}} objects
#' 
#' @param x an object of class \code{\link{ssm}}.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return Invisibly returns \code{NULL}.
#' 
#' @export
print.ssm <- function(x, ...) {
  if (!is.null(x$optim)) {
    if (any(!x$fit$is_variance)) {
      if (!is.null(x$fit$varb)) se <- sqrt(diag(x$fit$varb))
      else se <- rep(NA, length(x$optim$par))
      names(se) <- names(x$optim$par)
      nms <- names(x$fit$param)[!x$fit$is_variance]
      b <- cbind(x$optim$par[nms], se[nms])
      colnames(b) <- c("Estimate", "Std. Error")
      print(b)
      cat("\n")
    }
    print(x$fit$sigmas)
    cat("\nlog likelihood: ", -x$optim$value)
    cat("\nResidual standard error: ", sqrt(x$fit$sig2))
    cat("\nAIC:", x$fit$aic)
  } else {
    cat("State Space Model (Time Invariant)\n")
    cat("Observation Equation: ")
    if (x$cform) cat("z(t) = b'x(t) + u(t)\n")
    else cat("z(t) = b'x(t-1) + u(t)\n")  
    cat("State Equation: x(t) = Cx(t-1) + v(t)\n")
    cat("z: ", if (is.null(x$z.name)) "NULL" else x$z.name, "\n")
    cat("bc: ", x$bc, "\n")    
    if (x$m > 0) {
      A <- cbind(b = x$b, C = x$C)
      A <- cbind(rbind(rep(NA, ncol(A)), A), S = x$S)
      colnames(A) <- c("b", "C", rep("", x$m-1), "S", rep("", x$m))
      rownames(A) <- c("u", paste0("x", 1:x$m))
    } else A <- x$S
    print(A)
  }
  return(invisible(NULL))
}

#' Summary of fitted state space model
#'
#' \code{summary.ssm} generates a detailed summary of the results of the
#' estimation of a \code{\link{ssm}} object.
#'
#' @param object An object of class \code{\link{ssm}}.
#' @param ... additional arguments.
#' @return An object of class \code{summary.ssm}.
#' @export
summary.ssm <- function(object, ...) {
  stopifnot(is.ssm(object))
  if (is.null(object$optim)) {
    stop("Model must be fitted before computing summary")
  }  
  um1 <- .ssm2um(object)
  um1$bc <- FALSE
  if (object$bc) z <- log(object$z)
  else z <- object$z 
  N <- length(z)
  start <- start(object$z)
  end <- end(object$z)
  s <- frequency(object$z)
  
  if (!is.null(object$xreg)) 
    z <- z - object$xreg[1:N, ] %*% object$a
  w <- diffC(z, um1$nabla, FALSE)
  n <- length(w)
  res <- residuals.um(um1, z)
  mean.resid <- mean(res)
  rss <- var(res)*(length(res)-1)
  sd.resid <- sd(res)
  z.mean.resid <- mean.resid*sqrt(length(res))/sd.resid
  p.mean.resid <- pnorm(abs(z.mean.resid), lower.tail = F)*2
  
  ll <- -object$optim$value
  aic <- (-2.0*ll+2*length(object$optim$par))/n
  bic <- (-2*ll+log(n)*length(object$optim$par))/n
  
  Q1 <- Box.test(res, lag = um1$p+um1$q+1, type = "Ljung-Box", fitdf = um1$p+um1$q)
  Q2 <- Box.test(res, lag = as.integer(n/4)+um1$p+um1$q, type = "Ljung-Box",
                 fitdf = um1$p+um1$q)
  Q <- c(Q1 = Q1, Q2 = Q2)
  h.stat <- bartlett.test(res, g = cut(1:length(res), 3))  
  
  if (any(!object$fit$is_variance)) {
    b <- object$optim$par[!object$fit$is_variance]
    se <- sqrt(diag(object$fit$varb))[!object$fit$is_variance]
    z.ratio <- b/se
    p <- pnorm(abs(z.ratio), lower.tail = FALSE)*2
    X <- cbind(b, se, z.ratio, p)
    colnames(X) <- c("Estimate", "Std. Error", "z Value", "Pr(>|z|)")
    rownames(X) <- names(object$optim$par[!object$fit$is_variance])
  } else X <- NULL
  
  var <- unlist(object$param)
  var <- cbind(Variances = var, Ratio = var/object$fit$scale)
  rownames(var) <- names(object$fit$param[object$fit$is_variance])
  
  out <- list(z = object$z.name, N = N, n = n, logLik = ll, aic = aic, bic = bic,
              table = X, var = var, resid = res, sig2 = um1$sig2,
              mean.resid = mean.resid, sd.resid = sd.resid, 
              z.mean.resid = z.mean.resid, p.mean.resid = p.mean.resid,
              Q = Q, h.stat = h.stat, start = start, end = end, s = s)
  class(out) <- "summary.ssm"
  out
}

#' Print summary of fitted state space model
#'
#' \code{print.summary.ssm} prints a detailed summary of the results of the
#' estimation of a \code{\link{ssm}} object.
#'
#' @param x an object of class \code{summary.ssm}.
#' @param stats logical. If \code{TRUE}, additional diagnostic statistics are 
#' printed.
#' @param digits integer. Number of significant digits to display for numeric 
#' values. 
#' @param ... additional arguments.
#' @return Invisible \code{NULL}.
#' @export
print.summary.ssm <- function(x, stats = TRUE, digits = max(3L, getOption("digits") - 3L),  ...) {
  stopifnot(inherits(x, "summary.ssm")) 
  
  cat("\nTime series: ", x$z.name, "\n")

  if (!is.null(x$table)) {  
    cat("\nCoefficients:\n")
    printCoefmat(x$table, digits = digits, na.print = "NA")
    cat("\n")
  }

  cat("Error variances:\n")
  print(x$var, digits = digits)

  if (stats) {
    w.left <- 20
    w.right <- 10
    cat("\nStatistics:\n")
    cat(format("Total nobs", width = w.left, justify = "left"),
        format(x$N, digits = NULL, width = w.right, justify = "right"), 
        format("Effective nobs", width = w.left, justify = "left"),
        format(x$n, digits = NULL, width = w.right, justify = "right"), "\n")
    
    cat(format("log likelihood", width = w.left, justify = "left"),
        format(x$logLik, digits = digits, width = w.right, justify = "right"),
        format("Error variance", width = w.left, justify = "left"),
        format(x$sig2, digits = digits, width = w.right, justify = "right"), "\n")
  
    cat(format("Mean of residuals", width = w.left, justify = "left"),
        format(x$mean.resid, digits = digits, width = w.right, justify = "right"),
        format("SD of residuals", width = w.left, justify = "left"),
        format(x$sd.resid, digits = digits, width = w.right, 
               justify = "right"), "\n")
    
    cat(format("z-test for residuals", width = w.left, justify = "left"),
        format(x$z.mean.resid, digits = digits, width = w.right,
               justify = "right"),
        format("p-value", width = w.left, justify = "left"),
        format(x$p.mean.resid, digits = digits, width = w.right,
               justify = "right"), "\n")
    
    
    label <- paste0("Ljung-Box Q(", x$Q$Q1.parameter, ") st.")
    cat(format(label, width = w.left, justify = "left"),
        format(x$Q$Q1.statistic, digits = digits, width = w.right,
               justify = "right"),
        format("p-value", width = w.left, justify = "left"),
        format(x$Q$Q1.p.value, digits = digits, width = w.right,
               justify = "right"), "\n")
    
    label <- paste0("Ljung-Box Q(", x$Q$Q2.parameter, ") st.")
    cat(format(label, width = w.left, justify = "left"),
        format(x$Q$Q2.statistic, digits = digits, width = w.right,
               justify = "right"),
        format("p-value", width = w.left, justify = "left"),
        format(x$Q$Q2.p.value, digits = digits, width = w.right,
               justify = "right"), "\n")
    
    cat(format("Barlett H(3) stat.", width = w.left, justify = "left"),
        format(x$h.stat$statistic, digits = digits, width = w.right,
               justify = "right"),
        format("p-value", width = w.left, justify = "left"),
        format(x$h.stat$p.value, digits = digits, width = w.right, 
               justify = "right"), "\n")  
    
    cat(format("AIC", width = w.left, justify = "left"),
        format(x$aic, digits = digits, width = w.right, justify = "right"),
        format("BIC", width = w.left, justify = "left"),
        format(x$bic, digits = digits, width = w.right, justify = "right"), "\n")
    
  }
  
  invisible(NULL)
}

#' Residuals of fitted state space models
#'
#' \code{residuals.ssm} generates the residuals of a fitted \code{\link{ssm}}
#' object.
#'
#' @param object an object of class \code{\link{ssm}}.
#' @param method character. Either "exact" or "conditional" residuals.
#' @param ... additional arguments.
#' @return A \code{ts} object containing the residuals of the SSM.
#' @details These residuals are calculated by first converting the
#'   \code{\link{ssm}} object to a \code{\link{um}} object and then using the
#'   \code{\link{residuals.um}} function.
#' @examples
#' 
#' # Local level model
#' b <- 1
#' C <- as.matrix(1)
#' ssm1 <- ssm(Nile, b, C, S = diag(c(irr = 15127.7, lvl = 1453.2)))
#' u <-residuals(ssm1)
#' 
#' @export
residuals.ssm <- function(object, method = c("exact", "cond"), ...) {
  stopifnot(is.ssm(object))
  stopifnot(!is.null(object$z))
  method <- match.arg(method)
  um1 <- .ssm2um(object)
  um1$bc <- FALSE
  if (object$bc) z <- log(object$z)
  else z <- object$z 
  N <- length(z)
  if (!is.null(object$xreg)) 
    z <- z - object$xreg[1:N, ] %*% object$a
  res <- residuals.um(um1, z, method = method, ...)
  res <- ts(res, end = end(object$z), frequency = frequency(object$z))
  res
}

#' Log-likelihood of a SS model
#'
#' \code{logLik.ssm} computes the exact or conditional log-likelihood of a state
#' space model.
#'
#' @param object an object of class \code{\link{ssm}}.
#' @param method character. Either "exact" or "conditional" maximum likelihood.
#' @param ... additional parameters.
#'
#' @return The log-likelihood value.
#' 
#' @examples
#' 
#' # Local level model
#' b <- 1
#' C <- as.matrix(1)
#' ssm1 <- ssm(Nile, b, C, S = diag(c(irr = 15127.7, lvl = 1453.2)))
#' logLik(ssm1)
#'  
#' @export
logLik.ssm <-function(object, method = c("exact", "cond"), ...) {
  stopifnot(is.ssm(object))
  stopifnot(!is.null(object$z))
  method <- match.arg(method)
  uca <- .ssm2ucarima(object)
  if (object$bc) z <- log(object$z)
  else z <- object$z 
  N <- length(z)
  if (!is.null(object$xreg)) 
    z <- z - object$xreg[1:N, ] %*% object$a
  w <- diffC(z, uca$nabla, FALSE)
  s2 <- 1
  ll <- llrfC(w, uca$nabla, object$b, object$C, object$S, s2, object$cform)
  ll
}

#' AIC for fitted state space models
#'
#' @param object an object of class \code{\link{ssm}}.
#' @param k numeric, penalty per parameter (default is 2).
#' @param ... additional arguments.
#'
#' @return The AIC value, or \code{NULL} if model is not fitted.
#'
#' @export
AIC.ssm <- function(object, k = 2, ...) {
  stopifnot(is.ssm(object))
  if (is.null(object$optim)) return(NULL)
  ll <- -object$optim$value
  b <- object$optim$par
  aic <- (-2.0*ll+k*length(b))/object$fit$n
}

#' @rdname diagchk
#' @param mdl an object of class \code{\link{ssm}}.
#' @param lag.max integer; maximum number of lags for ACF/PACF.
#' @param lags.at numeric vector; specific lags in ACF/PACF plots.
#' @param freq.at numeric vector; specific frequencies in (cum) periodogram plot.
#' @param std logical; if TRUE standardized residuals are used.
#' @param ... additional arguments passed to \code{\link{residuals.ssm}}.
#'    
#' @examples
#' 
#' # Local level model
#' b <- 1
#' C <- as.matrix(1)
#' ssm1 <- ssm(Nile, b, C, S = diag(c(irr = 15127.7, lvl = 1453.2)))
#' diagchk(ssm1)
#' @export
diagchk.ssm <- function(mdl, lag.max = NULL, lags.at = NULL, freq.at = NULL,
                        std = TRUE, ...) {
  stopifnot(is.ssm(mdl))
  u <- residuals.ssm(mdl, ...)
  ide(u, graphs = c("plot", "hist", "acf", "pacf", "cpgram"), ylab = "u",
      lag.max = lag.max, lags.at = lags.at, freq.at = freq.at,
      std = std, ...)
}

#' Predict method for state space models
#'
#' \code{predict.ssm} generates forecasts from fitted state space models by 
#' converting to univariate ARIMA or transfer function models and using their 
#' prediction methods.
#'
#' @param object an object of class \code{\link{ssm}}.
#' @param ... additional arguments passed to \code{\link{predict.um}} (for models
#'   without regressors) or \code{\link{predict.tfm}} (for models with regressors).
#'   See their documentation for available options such as forecast horizon,
#'   confidence levels, etc.
#'
#' @return An object of class \code{\link{predict.um}} or \code{\link{predict.tfm}}.
#'
#' @seealso \code{\link{predict.um}}, \code{\link{predict.tfm}}
#'
#' @export
predict.ssm <- function(object, ...) {
  stopifnot(is.ssm(object))
  um1 <- .ssm2um(object)
  if (is.null(object$xreg))
    predict.um(um1, z = object$z, ...)
  else {
    tfm1 <- tfm(xreg = object$xreg, noise = um1, fit = FALSE, new.name = FALSE)
    tfm1$param[1:ncol(object$xreg)] <- object$a
    predict.tfm(tfm1, object$z, ...)
  }
}

#' @rdname decomp
#' @param mdl an object of class \code{\link{ssm}}.
#' @param tol numeric tolerance for classifying eigenvalues.
#' @param ... additional arguments passed to \code{\link{ks}}.
#' @export
decomp.ssm <- function(mdl, tol = 1e-5, ...) {
  stopifnot(is.ssm(mdl))
  if (!mdl$cform)
    mdl <- .switchSSMform(mdl)
  ks1 <- ks(mdl, ...)
  nms <- c("series")
  X <- cbind(ks1$ztilde)
  m <- c()
  if (inherits(mdl, "ucm")) {
    for (i in 1:(mdl$k-1)) 
      m <- c(m, mdl$uc[[i]]$m)
    nms <- c(nms, names(mdl$uc))
  } else {
    r <- eigen(mdl$C)$values
    l <- sum(abs(1-Mod(r)) > tol)
    if (l > 0) {
      m <- c(m, l)
      nms <- c(nms, "trans")
    }
    l <- sum(abs(1 - Re(r)) < tol & abs(Im(r)) < tol)
    if (l > 0) {
      m <- c(m, l)
      nms <- c(nms, "trend")
    }
    l <- mdl$m - sum(m) 
    if (l > 0) {
      m <- c(m, l)
      nms <- c(nms, "seas")
    }
    nms <- c(nms, "irreg")
  }
  
  l1 <- 1
  for (l2 in m) {
    X <- cbind(X, ks1$X[, l1:(l1+l2-1), drop = FALSE] %*% mdl$b[l1:(l1+l2-1)])
    l1 <- l1 + l2
  }
  irreg <- X[, 1]
  for (i in 2:ncol(X))
    irreg <- irreg - X[, i]
  X <- cbind(X, irreg)
  colnames(X) <- nms
  return(X)
}

#' @rdname calendar
#' @param mdl an object of class \code{\link{ssm}}.
#' @param ... additional arguments passed to \code{\link{fit.ssm}}.
#' @export
calendar.ssm <- function(mdl, form = c("dif", "td", "td7", "td6", "wd"),
           ref = 0, lom = TRUE, lpyear = TRUE, easter = FALSE, len = 4, 
           easter.mon = FALSE, n.ahead = 0, p.value = 1, envir = NULL, ...) {
  stopifnot(is.ssm(mdl))
  stopifnot(is.ts(mdl$z) && frequency(mdl$z) == 12)
  xreg <- CalendarVar(mdl$z, form, ref, lom, lpyear, easter, len, 
                      easter.mon, n.ahead)
  if (is.null(mdl$xreg)) mdl$xreg <- xreg
  else mdl$xreg <- cbind(xreg, mdl$xreg)
  fit.ssm(mdl, ...)
}

#' @rdname outliers
#' @param mdl an object of class \code{\link{ssm}}.
#' @param ... additional arguments passed to \code{\link{fit.ssm}}.
#' @export
outliers.ssm <- function(mdl, types = c("AO", "LS", "TC"), dates = NULL, c = 3, 
                         calendar = FALSE, easter = FALSE, 
                         resid = c("exact", "cond"), n.ahead = 0, p.value = 1, 
                         ...) {
  stopifnot(is.ssm(mdl))
  um1 <- .ssm2um(mdl)
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

#' Extract coefficients from UCM objects
#'
#' @param object an object of class \code{\link{ucm}}.
#' @param ... currently unused.
#'
#' @return A named numeric vector of model coefficients.
#'
#' @export
coef.ucm <- function(object, ...) {
  coef1 <- as.numeric(object$param)
  names(coef1) <- names(object$param)
  coef1 <- c(object$a, coef1)
  return(coef1)
}

#' @rdname autocov
#' @param mdl an object of class \code{\link{ssm}}.
#' @param arma logical. If TRUE, the autocovariances for the stationary ARMA
#'   model of the reduced form are computed. Otherwise, the autocovariances are
#'   only computed for the MA part.
#' @param varphi logical. If TRUE, the varphi polynomial of the reduced form is
#'   also returned.
#' @param tol tolerance to check if a root is close to one.    
#' 
#' @examples
#' # Local level model
#' ssm1 <- ssm(b = 1, C = 1, S = diag(c(irr = 0.8, lvl = 0.04)))
#' autocov(ssm1)
#'
#' @export
#' 
autocov.ssm <- function(mdl, lag.max = NULL, arma = TRUE, varphi = FALSE,
                        tol = 1e-4, ...) {
  stopifnot(is.ssm(mdl))
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

#' @rdname noise
#' @export
noise.ssm <- function(mdl, diff = TRUE, exp = FALSE, ...) {
  stopifnot(is.ssm(mdl))
  if (mdl$bc) z <- log(mdl$z)
  else z <- mdl$z
  N <- length(z)
  if (!is.null(mdl$xreg)) 
    z <- z - as.matrix(mdl$xreg[1:N, ]) %*% mdl$a
  if (diff) {
    uca <- .ssm2ucarima(mdl)
    z <- diffC(z, uca$nabla, FALSE)
  } else if(exp && mdl$bc) {
    z <- exp(z)
  }
  z <- ts(z, end = end(mdl$z), frequency = frequency(z))
}


#' Initialization of Kalman filter
#'
#' \code{init_kf} computes the starting values x0 and P0 by generalized least
#' squares using the first n observations.
#' @param mdl an object of class \code{ssm}.
#' @param z optional time series if it differs from model series.
#' @param n integer, number of observations used to estimate the initial
#'   conditions. If n < d (dimension of state vector), it defaults to length(z).
#'   
#' @return A list with components:
#' \item{x0}{initial state vector estimate}
#' \item{P0}{covariance matrix of the initial state estimate}
#' 
#' @export
init_kf <- function(mdl, z = NULL, n = 0) {
  stopifnot(is.ssm(mdl))
  cform <- mdl$cform
  if (cform)
    mdl <- .switchSSMform(mdl)
  d <- length(mdl$b)
  if (is.null(z))
    z <- mdl$z
  if (n < d || n > length(z))
    n <- length(z)
  z <- z[1:n]
  if (any(is.na(z[1:n]))) 
    stop("Initial observations contain NA values")  
  if (mdl$bc && min(z) > 0) z <- log(z)
  if (!is.null(mdl$xreg))
    z <- z - mdl$xreg[1:n,] %*% mdl$a 
  x <- mdl$b
  X <- x
  if (n > 1) {
    for (i in 2:n) {
      x <- x %*% mdl$C
      X <- rbind(X, x)
    }
  } else X <- as.matrix(X)
  x0 <- matrix(0, nrow = d, ncol = 1)
  P0 <- matrix(0, nrow = d, ncol = d)
  v <- matrix(0, nrow = n, ncol = 1)
  s2v <- matrix(0, nrow = n, ncol = 1)
  if(!kf0C(z, mdl$b, mdl$C, mdl$S, x0, P0, v, s2v))
    stop("Kalman filter algorithm failed")
  z <- v/sqrt(s2v)
  X <- sapply(1:d, function(j) { 
    kf0C(X[,j], mdl$b, mdl$C, mdl$S, x0, P0, v, s2v)
    v/sqrt(s2v)
  })
  XtX <- t(X) %*% X
  iXX <- tryCatch({
    solve(XtX)
  }, error = function(e) {
    if (det(XtX) < .Machine$double.eps) {
      warning("Design matrix is singular, using pseudo-inverse")
      MASS::ginv(XtX)
    } else {
      stop("Error inverting matrix: ", e$message)
    }
  })
  Xy <- t(X) %*% z
  b <- iXX %*% Xy
  s2 <- sum((z - X%*%b)^2)/(n-d)
  vb <- s2*iXX
  if (n == d) vb <- diag(10^4, n, n)
  list(x0 = b, P0 = vb)
}

#' Kalman filter for SS models
#'
#' \code{kf} computes the innovations and the conditional states with the Kalman
#' filter algorithm.
#'
#' @param mdl an object of class \code{ssm}.
#' @param z time series to be filtered when it differs from the model series.
#' @param x0 initial state vector.
#' @param P0 covariance matrix of x1.
#' @param filtered logical. If TRUE, the filtered states x_\{t|t\} and their
#'   covariance matrices P_\{t|t\} are returned. Otherwise, the forecasted states 
#'   x_\{t|t-1\} and thier covariance matrices P_\{t|t-1\} are returned.
#' @param ... additional arguments.
#' @return A list with the innovations, the conditional states and their
#'   covariance matrices.
#' @export
kf <- function(mdl, z = NULL, x0 = NULL, P0 = NULL, filtered = FALSE, ...) {
  stopifnot(is.ssm(mdl))
  cform <- mdl$cform
  if (cform) 
    mdl <- .switchSSMform(mdl)
  d <- length(mdl$b)
  if (is.null(x0)||is.null(P0)) {
    ic <- init_kf(mdl, z)
    if (is.null(x0)) 
      x0 <- ic[[1]]
    if (is.null(P0)) 
      P0 <- ic[[2]]
  }
  if (is.null(z))
    z <- mdl$z
  end <- end(z)
  s <- frequency(z)
  n <- length(z)
  if (mdl$bc) z <- log(z)
  if (!is.null(mdl$xreg))
    z <- z - (mdl$xreg %*% mdl$a) 
  X <- matrix(0, nrow = n, ncol = d)
  PX <- matrix(0, nrow = d*n, ncol = d)
  v <- matrix(0, nrow = n, ncol = 1)
  s2v <- matrix(0, nrow = n, ncol = 1)
  if(!kfC(z, mdl$b, mdl$C, mdl$S, x0, P0, v, s2v, X, PX, filtered))
    stop("Kalman filter failed")
  logdet <- mean(log(s2v))
  s2 <- sum(v^2/s2v)/n
  ll <- -0.5*n*(1 + log(2*pi) + logdet + log(s2))
  v <- ts(v, end = end, frequency = s)
  X <- ts(X, end = end, frequency = s)
  list(logLik = ll, se = sqrt(s2), X = X, PX = PX, u = v, s2u = s2v)
}

#' Kalman smoother for SS models
#'
#' \code{ks} computes smoothed states and their covariance matrices.
#'
#' @param mdl an object of class \code{ssm}.
#' @param x0 initial state vector.
#' @param P0 covariance matrix of x0.
#' @export
ks <- function(mdl, x0 = NULL, P0 = NULL) {
  stopifnot(is.ssm(mdl))
  cform <- mdl$cform
  if (cform) 
    mdl <- .switchSSMform(mdl)
  if (is.null(x0)||is.null(P0)) {
    ic <- init_kf(mdl)
    if (is.null(x0)) 
      x0 <- ic[[1]]
    if (is.null(P0)) 
      P0 <- ic[[2]]
  }
  d <- length(mdl$b)
  z <- mdl$z
  n <- length(z)
  if (mdl$bc) z <- log(z)
  if (!is.null(mdl$xreg))
    z <- z - mdl$xreg %*% mdl$a 
  X <- matrix(0, nrow = n, ncol = d)
  PX <- matrix(0, nrow = d*n, ncol = d)
  ksC(z, mdl$b, mdl$C, mdl$S, x0, P0, X, PX)
  X <- ts(X, end = end(mdl$z), frequency = frequency(mdl$z))
  list(ztilde = z, X = X, PX = PX, x0 = x0, P0 = P0)
}

#
# Internal/Hidden functions for ssm and related objects
#

#' Leverrier-Faddeev Algorithm
#'
#' \code{.LeverrierFaddeev} implements the Leverrier-Faddeev algorithm to
#' calculate the coefficients of the characteristic polynomial of the transition
#' matrix (\code{C}) in a state-space model (SSM). The algorithm also computes
#' the adjoint matrix of \code{C}.
#'
#' @param b numeric vector, the coefficients of the observation equation in the
#'   SSM.
#' @param C numeric matrix, the transition matrix in the state-space model.
#'
#' @return A list with two elements:
#'   \item{p}{Numeric vector containing the coefficients of the
#'     characteristic polynomial.}
#'   \item{A}{Numeric matrix representing the adjoint of \code{C}.}
#'
#' @details This function is used internally to obtain the reduced form of a
#'   state-space model.
#'  
#' @examples
#' # Example: Compute the characteristic polynomial and adjoint matrix
#' C <- matrix(c(1, 0, 1, 1), nrow = 2)
#' b <- c(1, 0)
#' lf <- .LeverrierFaddeev(b, C)
#' print(lf$p)
#' print(lf$A)
#'
#' @noRd
.LeverrierFaddeev <- function(b, C) {
  n <- nrow(C)
  m <- ncol(C)
  if (n != m)
    stop("Argument 'C' must be a square matrix.")
  if (length(b) != n)
    stop("Length of 'b' must equal number of rows in 'C'.")  
  
  p <- rep(1, n + 1)
  I <- diag(1, n)
  C1 <- I
  A <- matrix(0, n+1, n)
  A[1, ] <- b
  
  for (k in 2:(n+1)) {
    p[k] <- -sum(diag(C %*% C1)) / (k - 1)
    C1 <- C %*% C1 + p[k] * I
    A[k, ] <- b %*% C1
  }
  
  A <- A[-(n + 1), , drop = FALSE]
  list(p = p, A = A)
}

#' UCARIMA form for SS models
#'
#' \code{.ssm2ucarima} provides the UCARIMA representation of a SS model.
#'
#' @param mdl an object of class \code{\link{ssm}}.
#' @param tol numeric tolerance for unit roots.
#' 
#' @return A list with components:
#'   \code{phi} numeric vector, the AR polynomial;
#'   \code{nabla} numeric vector, the I polynomial;
#'   \code{MA} numeric matrix, the MA polynomials of each error.
#'  
#' @noRd
.ssm2ucarima <- function(mdl, tol = 1e-4) {
  stopifnot(inherits(mdl, "ssm"))
  
  lf <- .LeverrierFaddeev(mdl$b, mdl$C)
  r <- polyroot(lf$p)
  
  unit_roots <- abs(abs(r) - 1) <= tol
  if (all(!unit_roots)) {
    phi <- lf$p
    nabla <- 1
  } else if (all(unit_roots)) {
    phi <- 1
    nabla <- lf$p
  } else {
    phi <- roots2lagpol(r[!unit_roots], lp = FALSE)
    nabla <- roots2lagpol(r[unit_roots], lp = FALSE)
  }
  
  MA <- polymultC(phi, nabla)
  if (mdl$cform)
    MA <- cbind(MA, rbind(lf$A, rep(0, ncol(lf$A))))
  else
    MA <- cbind(MA, rbind(rep(0, ncol(lf$A)), lf$A))
  phi <- as.numeric(phi)
  nabla <- as.numeric(nabla)
  list(phi = phi, nabla = nabla, MA = t(MA))
}

#' ARIMA reduced form for SS models
#'
#' \code{.ssm2um} converts SSM models into univariate (ARIMA) models.
#'
#' @param mdl an object of class \code{\link{ssm}}.
#' @param tol tolerance to check if a root is close to one.
#' @param ... additional arguments.
#'
#' @return An object of class \code{\link{um}}.
#'
#' @noRd
.ssm2um <- function(mdl, tol = 1e-5, ...) {
  stopifnot(is.ssm(mdl))
  if (any(diag(mdl$S) < 0)) 
    stop("inadmissible model: negative variances")
  
  # Calculate autocovariances and Wold polynomial
  lst <- autocov.ssm(mdl, arma = FALSE, varphi = TRUE, lag.max = length(mdl$b))
  ma <- wold.pol(lst$g, type = "c", tol = tol)
  r <- polyroot(lst$varphi)
  
  # Classify roots according to distance from unit circle
  unit_roots <- abs(abs(r) - 1) <= tol
  
  # Build AR component
  if (any(!unit_roots)) {
    ar <- roots2lagpol(r[!unit_roots])$Pol
    ar <- as.lagpol(ar, coef.name = "ar")
  } else {
    ar <- NULL
  }
  
  # Combine with existing AR if present
  if (!is.null(mdl$ar)) {
    if (!is.null(ar))
      ar <- c(mdl$ar, list(ar))
    else
      ar <- mdl$ar
  } 
  
  # Build integration component
  if (any(unit_roots)) {
    i_poly <- roots2lagpol(r[unit_roots])$Pol
    i <- as.lagpol(i_poly, coef.name = "i")
  } else {
    i <- NULL
  }  
  
  # Create univariate model
  um1 <- um(bc = mdl$bc, ar = ar, i = i, 
            ma = as.lagpol(ma[-1], coef.name = "ma"),
            sig2 = ma[1], warn = FALSE)
  um1$z <- mdl$z.name
  um1 <- factorize(um1, tol = tol, ...)
  um1
}

#' Convert between contemporaneous and lagged forms of state space models
#' 
#' \code{.switchSSMform} converts a SSM from contemporaneous to lagged form and
#' vice versa.
#'
#' @param mdl an object of class \code{\link{ssm}}.
#' 
#' @return The modified SSM object.
#'
#' @noRd
.switchSSMform <- function(mdl) {
  stopifnot(inherits(mdl, "ssm"))
  if (mdl$cform) {
    s11 <- mdl$S[1, 1] + sum((mdl$b %*% mdl$S[-1, -1])*mdl$b) + 
      2*sum(mdl$b*mdl$S[1, -1])
    mdl$S[-1, 1] <- mdl$S[1, -1] <- mdl$S[1, -1] + (mdl$b %*% mdl$S[-1, -1])
    mdl$S[1, 1] <- s11
    mdl$b <- as.numeric(mdl$b %*% mdl$C)
    mdl$cform <- FALSE
  } else {
    b <- tryCatch({
      mdl$b %*% solve(mdl$C)
    }, error = function(e) {
      stop("Error inverting C matrix: ", e$message)
    })
    s11 <- mdl$S[1, 1] + sum((b %*% mdl$S[-1,-1])*b) - 2*sum(b*mdl$S[1, -1])
    mdl$S[-1, 1] <- mdl$S[1, -1] <- mdl$S[1, -1] - (b %*% mdl$S[-1, -1])
    mdl$S[1, 1] <- s11
    mdl$b <- as.numeric(b)
    mdl$cform <- TRUE
  }
  mdl$S[abs(mdl$S) < .Machine$double.eps] <- 0
  mdl
}

#' ARIMA reduced form for unobserved components
#'
#' \code{.uc2um} converts uc objects into univariate (ARIMA) models.
#'
#' @param x an object of class \code{\link{uc}}.
#' @param tol tolerance to check if a root is close to one.    
#' @param ... other arguments.
#'
#' @return An object of class \code{\link{um}}.
#'
#' @noRd
.uc2um <- function(x, tol = 1e-4, ...) {
  stopifnot(is.uc(x))
  if (x$m == 0)
    return(um(sig2 = x$S))
  lf <- .LeverrierFaddeev(x$b, x$C)
  K <- nrow(lf$A)
  g <- rep(0, K)
  for (k in 0:(K-1)) {
    sum <- 0
    for (i in (k+1):K) {
      sum <- sum + sum( (lf$A[i-k, ] %*% x$S) * lf$A[i, ] )
    }
    g[k+1] <- sum
  }
  ma <- wold.pol(g, type = "c")
  sig2 <- ma[1]
  ma <- as.lagpol(ma[-1], coef.name = "theta")
  
  r <- polyroot(lf$p)
  unit_roots <- abs(abs(r) - 1) <= tol
  if (any(!unit_roots)) {
    ar <- roots2lagpol(r[!unit_roots])$Pol
    ar <- as.lagpol(ar, coef.name = "ar")
  } else {
    ar <- NULL
  }
  
  if (any(unit_roots)) {
    i <- roots2lagpol(r[unit_roots])$Pol
    i <- as.lagpol(i, coef.name = "i")
  } else {
    i <- NULL
  }
  um1 <- um(ar = ar, i = i, ma = ma, sig2 = sig2)
  factorize(um1, ...)  
}

#' Update unobserved component
#'
#' \code{.update_uc} updates the parameters of an object of class
#' \code{\link{uc}}.
#' 
#' @param x an object of class \code{\link{uc}}.
#' @param param list with the new parameter values.
#' 
#' @return An object of class \code{\link{uc}}.
#' 
#' @noRd
.update_uc <- function(x, param) {
  nms <- names(x$param)
  nms <- nms[nms %in% names(param)]
  if (length(nms) == 0)
    return(x)    
  x$param[nms] <- param[nms] 
  x$is.adm <- TRUE
  
  if (!is.null(x$.b)) {
    x$b <- sapply(x$.b, function(z) unname(eval(z, envir = x$param)))
  } 
  
  if (!is.null(x$.C)) {
    C <- sapply(x$.C, function(z) unname(eval(z, envir = x$param)))
    x$C <- matrix(C, x$m, x$m)
    if (x$m == 1) {
      if (abs(x$C[1, 1]) >= 1) x$is.adm <- FALSE
    } else {
      if (any(abs(eigen(x$C)$values) > 1)) x$is.adm <- FALSE
    }
  } 
  
  if (!is.null(x$.S)) {
    S <- sapply(x$.S, function(z) unname(eval(z, envir = x$param)))
    if (x$m == 0) x$S <- matrix(S, 1, 1)    
    else if (x$lS <= x$m) x$S <- diag(c(S, rep(0, x$m - x$lS)), x$m, x$m)
    else x$S <- matrix(S, x$m, x$m)
    if (any(diag(x$S) < 0)) x$is.adm <- FALSE 
  } 
  return(x)
}

#' Update UC model
#'
#' \code{.update_ucm} updates the parameters of an object of class
#' \code{\link{ucm}}.
#' 
#' @param x an object of class \code{\link{ucm}}.
#' @param param list with the new parameter values.
#' 
#' @return An object of class \code{\link{ucm}}.
#' 
#' @noRd
.update_ucm <- function(x, param) {
  nms <- names(x$param)
  nms <- nms[nms %in% names(param)]
  if (length(nms) == 0)
    return(x)    
  x$param[nms] <- param[nms] 
  x$is.adm <- TRUE
  for (i in 1:x$k) {
    x$uc[[i]] <- .update_uc(x$uc[[i]], x$param)
    if (!x$uc[[i]]$is.adm) x$is.adm <- FALSE
  }  
  
  j1 <- 1
  for (i in 1:(x$k - 1)) {
    j2 <- j1 + x$uc[[i]]$m
    if (!is.null(x$uc[[i]]$.b))
      x$b[j1:(j2-1)] <- x$uc[[i]]$b
    if (!is.null(x$uc[[i]]$.C))
      x$C[j1:(j2-1), j1:(j2-1)] <- x$uc[[i]]$C
    x$S[(j1+1):j2, (j1+1):j2] <- x$uc[[i]]$S
    j1 <- j2
  }
  x$S[1, 1] <- x$uc[[x$k]]$S
  return(x)
}

#' @noRd
is.ssm <- function(x) {
  return(inherits(x, "ssm"))
}

#' @noRd
is.uc <- function(x) {
  inherits(x, "uc")
}
