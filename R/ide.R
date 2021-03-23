
#' Identification plots
#' 
#' \code{ide} displays graphs useful to identify a tentative ARIMA 
#' model for a time series.
#'
#' @param Y Univariate or multivariate time series.
#' @param transf Data transformations, list(bc = F, d = 0, D = 0, S = F),
#' where bc is the Box-Cox logarithmic transformation, d and D are the
#' number of nonseasonal and seasonal differences, and S is the annual 
#' sum operator.
#' @param order.polreg an integer indicating the order of a polynomial trend.
#' @param lag.max number of autocorrelations.
#' @param wn.bands logical. If TRUE confidence intervals for sample 
#' autocorrelations are computed assuming a white noise series. 
#' @param graphs graphs to be shown: plot, hist, acf, pacf, pgram, 
#' cpgram (cummulative periodogram), rm (range-median).
#' @param set.layout logical. If TRUE the layout is set by the function,
#' otherwise it is set by the user.
#' @param byrow logical. If TRUE the layout is filled by rows, otherwise it 
#' is filled by columns.
#' @param main title of the graph.
#' @param ... additional arguments.
#' @return NULL
#' 
#' @examples
#' Y <- AirPassengers
#' ide(Y, graphs = c("plot", "rm"))
#' ide(Y, transf = list(list(bc = TRUE, S = TRUE), list(bc = TRUE, d = 1, D = 1)))
#' 
#' @export
ide <- function(Y, transf = list(), order.polreg = 0, lag.max = NULL, 
                wn.bands = TRUE, graphs = c("plot", "acf", "pacf"), 
                set.layout = TRUE, byrow = TRUE, main = "", ...) {
  

  ylab <- deparse(substitute(Y))
  if (!exists(ylab, envir = .GlobalEnv)) ylab <- "y"
  
  args <- list(...)
  graphs <- tolower(graphs)
  graphs <- unique(graphs)
  graphs <- match.arg(graphs, c("plot", "hist", "acf", "pacf", "pgram", "cpgram", "rm", "sdm"), several.ok = TRUE)
  n.graphs <- length(graphs)
  n.transf <- length(transf)
  if (n.transf == 0) {
    n.transf <- 1
    transf <- list(list())
  } else {
    nm <- names(transf)
    if (is.null(nm)) {
      if(!all(sapply(transf, is.list))) stop("wrong list of transformations")
    } else {
      if (all(nm %in% c("bc", "d", "D", "S", "i", "std"))) {
        n.transf <- 1
        transf <- list(transf)
      }
    }
  }

  if (is.matrix(Y)) n.ser <- ncol(Y)
  else if (is.list(Y)) n.ser <- length(Y)
  else if (is.numeric(Y)) n.ser <- 1
  else stop("Y must be a ts object")  

  if (set.layout) {
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))  
    
    if (n.ser == 1 && !is.list(Y)) {
      if (n.transf == 1) {
        if (graphs[1] ==  "plot") {
          if (n.graphs ==  1) m <- matrix(1, nrow = 1, ncol = 1, byrow = TRUE)
          else if (n.graphs ==  2) m <- matrix(c(1, 1, 2), nrow = 1, ncol = 3, byrow = TRUE)
          else if (n.graphs ==  3) m <- matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2, byrow = byrow)
          else if (n.graphs ==  4) m <- matrix(c(1, 1, 1, 2, 3, 4), nrow = 2, ncol = 3, byrow = TRUE)
          else if (n.graphs ==  5) m <- matrix(c(1, 1, 2, 3, 4, 5), nrow = 2, ncol = 3, byrow = TRUE)
          else if (n.graphs ==  6) m <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 4, byrow = TRUE)
          else if (n.graphs ==  7) m <- matrix(c(1, 1, 1, 1, 2, 3, 4, 5, 6), nrow = 2, ncol = 4, byrow = TRUE)
          else stop("error: too many graphs")
        } else {
          if (n.graphs < 4) m <- matrix(1:n.graphs, nrow = 1, ncol = n.graphs, byrow = TRUE)
          else m <- matrix(1:n.graphs, nrow = 2, ncol = n.graphs/2 +  (n.graphs%%2!= 0)*1, byrow = TRUE)
        } 
      } else {
        if (graphs[1] ==  "plot") {
          j <- -n.graphs + 1
          m <- sapply(1:n.transf, function(tr) {
            j <<- j + n.graphs
            c(j , j:(j+n.graphs-1))
          })
          if (byrow) m <- t(m)
        } else {
          m <- matrix(1:(n.graphs*n.transf), nrow = n.transf, ncol = n.graphs, byrow = TRUE)
          if (!byrow) m <- t(m)
        }
      }
    } else if (n.graphs > 1) { # n.ser > 1
      m <- matrix(0, nrow = n.ser*n.transf, ncol = n.graphs)
      k <- 1
      for (tr in 1:n.transf) {
        for (ser in 1:n.ser) {
          row <- (ser - 1)*n.transf + tr
          for (gr in 1:n.graphs) {
            m[row, gr] <- k
            k <- k + 1
          }
        }
      }
      if (!byrow) m <- t(m) 
      else if (graphs[1] == "plot") m <- cbind(m[, 1], m)
    } else {
      if (n.transf > 1) {
        m <- matrix(1:(n.ser*n.transf), nrow = n.ser, ncol = n.transf, byrow = TRUE)
        if (!byrow) m <- t(m)
      } else {
        if (n.ser < 4) {
          if (byrow) m <- matrix(1:n.ser, nrow = 1, ncol = n.ser)
          else m <- matrix(1:n.ser, nrow = n.ser, ncol = 1)
        } else m <- matrix(1:n.ser, nrow = 2, ncol = n.ser/2 +  (n.ser%%2!= 0)*1, byrow = TRUE)
      }
    }
    
    if (main !=  "") {
      par(oma = c(0, 0, 2.5, 1.5), mar = c(2.5,2.5,1.5,0.8), mgp = c(1.5,0.6,0))
    }
    else {
      par(mar = c(3, 3, 1.5, 0.5), mgp = c(1.5, 0.6, 0))
    }
    layout(m)
  }
  
  if (n.ser > 1) ylab1 <- ylab
  for (ser in 1:n.ser) {
    for (tr in 1:n.transf) {
      if (is.matrix(Y)) y <- Y[, ser]
      else if (is.list(Y)) y <- Y[[ser]]
      else y <- Y 
      if (n.ser > 1) ylab <- paste(ylab1, ser, sep = "") 
      if (!is.ts(y)) y <- ts(y)
      s <- frequency(y)
      
      bc <- FALSE; d <- 0L; D <- 0L; i <- NULL; S <- FALSE; std <- FALSE;  
      if (!is.null(transf[[tr]]$bc)) bc <- as.logical(transf[[tr]]$bc)
      if (!is.null(transf[[tr]]$d)) d <- as.integer(transf[[tr]]$d)
      if (!is.null(transf[[tr]]$D)) D <- as.integer(transf[[tr]]$D)
      if (!is.null(transf[[tr]]$S)) S <- as.logical(transf[[tr]]$S)
      if (!is.null(transf[[tr]]$i)) i <- transf[[tr]]$i
      if (!is.null(transf[[tr]]$std)) std <- as.logical(transf[[tr]]$std)
      
      if (bc & min(y)>0) y <- log(y)
      else bc <- FALSE
      if (d > 0) y <- diff(y, differences = d)
      if (D > 0) y <- diff(y, lag = s, differences = D)
      if (S & s > 1) y <- ts(diffC(y, rep(1, s), FALSE), end = end(y), frequency = frequency(y))

      if (!is.null(i)) {
        if (inherits(i, "lagpol")) i <- i$pol
        else if (is.lagpol.list(i)) i <- polyexpand(i)
        else stop("i must be a list of lag polynomials")
        y <- ts(diffC(y, i, FALSE), end = end(y), frequency = frequency(y))
        ylab <- "w"
      }  
    
      n <- length(y)
      stopifnot(n > 1)
      if (is.null(lag.max) ) {
        if (s > 1) lag.max <- min(3*s+3, n/4)
        else lag.max <- floor(n/4)
      } else if (lag.max < 1 | lag.max > n-1) {
        if (s > 1) lag.max <- min(3*s+3, n/4)
        else lag.max <- floor(n/4)
      }
      stopifnot(lag.max > 0)
    
      if (std) { 
        y <- (y - mean(y))/sd(y)
        maxy <- ceiling(max(y))
        miny <- floor(min(y))
        my <- 0
        sy <- 1
      } else {
        my <- mean(y)
        sy <- sd(y)
        maxy <- my + ceiling(max(y - my)/sy)*sy
        miny <- my - abs(floor(min(y - my)/sy))*sy
      }
    
      for (j in 1:n.graphs) { 
        if (graphs[j] ==  "plot") {
          plot(y, xlab = "t", ylab = tslabel(ylab, bc, d, D, s, S), type = "n", col = "black", ylim = c(miny, maxy))
          abline(h = seq(miny, maxy, sy), lty = 2, col = "gray")
          abline(h = my, col = "gray")  
          lines(y, type = "l")
          if (order.polreg > 0) {
            t <- time(y)
            N <- t[length(t)]
            X <- sapply(1:order.polreg, function(r) (t/N)^r)
            reg <- lm(y~X)
            points(t, reg$fitted.values, ty = "l", col = "gray", lty = 2)
          }
        }
        else if (graphs[j] ==  "acf") {
          rho <- stats::acf(y, lag.max = lag.max, plot = F)$acf[, ,]
          rho <- rho[-1]
          k <- 1:length(rho)
          if (max(abs(rho)) < 0.5) ylim = c(-0.5, 0.5)
          else ylim = c(-1, 1)
          plot(k, rho, type = "n" ,xlab = "lag", ylab = "ACF", ylim = ylim, xaxt = "n")
          if (wn.bands) abline(h = c(-1.96/sqrt(n), 1.96/sqrt(n)), lty = 2, col = "gray")
          else {
            se <- sapply(k, function(x) {
              sqrt((1+2*sum(rho[1:x]^2))/n)
            })
            lines(k, 1.96*se, lty = 2, col = "gray")
            lines(k, -1.96*se, lty = 2, col = "gray")
          }
          abline(h = 0)  
          #
          if (s>1 & lag.max > s) { 
            abline(v = seq(s, lag.max, s), lty = 2, col = "gray" )
            axis(1, at = seq(s, lag.max, s))
          } else if (lag.max > 5) {
            abline(v = seq(5, lag.max, 5), lty = 2, col = "gray" )
            axis(1, at = seq(5, lag.max, 5))
          }
          lines(k, rho, type = "h")  
        }
        else if (graphs[j] ==  "pacf") {
          phi <- stats::pacf(y, lag.max = lag.max, plot = F)$acf[, ,]
          k <- 1:length(phi)
          if (max(abs(phi)) < 0.5) ylim = c(-0.5, 0.5)
          else ylim = c(-1, 1)
          plot(k, phi, type = "n", xlab = "lag", ylab = "PACF", ylim = ylim, xaxt = "n")
          abline(h = 0)
          abline(h = c(-1.96/sqrt(n), 1.96/sqrt(n)), lty = 2, col = "gray")  
          if (s>1 & lag.max > s) {
            abline(v = seq(s, lag.max, s), lty = 2, col = "gray" )
            axis(1, at = seq(s, lag.max, s))
          } else if (lag.max > 5) {
            abline(v = seq(5, lag.max, 5), lty = 2, col = "gray" )
            axis(1, at = seq(5, lag.max, 5))
          }
          lines(k, phi, type = "h")
        } else if (graphs[j] ==  "hist") {
          fy <- density(y)
          x <- seq(miny, maxy, length = 100)
          fx <- dnorm(x, mean(y), sd(y))
          if (j ==  2 & graphs[1] ==  "plot") {
            plot(fy$y, fy$x, ylim = c(miny, maxy), xlim = c(0, max(max(fy$y), max(fx))), type = "l",
                 ylab = tslabel(ylab, bc, d, D, s, S), xlab = "Density")
            lines(fx, x, col = "gray")
          } else {
            plot(x, fx, xlim = c(miny, maxy), ylim = c(0, max(max(fy$y), max(fx))), type = "l",
                 xlab = tslabel(ylab, bc, d, D, s, S), ylab = "Density", col = "gray")
            lines(fy$x, fy$y)
          }
        } else if (graphs[j] ==  "pgram") {
          P <- pgramC(y, FALSE)
          plot(P[, 1], P[, 2], type = "l", xlim = c(0, 0.5), ylab = "Periodogram", xlab = "Frequency")
          #lines(x = c(0, 6), z = c(0, 1), col = "gray", lty = 1)
          title(ylab = "Periodogram")
        } else if (graphs[j] ==  "cpgram") {
          P <- pgramC(y, TRUE)
          plot(P[, 1], P[, 2], type = "l", ylim = c(0, 1), xlim = c(0, 0.5), ylab = "Cum. per.", xlab = "Frequency")
          abline(a = 0, b = 2, col = "gray", lty = 1)
          q <- ifelse(n%%2 ==  0, (n-2)/2, (n-1)/2)  
          abline(a = 1.36/sqrt(q), b = 2, col = "gray", lty = 2)
          abline(a = -1.36/sqrt(q), b = 2, col = "gray", lty = 2)
        } else if (graphs[j] ==  "rm" || graphs[j] ==  "sdm") {
          if (is.null(args[["ng"]])) {
            if (s > 1) ng <- n/(2*s)
            else ng <- n/10
          } else ng <- args[["ng"]]
          if (ng < 2) ng <- 2
          t <- rep(1:ng, each = n%/%ng)
          n1 <- length(t)
          if (n1 < n) t <- c(t, rep(ng, n - n1))
          t <- t[1:n]
          rm <- graphs[j] ==  "rm"
          m <- sapply(1:ng, function(x) {
            if (rm) median(y[t == x])
            else mean(y[t == x])
          })
          r <- sapply(min(t):max(t), function(x) {
            if (rm) (max(y[t == x]) - min(y[t == x]))
            else sd(y[t == x])
          })
          if (rm) {Ylab = "range"; Xlab = "median"}
          else {Ylab = "st. dev."; Xlab = "mean"}
          plot(m, r, ylab = Ylab, xlab = Xlab)
          abline(lm(r~m), col = "gray", lty = 2)
        } 
      }
    }
      
    if (main !=  "") 
      title(main = main, outer = TRUE, cex.main = 1.10)
  }

  invisible(NULL)
  
}

tslabel <- function(ylab, bc, d, D, s, S) {
  
  if (S & s < 2) S <- FALSE
  
  if (S & d > 0) {
    S <- FALSE
    d <- d - 1
    D <- D + 1
  }

  if (bc) {
    if (d > 0) {
      if (D > 0) {
        if (d > 1) {
          if (D > 1) substitute(nabla^d*nabla^D[s]*log*(ylab[t]), list(ylab=ylab, d = d, D = D, s = s))
          else substitute(nabla^d*nabla[s]*log*(ylab[t]), list(ylab=ylab, d = d,s = s))
        } else {
          if (D > 1) substitute(nabla*nabla^D[s]*log*(ylab[t]), list(ylab=ylab, D = D, s = s))
          else substitute(nabla*nabla[s]*log*(ylab[t]), list(ylab=ylab, s = s))
        }
      } else {
        if (d > 1) substitute(nabla^d*log*(ylab[t]), list(ylab=ylab, d = d))
        else substitute(nabla*log*(ylab[t]), list(ylab=ylab))
      }
    } else if (D > 0) {
      if (S) {
        if (D > 1) substitute(nabla[s]^D*S[s]*log*(ylab[t]), list(ylab=ylab, D = D, s = s))
        else substitute(nabla[s]*S[s]*log*(ylab[t]), list(ylab=ylab, s = s))   
      } else {
        if (D > 1) substitute(nabla[s]^D*log*(ylab[t]), list(ylab=ylab, D = D, s = s))
        else substitute(nabla[s]*log*(ylab[t]), list(ylab=ylab, s = s))
      }
    } else if (S) {
      substitute(S[s]*log*(ylab[t]), list(ylab=ylab, s = s))
    } else {
      substitute(log*(ylab[t]), list(ylab=ylab))
    }
  } else {
    if (d > 0) {
      if (D > 0) {
        if (d > 1) {
          if (D > 1) substitute(nabla^d*nabla^D[s]*ylab[t], list(ylab=ylab, d = d, D = D, s = s))
          else substitute(nabla^d*nabla[s]*ylab[t], list(ylab=ylab, d = d,s = s))
        } else {
          if (D > 1) substitute(nabla*nabla^D[s]*ylab[t], list(ylab=ylab, D = D, s = s))
          else substitute(nabla*nabla[s]*ylab[t], list(ylab=ylab, s = s))
        }
      } else {
        if (d > 1) substitute(nabla^d*ylab[t], list(ylab=ylab, d = d))
        else substitute(nabla*ylab[t], list(ylab=ylab))
      }
    } else if (D > 0) {
      if (S) {
        if (D > 1) substitute(nabla[s]^D*S[s]*ylab[t], list(ylab=ylab, D = D, s = s))
        else substitute(nabla[s]*S[s]*ylab[t], list(ylab=ylab, s = s))   
      } else {
        if (D > 1) substitute(nabla[s]^D*ylab[t], list(ylab=ylab, D = D, s = s))
        else substitute(nabla[s]*ylab[t], list(ylab=ylab, s = s))
      }
    } else if (S) {
      substitute(S[s]*ylab[t], list(ylab=ylab, s = s))
    } else {
      substitute(ylab[t], list(ylab=ylab))
    }
  }
}


