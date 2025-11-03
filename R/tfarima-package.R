#' Transfer Function and ARIMA Models
#'
#' The \pkg{tfarima} package provides classes and methods to build customized
#' transfer function and ARIMA models with multiple operators and parameter
#' restrictions. It includes functions for model identification, estimation using
#' exact or conditional maximum likelihood, diagnostic checking, automatic outlier
#' detection, calendar effects, forecasting, and seasonal adjustment.
#'
#' The current version extends the functionality by incorporating the estimation
#' of unobserved components in ARIMA models through the UCARIMA representation
#' and structural time series models. 
#'
#' @name tfarima-package
#' @aliases tfarima
#' @author
#' Jose Luis Gallego \email{jose.gallego@@unican.es}
#'
#' @references
#' Bell, W. R. and Hillmer, S. C. (1983).
#' Modeling Time Series with Calendar Variation.
#' \emph{Journal of the American Statistical Association}, 78(383), 526–534.  
#'
#' Box, G. E. P., Jenkins, G. M., Reinsel, G. C., and Ljung, G. M. (2015).
#' \emph{Time Series Analysis: Forecasting and Control}.
#' John Wiley & Sons, Hoboken.  
#'
#' Box, G. E. P., Pierce, D. A., and Newbold, D. A. (1987).
#' Estimating Trend and Growth Rates in Seasonal Time Series.
#' \emph{Journal of the American Statistical Association}, 82(397), 276–282.  
#'
#' Box, G. E. P. and Tiao, G. C. (1975).
#' Intervention Analysis with Applications to Economic and Environmental Problems.
#' \emph{Journal of the American Statistical Association}, 70(349), 70–79.  
#'
#' Chen, C. and Liu, L. (1993).
#' Joint Estimation of Model Parameters and Outlier Effects in Time Series.
#' \emph{Journal of the American Statistical Association}, 88(421), 284–297.  
#'
#' Thompson, H. E. and Tiao, G. C. (1971).
#' Analysis of Telephone Data: A Case Study of Forecasting Seasonal Time Series.
#' \emph{Bell Journal of Economics}, 2(2), 515–541.
#'
#' @importFrom Rcpp evalCpp
#' @importFrom numDeriv jacobian
#' @importFrom MASS ginv
#' @importFrom stats acf AIC as.ts bartlett.test BIC Box.test cycle density  
#' @importFrom stats dnorm end frequency is.ts is.mts lm logLik median na.pass  
#' @importFrom stats optim optimize plot.ts pnorm predict printCoefmat qnorm   
#' @importFrom stats resid residuals rnorm sd start time ts tsdiag update var 
#' @importFrom stats window rbinom toeplitz
#' @importFrom graphics abline axis hist layout lcm legend lines mtext par
#' @importFrom graphics plot plot.new points segments text title 
#' @importFrom utils tail
#' @useDynLib tfarima, .registration = TRUE
NULL
