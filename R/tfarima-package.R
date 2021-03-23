#' Transfer Function and ARIMA Models.
#'
#' The tfarima package provides classes and methods to build customized transfer
#' function and ARIMA models with multiple operators and parameter restrictions.
#' The package also includes functions for model identification, model
#' estimation (exact or conditional maximum likelihood), model diagnostic
#' checking, automatic outlier detection, calendar effects, forecasting and
#' seasonal adjustment.
#'
#' @name tfarima-package
#' @aliases tfarima
#' @author Jose Luis Gallego \email{jose.gallego@@unican.es}
#'
#' @references
#'
#' Bell, W.R. and Hillmer, S.C. (1983) Modeling Time Series with Calendar
#' Variation, Journal of the American Statistical Association, Vol. 78, No. 383,
#' pp. 526-534.
#'
#' Box, G.E., Jenkins, G.M., Reinsel, G.C. and Ljung, G.M. (2015) Time Series
#' Analysis: Forecasting and Control. John Wiley & Sons, Hoboken.
#'
#' Box, G.E.P.,  Pierce, D.A. and Newbold, D. A. (1987) Estimating Trend and
#' Growth Rates in Seasonal Time Series, Journal of the American Statistical
#' Association, Vol. 82, No. 397, pp. 276-282.
#'
#' Box, G.E.P. and  Tiao,  G.C. (1975) “Intervention Analysis with Applications
#' to Economic and Environmental Problems”, Journal of the American Statistical
#' Association, Vol. 70, No. 349, pp. 70-79.
#'
#' Chen, C. and Liu, L. (1993) Joint Estimation of Model Parameters and Outlier
#' Effects in Time Series, Journal of the American Statistical Association, Vol.
#' 88, No. 421, pp. 284-297
#'
#' Thompson, H. E. and Tiao, G. C. (1971) "Analysis of Telephone
#' Data: A Case Study of Forecasting Seasonal Time Series," Bell Journal of
#' Economics, The RAND Corporation, vol. 2(2), pages 515-541, Autumn.
#'
#' @keywords package
#' @docType package
#' @importFrom Rcpp evalCpp
#' @importFrom numDeriv jacobian
#' @importFrom stats Box.test acf as.ts bartlett.test cycle density dnorm end
#' @importFrom stats frequency is.ts lm median na.pass optim plot.ts pnorm
#' @importFrom stats printCoefmat qnorm resid residuals rnorm sd start time
#' @importFrom stats ts tsdiag update var window
#' @importFrom graphics plot abline axis layout lcm legend lines mtext par
#' @importFrom graphics plot.new points text title
#' @useDynLib tfarima, .registration=TRUE
NULL