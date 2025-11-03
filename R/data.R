#' Retail Sales of Variety Stores (U.S. Bureau of the Census)
#'
#' 156 monthly observations from January 1967 to December 1979.
#'
#' @references 
#' Chen, C. and Liu, L. (1993) Joint Estimation of Model Parameters
#' and Outlier Effects in Time Series, Journal of the American Statistical
#' Association, Vol. 88, No. 421, pp. 284-297
#'
"rsales"
rsales <- NULL

#' Series C Chemical Process Temperature Readings: Every Minute.
#'
#' 226 observations.
#'
#' @references
#' Box, G.E., Jenkins, G.M., Reinsel, G.C. and Ljung, G.M. (2015) Time Series
#' Analysis: Forecasting and Control. John Wiley & Sons, Hoboken.
#' 
"seriesC"
seriesC <- NULL

#' Gas furnace data
#' 
#' Sampling interval 9 seconds; observations for 296 pairs of data points.
#' 
#' @format A object of class data.frame with 296 rows and 2 columns: 
#' \describe{
#' \item{X}{0.60-0.04 (input gas rate in cubir feet per minute.)}
#' \item{Y}{\% CO2 in outlet gas.}
#' }
#'
#' @references
#' Box, G.E., Jenkins, G.M., Reinsel, G.C. and Ljung, G.M. (2015) Time Series
#' Analysis: Forecasting and Control. John Wiley & Sons, Hoboken.
#' 
"seriesJ"
seriesJ <- NULL

#' Wisconsin Telephone Company
#'
#' Monthly data from January 1951 to October 1966.
#'
#' @format A object of class data.frame with 215 rows and 2 columns: 
#' \describe{
#'   \item{X}{Monthly outward station movements.} 
#'   \item{Y}{Montly inward station movements.} 
#' }
#'
#' @references 
#' 
#' Thompson, H. E. and Tiao, G. C. (1971) "Analysis of Telephone
#' Data: A Case Study of Forecasting Seasonal Time Series," Bell Journal of
#' Economics, The RAND Corporation, vol. 2(2), pages 515-541, Autumn.
#'
#' @source \url{https://drive.google.com/file/d/1LP8aMIQewMrxgOlrg9rN3eWHhZuUsY8K/view?usp=sharing}
"Wtelephone"
Wtelephone <- NULL

#' Monthly Retail Sales: Building Material and Supplies Dealers (NAICS 4441)
#'
#' Monthly U.S. retail sales for building material and supplies dealers
#' (NAICS code 4441), in millions of dollars, seasonally adjusted.
#' Source: U.S. Census Bureau, Monthly Retail Trade Survey (via FRED).
#'
#' @format A time series object of class \code{ts} with frequency 12,
#' starting in January 1992 and ending in December 2024.
#' @references U.S. Census Bureau, Monthly Retail Trade Survey.
#' FRED series code: MRTSSM4441USS.
"BuildingMat"
BuildingMat <- NULL
