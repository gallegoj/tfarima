% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/um.R
\name{add.um}
\alias{add.um}
\title{Addition or substraction of univariate (ARIMA) models}
\usage{
add.um(um1, um2, add = TRUE, tol = 1e-05)
}
\arguments{
\item{um1, um2}{Two "um" S3 objects.}

\item{add}{logical. If FALSE, the second model is substracted from the first
one.}

\item{tol}{tolerance to check if a value is null.}
}
\value{
A "um" S3 object.
}
\description{
\code{add.um} creates a univariate (ARIMA) model from the addition or
substraction of two univariate (arima) models.
}
\note{
The + and - operators can also be used to add or substract ARIMA models.
}
\examples{
um1 <- um(i = "(1 - B)", ma = "(1 - 0.8B)")
um2 <- um(i = "(1 - B12)", ma = "(1 - 0.8B^12)")
um3 <- add.um(um1, um2)
um4 <- um3 - um2
}
