% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stsm.R
\name{bsm}
\alias{bsm}
\title{Basic Structural Time Series models}
\usage{
bsm(
  y,
  bc = FALSE,
  seas = c("hd", "ht", "hs"),
  par = c(irr = 0.75, lvl = 1, slp = 0.05, seas = 0.075),
  fixed = c(lvl = TRUE),
  xreg = NULL,
  fit = TRUE,
  updmdl = NULL,
  ...
)
}
\arguments{
\item{y}{an object of class \code{ts}, with frequency 4 or 12.}

\item{bc}{logical. If TRUE logs are taken.}

\item{seas}{character, type of seasonality (Harvey-Durbin (hd), Harvey-Todd
(ht), Harrison-Steven (ht))}

\item{par}{real vector with the error variances of each unobserved
component.}

\item{fixed}{logical vector to fix parameters.}

\item{xreg}{matrix of regressors.}

\item{fit}{logical. If TRUE, model is fitted.}

\item{updmdl}{function to update the parameters of the BSM.}

\item{...}{additional arguments.}
}
\value{
An object of class \code{stsm}.
}
\description{
\code{bsm} creates/estimates basic structural models for seasonal time
series.
}
\examples{

bsm1 <- bsm(AirPassengers, bc = TRUE)

}
\references{
Durbin, J. and Koopman, S.J. (2012) Time Series Analysis by State Space
Methods, 2nd ed., Oxford University Press, Oxford.
}
