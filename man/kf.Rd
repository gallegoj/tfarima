% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stsm.R
\name{kf}
\alias{kf}
\alias{kf.stsm}
\title{Kalman filter for STS models}
\usage{
kf(mdl, ...)

\method{kf}{stsm}(mdl, y = NULL, x1 = NULL, P1 = NULL, filtered = FALSE, ...)
}
\arguments{
\item{mdl}{an object of class \code{stsm}.}

\item{...}{additional arguments.}

\item{y}{time series to be filtered when it differs from the model series.}

\item{x1}{initial state vector.}

\item{P1}{covariance matrix of x1.}

\item{filtered}{logical. If TRUE, the filtered states x_{t|t} and their covariances
matrices P_{t|t} are returned. Otherwise, x_{t|t-1} and P_{t|t-1} are}
}
\value{
An list with the innovations, the conditional states and their
  covariance matrices.
}
\description{
\code{kf} computes the innovations and the conditional states with the Kalman
filter algorithm.
}
