% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tfm.R
\name{varsel}
\alias{varsel}
\alias{varsel.tfm}
\title{Variable selection}
\usage{
varsel(tfm, ...)

\method{varsel}{tfm}(tfm, y = NULL, p.value = 0.1, envir = NULL, ...)
}
\arguments{
\item{tfm}{a \code{tfm} object.}

\item{...}{other arguments.}

\item{y}{a "ts" object.}

\item{p.value}{probability value to decide whether or not to omit an input.}

\item{envir}{environment in which the function arguments are evaluated.
If NULL the calling environment of this function will be used.}
}
\value{
A \code{tfm} object or a "um" if no input is significant at that level.
}
\description{
\code{varsel} omits non-significant inputs from a transfer function model.
}
