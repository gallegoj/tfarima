% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stsm.R
\name{rform}
\alias{rform}
\alias{rform.stsm}
\title{Reduce form for STS model}
\usage{
rform(mdl, ...)

\method{rform}{stsm}(mdl, tol = 1e-04, ...)
}
\arguments{
\item{mdl}{an object of class \code{stsm}.}

\item{...}{other arguments.}

\item{tol}{tolerance to check if a root is close to one.}
}
\value{
An object of class \code{um}.
}
\description{
\code{rform} finds the reduce form for a STS model.
}
\examples{

b <- 1
C <- as.matrix(1)
stsm1 <- stsm(b = b, C = C, Sv = c(lvl = 1469.619), s2u = c(irr = 15103.061))
rf1 <- rform(stsm1)
nabla(rf1)
theta(rf1)
}
