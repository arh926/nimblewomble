% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wombling_matern2.R
\name{wombling_matern2}
\alias{wombling_matern2}
\title{Posterior samples for wombling measures from the Matern kernel
 with \eqn{\nu=5/2}}
\usage{
wombling_matern2(coords, curve, dists, tvec, umat, z, phi, sigma2)
}
\arguments{
\item{coords}{coordinates}

\item{curve}{curve coordinates}

\item{dists}{distance matrix}

\item{tvec}{vector of t's}

\item{umat}{matrix of u's}

\item{z}{posterior samples of \eqn{Z(s)}}

\item{phi}{posterior samples of \eqn{\phi}}

\item{sigma2}{posterior samples of \eqn{\sigma^2}}
}
\description{
For internal use only.
}
\examples{
\dontrun{
#####################
# Internal use only #
#####################
# Example usage inside of nimblewomble::spwombling(...)
WM2 = compileNimble(wombling_matern2)
wmeasure = WM2(coords = coords,
              curve = curve,
              dists = distM,
              tvec = tvec,
              umat = umat,
              z = z,
              phi = phi,
              sigma2 = sigma2)
}
}
\author{
Aritra Halder <aritra.halder@drexel.edu>, \cr
Sudipto Banerjee <sudipto@ucla.edu>
}
\keyword{wombling_matern2}
