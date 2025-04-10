% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gamma1n2_gauss.R
\name{gamma1n2.gauss}
\alias{gamma1n2.gauss}
\title{Cross-covariance terms for the posterior distribution of wombling measures
for Matern \eqn{\nu\to\infty}, the squared exponential kernel.}
\usage{
gamma1n2.gauss(coords, t, u, s0, phi)
}
\arguments{
\item{coords}{coordinates}

\item{t}{value of t}

\item{u}{vector of u}

\item{s0}{starting point on curve \eqn{s_0}}

\item{phi}{posterior sample of \eqn{\phi}}
}
\description{
For internal use only. Performs one-dimensional quadrature using integral as
a limit of a sum.
}
\examples{
\dontrun{
#####################
# Internal use only #
#####################
# Example usage inside nimblewomble::wombling_gaussian(...)
gamma1n2.gauss(coords = coords[1:ncoords, 1:2], t = tvec[j],
               u = umat[j, 1:2], s0 = curve[j, 1:2], phi = phi[i])
}
}
\author{
Aritra Halder <aritra.halder@drexel.edu>, \cr
Sudipto Banerjee <sudipto@ucla.edu>
}
\keyword{gamma1n2.gauss}
