% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zbeta_matern2.R
\name{zbeta_matern2}
\alias{zbeta_matern2}
\title{Posterior samples of spatial effects and intercept for the Matern kernel with
 \eqn{\nu=5/2}}
\usage{
zbeta_matern2(y, dists, phi, sigma2, tau2)
}
\arguments{
\item{y}{response}

\item{dists}{distance matrix derived from coordinates}

\item{phi}{posterior samples of \eqn{\phi}}

\item{sigma2}{posterior samples of \eqn{\sigma^2}}

\item{tau2}{posterior samples of \eqn{\tau^2}}
}
\description{
For internal use only.
}
\examples{
\dontrun{
#####################
# Internal use only #
#####################
# Example usage inside of nimblewomble::zbeta_samples(...)
zbM2 = compileNimble(zbeta_matern2)
zb.samples = zbM2(y = y, dists = dists, phi = phi, sigma2 = sigma2,
                 tau2 = tau2)
}
}
\author{
Aritra Halder <aritra.halder@drexel.edu>, \cr
Sudipto Banerjee <sudipto@ucla.edu>
}
\keyword{zbeta_matern2}
