% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/materncov1.R
\name{materncov1}
\alias{materncov1}
\title{Matern Covariance kernel with \eqn{\nu = 3/2}}
\usage{
materncov1(dists, phi, sigma2, tau2)
}
\arguments{
\item{dists}{distance matrix}

\item{phi}{spatial range}

\item{sigma2}{spatial variance}

\item{tau2}{nugget variance}
}
\description{
Computes the Matern covariance matrix with fractal parameter \eqn{\nu = 3/2}.
Has the option to compute \eqn{\Sigma_{d\times d} + \tau^2 I_d}.
}
\examples{
\dontrun{
#####################
# Internal use only #
#####################
# Used across multiple functions
# Example usage
materncov1(dists = dists[1:ncoords, 1:ncoords], phi = phi[i], sigma2 = 1, tau2 = 0)
}
}
\author{
Aritra Halder <aritra.halder@drexel.edu>, \cr
Sudipto Banerjee <sudipto@ucla.edu>
}
\keyword{materncov1}
