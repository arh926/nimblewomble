% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spwombling.R
\name{spwombling}
\alias{spwombling}
\title{Posterior samples for wombling measures}
\usage{
spwombling(
  coords = NULL,
  curve = NULL,
  model = NULL,
  kernel = c("matern1", "matern2", "gaussian")
)
}
\arguments{
\item{coords}{coordinates}

\item{curve}{coordinates of the curve for wombling}

\item{model}{posterior samples of \eqn{Z(s)}, \eqn{\phi}, \eqn{\sigma^2}}

\item{kernel}{choice of kernel; must be one of "matern1", "matern2", "gaussian"}
}
\description{
Posterior samples for wombling measures
}
\keyword{sprates}
