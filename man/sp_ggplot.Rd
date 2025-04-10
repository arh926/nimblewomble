% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sp_plot.R
\name{sp_ggplot}
\alias{sp_ggplot}
\title{Spatial Plot Function}
\usage{
sp_ggplot(
  data_frame = NULL,
  sp = FALSE,
  shape = NULL,
  legend.key.height = 0.7,
  legend.key.width = 0.4,
  text.size = 10,
  point.size = 0.7,
  clr.pt = "black",
  palette = "Spectral",
  extend = TRUE,
  title = NULL,
  bound.box = NULL
)
}
\arguments{
\item{data_frame}{data frame consisting of coordinates and data}

\item{sp}{logical parameter indicating whether to make a spatial plot}

\item{shape}{if sp = TRUE shape file should be provided (should be an sf object)}

\item{legend.key.height}{height of legend (defaults to .7)}

\item{legend.key.width}{width of legend (defaults to .4)}

\item{text.size}{size of legend text (defaults to 10)}

\item{point.size}{size of points to be plotted (defaults to 0.7)}

\item{clr.pt}{color of point to be plotted (defaults to black)}

\item{palette}{(optional) color palette}

\item{extend}{logical parameter indicating whether to extend the interpolation (defaults to TRUE)}

\item{title}{title of the plot (defaults to NULL)}

\item{bound.box}{bounding box for spatial maps (leave as NULL if not known)}
}
\description{
Spatial Plot Function
}
\examples{
require(nimblewomble)

set.seed(1)
# Generated Simulated Data
N = 1e2
tau = 1
coords = matrix(runif(2 * N, -10, 10), ncol = 2)
colnames(coords) = c("x", "y")
y = rnorm(N, mean = 20 * sin(sqrt(coords[, 1]^2  + coords[, 2]^2)), sd = tau)

sp_ggplot(data_frame = data.frame(coords, z = y))
}
\author{
Aritra Halder <aritra.halder@drexel.edu>, \cr
Sudipto Banerjee <sudipto@ucla.edu>
}
\keyword{sp_plot}
