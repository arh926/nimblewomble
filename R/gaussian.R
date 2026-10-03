#' Squared Exponential Covariance kernel
#'
#' Computes the Matern covariance matrix with fractal parameter
#' \eqn{\nu \to \infty}.
#'
#' Has the option to compute \eqn{\Sigma_{d\times d} + \tau^2 I_d}.
#'
#' @param dists distance matrix
#' @param phi spatial range
#' @param sigma2 spatial variance
#' @param tau2 nugget variance
#' @returns A matrix of covariance terms. For internal use only.
#' @keywords materncov2
#' @importFrom nimble nimbleFunction nimDim nimMatrix
#' @examples
#' \donttest{
#' #####################
#' # Internal use only #
#' #####################
#' # Used across multiple functions
#' # Example usage
#' library(nimble)
#' library(nimblewomble)
#'
#' set.seed(1)
#' # Generated Simulated Data
#' N = 1e2
#' tau = 1
#' coords = matrix(runif(2 * N, -10, 10), ncol = 2)
#' colnames(coords) = c('x', 'y')
#' dists = as.matrix(dist(coords))
#'
#' cGaussian  = compileNimble(gaussian)
#' cGaussian(dists = dists[1:N, 1:N], phi = 1, sigma2 = 1, tau2 = 0)
#' }
#' @author Aritra Halder <aritra.halder@drexel.edu>, \cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @export
gaussian <- nimble::nimbleFunction(run = function(
  dists = double(2),
  phi = double(0),
  sigma2 = double(0),
  tau2 = double(0)
) {
  returnType(double(2))
  n <- dim(dists)[1]

  result <- matrix(nrow = n, ncol = n, init = FALSE)

  # Vectorized form: nimble compiles this to BLAS-backed matrix ops, dramatically faster
  # than the previous n^2 scalar for-loop.
  result[1:n, 1:n] <- sigma2 * exp(-phi^2 * dists[1:n, 1:n]^2) + tau2 * diag(n)

  return(result)
})
