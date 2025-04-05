#' Posterior samples of spatial effects and intercept for the squared
#' exponential kernel
#'
#' For internal use only.
#'
#' @param y response
#' @param dists distance matrix derived from coordinates
#' @param phi posterior samples of \eqn{\phi}
#' @param sigma2 posterior samples of \eqn{\sigma^2}
#' @param tau2 posterior samples of \eqn{\tau^2}
#' @keywords zbeta_gaussian
#' @importFrom nimble nimbleFunction nimMatrix inverse rmnorm_chol
#' @examples
#' \dontrun{
#' #####################
#' # Internal use only #
#' #####################
#' # Example usage inside of nimblewomble::zbeta_samples(...)
#' zbG = compileNimble(zbeta_gaussian)
#' zb.samples = zbG(y = y, dists = dists, phi = phi, sigma2 = sigma2,
#'                  tau2 = tau2)
#' }
#' @author Aritra Halder <aritra.halder@drexel.edu>, \cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @export
zbeta_gaussian <- nimble::nimbleFunction(
  run = function(y = double(1),
                 dists = double(2),
                 phi = double(1),
                 sigma2 = double(1),
                 tau2 = double(1)){
    returnType(double(2))
    n <- length(y)
    nmcmc <- length(phi)

    result <- matrix(ncol = n + 1, nrow = nmcmc, init = FALSE)

    Sigma <- matrix(nrow = n, ncol = n, init = FALSE)
    mu <- matrix(nrow = n, ncol = 1, init = FALSE)

    for(i in 1:nmcmc){
      Sigma[1:n, 1:n] <- inverse(inverse(gaussian(dists = dists[1:n,1:n],
                                                  phi = phi[i],
                                                  sigma2 = 1,
                                                  tau2 = 0))/sigma2[i] + diag(n)/tau2[i])
      mu[1:n, 1] <- Sigma[1:n, 1:n] %*% y[1:n]/tau2[i]

      result[i, 1:n] <- rmnorm_chol(1,
                                    mean = mu[1:n, 1],
                                    cholesky = chol(Sigma[1:n, 1:n]),
                                    prec_param = FALSE)
      result[i, (n + 1)] <- mean(result[i, 1:n])
      result[i, 1:n] <- result[i, 1:n] - result[i, (n + 1)]
    }
    return(result)
  })
