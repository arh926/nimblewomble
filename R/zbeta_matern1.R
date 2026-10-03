#' Posterior samples of spatial effects and intercept for the Matern kernel
#' with \eqn{nu=3/2}
#'
#' For internal use only.
#'
#' @param y response
#' @param dists distance matrix derived from coordinates
#' @param phi posterior samples of \eqn{\phi}
#' @param sigma2 posterior samples of \eqn{\sigma^2}
#' @param tau2 posterior samples of \eqn{\tau^2}
#' @returns A matrix of spatial effects and intercept. For internal use only.
#' @keywords zbeta_matern1
#' @importFrom nimble nimbleFunction nimMatrix inverse rmnorm_chol
#' @examples
#' \dontrun{
#' #####################
#' # Internal use only #
#' #####################
#' # Example usage inside of nimblewomble::zbeta_samples(...)
#' zbM1 = compileNimble(zbeta_matern1)
#' zb.samples = zbM1(y = y, dists = dists, phi = phi, sigma2 = sigma2,
#'                  tau2 = tau2)
#' }
#' @author Aritra Halder <aritra.halder@drexel.edu>, \cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @export
zbeta_matern1 <- nimble::nimbleFunction(run = function(
  y = double(1),
  dists = double(2),
  phi = double(1),
  sigma2 = double(1),
  tau2 = double(1)
) {
  returnType(double(2))
  n <- length(y)
  nmcmc <- length(phi)

  result <- matrix(ncol = n + 1, nrow = nmcmc, init = FALSE)

  Sigma <- matrix(nrow = n, ncol = n, init = FALSE)
  M <- Sigma
  M_inv <- Sigma
  mu <- matrix(nrow = n, ncol = 1, init = FALSE)

  for (i in 1:nmcmc) {
    # Woodbury rewrite of the posterior variance Sigma = ( K^{-1}/sigma2 + I/tau2 )^{-1}
    # = tau2*I - sigma2*tau2^2 * (K + sigma2*tau2*I)^{-1} One solve on the
    # well-conditioned SPD M = K + sigma2*tau2*I, instead of the original double inverse
    # (invert K, then invert the sum) which both cost O(n^3) and amplified K's
    # ill-conditioning when phi was small.
    M[1:n, 1:n] <- materncov1(
      dists = dists[1:n, 1:n],
      phi = phi[i],
      sigma2 = 1,
      tau2 = 0
    )
    for (k in 1:n) {
      M[k, k] <- M[k, k] + sigma2[i] * tau2[i]
    }
    M_inv[1:n, 1:n] <- inverse(M[1:n, 1:n])

    Sigma[1:n, 1:n] <- -sigma2[i] * tau2[i]^2 * M_inv[1:n, 1:n]
    for (k in 1:n) {
      Sigma[k, k] <- Sigma[k, k] + tau2[i]
    }

    mu[1:n, 1] <- Sigma[1:n, 1:n] %*% y[1:n] / tau2[i]

    result[i, 1:n] <- rmnorm_chol(
      1,
      mean = mu[1:n, 1],
      cholesky = chol(Sigma[1:n, 1:n]),
      prec_param = FALSE
    )
    result[i, (n + 1)] <- mean(result[i, 1:n])
    result[i, 1:n] <- result[i, 1:n] - result[i, (n + 1)]
  }
  return(result)
})
