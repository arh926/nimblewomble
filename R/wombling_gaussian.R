#' Posterior samples for wombling measures for the squared exponential kernel
#'
#' For internal use only.
#'
#' @param coords coordinates
#' @param curve curve coordinates
#' @param dists distance matrix
#' @param tvec vector of t's
#' @param umat matrix of u's
#' @param z posterior samples of \eqn{Z(s)}
#' @param phi posterior samples of \eqn{\phi}
#' @param sigma2 posterior samples of \eqn{\sigma^2}
#' @returns A matrix of wombling measures. For internal use only.
#' @keywords wombling_gaussian
#' @importFrom nimble nimbleFunction nimDim nimMatrix rmnorm_chol inverse
#' @examples
#' \dontrun{
#' #####################
#' # Internal use only #
#' #####################
#' # Example usage inside of nimblewomble::spwombling(...)
#' WG = compileNimble(wombling_gaussian)
#' wmeasure = WG(coords = coords,
#'               curve = curve,
#'               dists = distM,
#'               tvec = tvec,
#'               umat = umat,
#'               z = z,
#'               phi = phi,
#'               sigma2 = sigma2)
#' }
#' @author Aritra Halder <aritra.halder@drexel.edu>, \cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @export
wombling_gaussian <- nimble::nimbleFunction(run = function(
  coords = double(2),
  curve = double(2),
  dists = double(2),
  tvec = double(1),
  umat = double(2),
  z = double(2),
  phi = double(1),
  sigma2 = double(1)
) {
  # specify return type
  returnType(double(2))
  # number of coordinates
  ncoords <- dim(coords)[1]
  # number of posterior samples
  nmcmc <- length(phi)
  # curve length
  ncurve <- length(tvec)

  # initialize storage:
  result <- matrix(nrow = nmcmc, ncol = 2 * ncurve, init = FALSE)

  Sigma <- matrix(nrow = ncoords, ncol = ncoords, init = FALSE)
  Sigma.inv <- Sigma
  gamma <- matrix(nrow = ncoords, ncol = 4, init = FALSE)
  tmp <- matrix(nrow = 2, ncol = ncoords, init = FALSE)
  k <- matrix(nrow = 2, ncol = 2, init = TRUE)
  mean.w <- matrix(nrow = 2, ncol = 1, init = FALSE)
  var.w <- matrix(nrow = 2, ncol = 2, init = FALSE)

  for (i in 1:nmcmc) {
    Sigma[1:ncoords, 1:ncoords] <- gaussian(
      dists = dists[1:ncoords, 1:ncoords],
      phi = phi[i],
      sigma2 = 1,
      tau2 = 0
    )
    # Hoisted out of the j loop: Sigma depends only on phi[i], so inverting it TWICE per
    # inner iteration was O(2 * ncurve * ncoords^3) wasted.
    Sigma.inv[1:ncoords, 1:ncoords] <- inverse(Sigma[1:ncoords, 1:ncoords])

    for (j in 1:ncurve) {
      # cross-covariance
      gamma[1:ncoords, 1:4] <- gamma1n2.gauss(
        coords = coords[1:ncoords, 1:2],
        t = tvec[j],
        u = umat[j, 1:2],
        s0 = curve[j, 1:2],
        phi = phi[i]
      )

      # closed-form
      k[1, 1] <- 2 *
        sqrt(pi) *
        phi[i] *
        tvec[j] *
        (2 * pnorm_nimble(sqrt(2) * phi[i] * tvec[j]) - 1)
      k[2, 2] <- k[1, 1] * 6 * phi[i]^2

      # Collapsed to one matmul: previously two separate solves against the same
      # Sigma^{-1}, one per gamma column. t(gamma[, 3:4]) is 2 x ncoords.
      tmp[1:2, 1:ncoords] <- t(gamma[1:ncoords, 3:4]) %*%
        Sigma.inv[1:ncoords, 1:ncoords]
      mean.w[1:2, 1] <- tmp[1:2, 1:ncoords] %*% z[i, 1:ncoords]
      # PATCH: gamma[, 1:2] -> gamma[, 3:4], as in wombling_matern2.R. tmp is built from
      # columns 3:4; gamma1n2_gauss.R sets gamma[, 3] = -gamma[, 1], gamma[, 4] =
      # +gamma[, 2].
      var.w[1:2, 1:2] <- sigma2[i] *
        (k[1:2, 1:2] - tmp[1:2, 1:ncoords] %*% gamma[1:ncoords, 3:4])

      result[i, (2 * j - 1):(2 * j)] <- rmnorm_chol(
        1,
        mean = mean.w[1:2, 1],
        cholesky = chol(var.w[1:2, 1:2]),
        prec_param = FALSE
      )
    }
  }
  return(result)
})
