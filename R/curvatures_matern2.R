#' Posterior samples of rates of change (gradients and curvatures) for the
#'  Matern kernel with \eqn{\nu=5/2}
#'
#' For internal use only.
#'
#' @param dists.1 distance matrix generated from coordinates
#' @param dists.2 distance of grid from coordinates
#' @param dists.3 delta = coordinate - grid
#' @param z posterior samples of \eqn{Z(s)}
#' @param phi posterior samples of \eqn{\phi}
#' @param sigma2 posterior samples of \eqn{\sigma^2}
#' @returns A matrix of posterior samples for the gradient and curvatures. For internal use only.
#' @keywords curvatures_matern2
#' @importFrom nimble nimDim nimMatrix inverse rmnorm_chol nimCat
#' @examples
#' \dontrun{
#' #####################
#' # Internal use only #
#' #####################
#' # Example usage inside of nimblewomble::sprates()
#'  CM2 = compileNimble(curvatures_matern2)
#'  sprates = CM2(dists.1 = distM,
#'               dists.2 = dist.2,
#'               dists.3 = dist.3,
#'               z = z,
#'               phi = phi,
#'               sigma2 = sigma2)
#' }
#' @author Aritra Halder <aritra.halder@drexel.edu>, \cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @export
curvatures_matern2 <- nimble::nimbleFunction(run = function(
  dists.1 = double(2),
  dists.2 = double(2),
  dists.3 = double(2),
  z = double(2),
  phi = double(1),
  sigma2 = double(1)
) {
  # specify return type
  returnType(double(2))
  # number of coordinates
  ncoords <- dim(dists.1)[1]
  # number of grid points
  ngrid <- dim(dists.2)[2]
  # number of posterior samples
  nmcmc <- length(phi)

  # initialize storage: only spatial gradient (along x and y) are estimable
  result <- matrix(ncol = 5 * ngrid, nrow = nmcmc, init = FALSE)

  Sigma <- matrix(nrow = ncoords, ncol = ncoords, init = FALSE)
  Sigma.inv <- Sigma
  nabla.K <- matrix(nrow = ncoords, ncol = 5, init = FALSE)
  nabla.K.t <- matrix(nrow = 5, ncol = ncoords, init = FALSE)
  tmp <- matrix(nrow = 5, ncol = ncoords, init = FALSE)

  V0 <- matrix(nrow = 5, ncol = 5, init = TRUE)
  Sigma.grad <- matrix(nrow = 5, ncol = 5, init = FALSE)
  mu.grad <- matrix(nrow = 5, ncol = 1, init = FALSE)

  for (i in 1:nmcmc) {
    Sigma[1:ncoords, 1:ncoords] <- materncov2(
      dists = dists.1[1:ncoords, 1:ncoords],
      phi = phi[i],
      sigma2 = 1,
      tau2 = 0
    )
    # Hoisted out of the j loop: Sigma depends only on phi[i], so inverting it inside the
    # inner loop was O(ngrid * ncoords^3) wasted per sample.
    Sigma.inv[1:ncoords, 1:ncoords] <- inverse(Sigma[1:ncoords, 1:ncoords])

    # Also hoisted: everything below depends only on phi[i], not on j.
    phi2 <- phi[i]^2
    phi4 <- phi2 * phi2
    V0[1, 1] <- 5 / 3 * phi2
    V0[2, 2] <- 5 / 3 * phi2
    V0[3, 3] <- 25 * phi4
    V0[3, 5] <- 25 / 3 * phi4
    V0[4, 4] <- 25 / 3 * phi4
    V0[5, 3] <- 25 / 3 * phi4
    V0[5, 5] <- 25 * phi4

    for (j in 1:ngrid) {
      # Reused per-j vectors — exp() previously fired 5x per grid point.
      dx <- dists.3[(2 * j - 1), 1:ncoords]
      dy <- dists.3[2 * j, 1:ncoords]
      exp.term <- exp(-phi[i] * sqrt(5) * dists.2[1:ncoords, j])
      lin.term <- 1 + sqrt(5) * phi[i] * dists.2[1:ncoords, j]

      nabla.K[1:ncoords, 1] <- -5 / 3 * phi2 * lin.term * exp.term * dx
      nabla.K[1:ncoords, 2] <- -5 / 3 * phi2 * lin.term * exp.term * dy
      nabla.K[1:ncoords, 3] <- -5 /
        3 *
        phi2 *
        (lin.term - 5 * phi2 * dx^2) *
        exp.term
      # PATCH: phi[i]^2 -> phi[i]^4. d^2K/dg_x dg_y = (25/3) phi^4 exp(-sqrt(5) phi d) dx
      # dy for Matern nu=5/2. Consistent with V0[4, 4] = 25/3 * phi^4 above. Verified by
      # finite-differencing the kernel.
      nabla.K[1:ncoords, 4] <- 25 / 3 * phi4 * exp.term * dx * dy
      nabla.K[1:ncoords, 5] <- -5 /
        3 *
        phi2 *
        (lin.term - 5 * phi2 * dy^2) *
        exp.term

      nabla.K.t[1:2, 1:ncoords] <- -t(nabla.K[1:ncoords, 1:2])
      nabla.K.t[3:5, 1:ncoords] <- t(nabla.K[1:ncoords, 3:5])

      tmp[1:5, 1:ncoords] <- nabla.K.t[1:5, 1:ncoords] %*%
        Sigma.inv[1:ncoords, 1:ncoords]
      # PATCH: right factor was nabla.K rather than t(nabla.K.t). Rows 1:2 of nabla.K.t
      # are negated and rows 3:5 are not, so with A = nabla.K and S = diag(-1,-1,1,1,1)
      # the old line computed V0 - S A' Sigma^-1 A instead of V0 - S A' Sigma^-1 A S.
      # That matrix is NOT SYMMETRIC: blocks (1:2, 3:5) and (3:5, 1:2) carry opposite
      # signs, and the gradient block is inflated rather than reduced. chol() reads only
      # one triangle, so this failed silently rather than erroring. Curvature variances
      # were unaffected.
      Sigma.grad[1:5, 1:5] <- sigma2[i] *
        (V0[1:5, 1:5] - tmp[1:5, 1:ncoords] %*% t(nabla.K.t[1:5, 1:ncoords]))
      mu.grad[1:5, 1] <- tmp[1:5, 1:ncoords] %*% z[i, 1:ncoords]

      result[i, (5 * j - 4):(5 * j)] <- rmnorm_chol(
        1,
        mean = mu.grad[1:5, 1],
        cholesky = chol(Sigma.grad[1:5, 1:5]),
        prec_param = FALSE
      )
    }
  }
  return(result)
})
