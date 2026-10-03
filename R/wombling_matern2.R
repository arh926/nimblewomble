#' Posterior samples for wombling measures from the Matern kernel
#'  with \eqn{\nu=5/2}
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
#' @keywords wombling_matern2
#' @importFrom nimble nimbleFunction nimDim nimMatrix rmnorm_chol inverse
#' @examples
#' \dontrun{
#' #####################
#' # Internal use only #
#' #####################
#' # Example usage inside of nimblewomble::spwombling(...)
#' WM2 = compileNimble(wombling_matern2)
#' wmeasure = WM2(coords = coords,
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
wombling_matern2 <- nimble::nimbleFunction(run = function(
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
    Sigma[1:ncoords, 1:ncoords] <- materncov2(
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
      gamma[1:ncoords, 1:4] <- gamma1n2.mcov2(
        coords = coords[1:ncoords, 1:2],
        t = tvec[j],
        u = umat[j, 1:2],
        s0 = curve[j, 1:2],
        phi = phi[i]
      )

      # closed-form PATCH: same gamma-convention issue as in wombling_matern1.R. Thm 2 is
      # stated with the dimensionless G(a, z), z = sqrt(5) phi t*, but gamma_int()
      # returns b^a G(a, t*/b) with b = 1/(sqrt(5) phi). Here the two terms carry a = 1
      # and a = 2, so they were short by DIFFERENT factors (b and b^2) -- the reason the
      # error is not a single overall constant. Folding b^-a in term by term: (2
      # sqrt(5)/3) phi t* [ G(1,z) + G(2,z) ] -> (10/3) phi^2 t* Gint(a=1) + (10
      # sqrt(5)/3) phi^3 t* Gint(a=2) 10 sqrt(5) phi^3 t* G(1,z) -> 50 phi^4 t* Gint(a=1)
      # As above, this corrects only the gamma convention, not the Thm 1 (t* - |x|)
      # weight.
      k[1, 1] <- 10 /
        3 *
        phi[i]^2 *
        tvec[j] *
        gamma_int(x = tvec[j], a = 1, b = 1 / sqrt(5) / phi[i]) +
        10 *
          sqrt(5) /
          3 *
          phi[i]^3 *
          tvec[j] *
          gamma_int(x = tvec[j], a = 2, b = 1 / sqrt(5) / phi[i])
      k[2, 2] <- 50 *
        phi[i]^4 *
        tvec[j] *
        gamma_int(x = tvec[j], a = 1, b = 1 / sqrt(5) / phi[i])

      # Collapsed to one matmul: previously two separate solves against the same
      # Sigma^{-1}, one per gamma column. t(gamma[, 3:4]) is 2 x ncoords.
      tmp[1:2, 1:ncoords] <- t(gamma[1:ncoords, 3:4]) %*%
        Sigma.inv[1:ncoords, 1:ncoords]
      mean.w[1:2, 1] <- tmp[1:2, 1:ncoords] %*% z[i, 1:ncoords]
      # PATCH: gamma[, 1:2] -> gamma[, 3:4]. tmp is built from columns 3:4, and
      # gamma1n2_mcov2.R sets gamma[, 3] = -gamma[, 1] but gamma[, 4] = +gamma[, 2].
      # Using columns 1:2 on the right therefore flipped var.w[1, 1] (gradient wombling
      # variance inflated), left var.w[2, 2] accidentally correct, and made the
      # off-diagonals antisymmetric: var.w[1, 2] = -var.w[2, 1]. chol() reads one
      # triangle only, so this produced silently wrong draws rather than an error. NB do
      # not instead switch tmp to columns 1:2 -- that symmetrises var.w but flips the
      # sign of mean.w[1].
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
