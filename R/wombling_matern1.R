#' Posterior samples for wombling measures from the Matern kernel
#'  with \eqn{\nu=3/2}
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
#' @keywords wombling_matern1
#' @importFrom nimble nimbleFunction nimDim nimMatrix inprod inverse
#' @importFrom stats rnorm
#' @examples
#' \dontrun{
#' #####################
#' # Internal use only #
#' #####################
#' # Example usage inside of nimblewomble::spwombling(...)
#' WM1 = compileNimble(wombling_matern1)
#' wmeasure = WM1(coords = coords,
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
wombling_matern1 <- nimble::nimbleFunction(run = function(
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
  result <- matrix(nrow = nmcmc, ncol = ncurve, init = FALSE)

  Sigma <- matrix(nrow = ncoords, ncol = ncoords, init = FALSE)
  Sigma.inv <- Sigma
  gamma <- matrix(nrow = ncoords, ncol = 2, init = FALSE)
  tmp <- matrix(nrow = 1, ncol = ncoords, init = FALSE)

  for (i in 1:nmcmc) {
    Sigma[1:ncoords, 1:ncoords] <- materncov1(
      dists = dists[1:ncoords, 1:ncoords],
      phi = phi[i],
      sigma2 = 1,
      tau2 = 0
    )
    # Hoisted out of the j loop: Sigma depends only on phi[i], so inverting it inside the
    # inner loop was O(ncurve * ncoords^3) wasted per sample.
    Sigma.inv[1:ncoords, 1:ncoords] <- inverse(Sigma[1:ncoords, 1:ncoords])

    for (j in 1:ncurve) {
      # cross-covariance
      gamma[1:ncoords, 1:2] <- gamma1.mcov1(
        coords = coords[1:ncoords, 1:2],
        t = tvec[j],
        u = umat[j, 1:2],
        s0 = curve[j, 1:2],
        phi = phi[i]
      )
      # closed form PATCH: Thm 2 of Halder & Banerjee (2025) reads k11 = 2 sqrt(3) sigma2
      # phi t* G(1, z) with z = sqrt(3) phi t*, where G is the DIMENSIONLESS lower
      # incomplete gamma G(a, z) = int_0^z u^(a-1) exp(-u) du (this is the form the
      # proof's intermediate expression 3 sigma2 phi^2 t* int_{-t*}^{t*} exp(-sqrt(3) phi
      # |x|) dx reduces to).  gamma_int() instead returns the DIMENSIONAL integral
      # int_0^t* s^(a-1) exp(-s/b) ds = b^a G(a, t*/b), so the original line was short by
      # a factor b = 1/(sqrt(3) phi).  Folding b^-1 in gives 2 sqrt(3) phi * sqrt(3) phi
      # = 6 phi^2.  NB this fixes only the gamma convention. The (t* - |x|) weight
      # missing from Thm 1 is a separate correction; k11 here still equals t* int f
      # rather than int (t*-|x|) f.
      k11 <- 6 *
        phi[i]^2 *
        tvec[j] *
        gamma_int(x = tvec[j], a = 1, b = 1 / sqrt(3) / phi[i])

      tmp[1, 1:ncoords] <- gamma[1:ncoords, 2] %*%
        Sigma.inv[1:ncoords, 1:ncoords]
      # inprod ensures they are scalars
      mean.w <- inprod(tmp[1, 1:ncoords], z[i, 1:ncoords])
      # PATCH: gamma[, 1] -> gamma[, 2]. tmp is built from gamma[, 2] on the line above,
      # but the quadratic form used gamma[, 1] = -gamma[, 2] (see gamma1_mcov1.R), so the
      # term was ADDED rather than subtracted and the posterior variance was inflated.
      # Both sides of the quadratic form must use the same cross-covariance vector.
      var.w <- sigma2[i] *
        (k11 - inprod(tmp[1, 1:ncoords], gamma[1:ncoords, 2]))

      # NIMBLE can only sample from N(0,1). Clamp var.w to 0 first: k11 - tmp %*% gamma
      # is a conditional variance and must be >=0 in exact arithmetic, but roundoff can
      # push it slightly negative when the curve nearly interpolates the training coords,
      # producing NaN downstream.
      result[i, j] <- sqrt(max(var.w, 0)) * rnorm(1) + mean.w
    }
  }
  return(result)
})
