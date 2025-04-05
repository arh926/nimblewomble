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
wombling_matern2 <- nimble::nimbleFunction(
  run = function(coords = double(2),
                 curve = double(2),
                 dists = double(2),
                 tvec = double(1),
                 umat = double(2),
                 z = double(2),
                 phi = double(1),
                 sigma2 = double(1)){
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
    gamma <- matrix(nrow = ncoords, ncol = 4, init = FALSE)
    tmp <- matrix(nrow = 2, ncol = ncoords, init = FALSE)
    k <- matrix(nrow = 2, ncol = 2, init = TRUE)
    mean.w <- matrix(nrow = 2, ncol = 1, init = FALSE)
    var.w <- matrix(nrow = 2, ncol = 2, init = FALSE)

    for(i in 1:nmcmc){
      Sigma[1:ncoords, 1:ncoords] <- materncov2(dists = dists[1:ncoords, 1:ncoords],
                                                phi = phi[i],
                                                sigma2 = 1,
                                                tau2 = 0)

      for(j in 1:ncurve){
        # cross-covariance
        gamma[1:ncoords, 1:4] <- gamma1n2.mcov2(coords = coords[1:ncoords, 1:2],
                                                t = tvec[j],
                                                u = umat[j, 1:2],
                                                s0 = curve[j, 1:2],
                                                phi = phi[i])

        # closed-form
        k[1, 1] <- 2 * sqrt(5)/3 * phi[i] * tvec[j] * (gamma_int(x = tvec[j],
                                                                 a = 1,
                                                                 b = 1/sqrt(5)/phi[i]) +
                                                         gamma_int(x = tvec[j],
                                                                   a = 2,
                                                                   b = 1/sqrt(5)/phi[i]))
        k[2, 2] <- 10 * sqrt(5) * phi[i]^3 * tvec[j] * gamma_int(x = tvec[j],
                                                                 a = 1,
                                                                 b = 1/sqrt(5)/phi[i])

        tmp[1, 1:ncoords] <- gamma[1:ncoords, 3] %*% inverse(Sigma[1:ncoords, 1:ncoords])
        tmp[2, 1:ncoords] <- gamma[1:ncoords, 4] %*% inverse(Sigma[1:ncoords, 1:ncoords])
        mean.w[1:2, 1] <- tmp[1:2, 1:ncoords] %*% z[i, 1:ncoords]
        var.w[1:2, 1:2] <- sigma2[i] * (k[1:2, 1:2] - tmp[1:2, 1:ncoords] %*% gamma[1:ncoords, 1:2])


        result[i, (2 * j - 1):(2 * j)] <- rmnorm_chol(1, mean = mean.w[1:2, 1], cholesky = chol(var.w[1:2, 1:2]), prec_param = FALSE)
      }
      if(i == 100) cat("Iteration:", "\t", i, "...", "\t")
      if((i > 100) & (i%%100 == 0)) cat(i , "...","\t")
    }
    return(result)
  })
