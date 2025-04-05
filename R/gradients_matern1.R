#' Posterior samples of rates of change (gradients) for the
#'  Matern kernel with \eqn{\nu=3/2}
#'
#' For internal use only.
#'
#' @param dists.1 distance matrix generated from coordinates
#' @param dists.2 distance of grid from coordinates
#' @param dists.3 delta = coordinate - grid
#' @param z posterior samples of \eqn{Z(s)}
#' @param phi posterior samples of \eqn{\phi}
#' @param sigma2 posterior samples of \eqn{\sigma^2}
#' @keywords gradients_matern1
#' @importFrom nimble nimbleFunction nimDim nimMatrix rmnorm_chol inverse
#' @examples
#' \dontrun{
#' #####################
#' # Internal use only #
#' #####################
#' # Example usage inside of nimblewomble::sprates()
#'  GM1 = compileNimble(gradients_matern1)
#'  sprates = GM!(dists.1 = distM,
#'               dists.2 = dist.2,
#'               dists.3 = dist.3,
#'               z = z,
#'               phi = phi,
#'               sigma2 = sigma2)
#' }
#' @author Aritra Halder <aritra.halder@drexel.edu>, \cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @export
gradients_matern1 <- nimble::nimbleFunction(
  run = function(dists.1 = double(2),
                 dists.2 = double(2),
                 dists.3 = double(2),
                 z = double(2),
                 phi = double(1),
                 sigma2 = double(1)){
    # specify return type
    returnType(double(2))
    # number of coordinates
    ncoords <- dim(dists.1)[1]
    # number of grid points
    ngrid <- dim(dists.2)[2]
    # number of posterior samples
    nmcmc <- length(phi)

    # initialize storage:
    # only spatial gradient (along x and y) are estimable
    result <- matrix(ncol =  2 * ngrid, nrow = nmcmc, init = TRUE)

    Sigma <- matrix(nrow = ncoords, ncol = ncoords, init = TRUE)
    nabla.K <- matrix(nrow = ncoords, ncol = 2, init = TRUE)
    nabla.K.t <- matrix(nrow = 2, ncol = ncoords, init = TRUE)
    tmp <- matrix(nrow = 2, ncol = ncoords, init = TRUE)

    V0 <- matrix(nrow = 2, ncol = 2, init = TRUE)
    Sigma.grad <- matrix(nrow = 2, ncol = 2, init = TRUE)
    mu.grad <- matrix(nrow = 2, ncol = 1, init = TRUE)


    for(i in 1:nmcmc){
      Sigma[1:ncoords, 1:ncoords] <- materncov1(dists = dists.1[1:ncoords,1:ncoords],
                                                phi = phi[i],
                                                sigma2 = 1,
                                                tau2 = 0)

      for(j in 1:ngrid){
        nabla.K[1:ncoords, 1] <- -3 * phi[i]^2 * exp(-phi[i] * sqrt(3) * dists.2[1:ncoords, j]) * dists.3[(2 * j - 1), 1:ncoords]
        nabla.K[1:ncoords, 2] <- -3 * phi[i]^2 * exp(-phi[i] * sqrt(3) * dists.2[1:ncoords, j]) * dists.3[(2 * j), 1:ncoords]
        V0[1:2, 1:2] <- 3 * phi[i]^2 * diag(2)
        nabla.K.t[1:2, 1:ncoords] <- -t(nabla.K[1:ncoords, 1:2])

        tmp[1:2, 1:ncoords] <- nabla.K.t[1:2, 1:ncoords] %*% inverse(Sigma[1:ncoords, 1:ncoords])
        Sigma.grad[1:2, 1:2] <- sigma2[i] * (V0[1:2, 1:2] - tmp[1:2, 1:ncoords] %*% nabla.K[1:ncoords, 1:2])
        mu.grad[1:2, 1] <- tmp[1:2, 1:ncoords] %*% z[i, 1:ncoords]


        result[i, (2 * j - 1):(2 * j)] <- rmnorm_chol(1,
                                                      mean = mu.grad[1:2, 1],
                                                      cholesky = chol(Sigma.grad[1:2, 1:2]),
                                                      prec_param = FALSE)
      }
      if(i == 500) cat("Iteration:", "\t", i, "...", "\t")
      if((i > 500) & (i%%500 == 0)) cat(i, "...","\t")
    }
    return(result)
  })
