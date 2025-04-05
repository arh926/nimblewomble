#' Posterior samples of rates of change for Matern with \eqn{\nu=5/2}
#'
#' For internal use only.
#' @param dists.1 distance matrix generated from coordinates
#' @param dists.2 distance of grid from coordinates
#' @param dists.3 delta = coordinate - grid
#' @param z posterior samples of \eqn{Z(s)}
#' @param phi posterior samples of \eqn{\phi}
#' @param sigma2 posterior samples of \eqn{\sigma^2}
#' @keywords curvatures_gaussian
#' @importFrom nimble nimbleFunction nimDim nimMatrix inverse rmnorm_chol nimCat
#' @export
curvatures_gaussian <- nimble::nimbleFunction(
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
    result <- matrix(ncol =  5 * ngrid, nrow = nmcmc, init = FALSE)

    Sigma <- matrix(nrow = ncoords, ncol = ncoords, init = FALSE)
    nabla.K <- matrix(nrow = ncoords, ncol = 5, init = FALSE)
    nabla.K.t <- matrix(nrow = 5, ncol = ncoords, init = FALSE)
    tmp <- matrix(nrow = 5, ncol = ncoords, init = FALSE)

    V0 <- matrix(nrow = 5, ncol = 5, init = TRUE)
    Sigma.grad <- matrix(nrow = 5, ncol = 5, init = FALSE)
    mu.grad <- matrix(nrow = 5, ncol = 1, init = FALSE)


    for(i in 1:nmcmc){
      Sigma[1:ncoords, 1:ncoords] <- gaussian(dists = dists.1[1:ncoords,1:ncoords],
                                              phi = phi[i],
                                              sigma2 = 1,
                                              tau2 = 0)

      for(j in 1:ngrid){
        nabla.K[1:ncoords, 1] <- -2 * phi[i]^2 * exp(-phi[i]^2 * dists.2[1:ncoords, j]^2) * dists.3[(2 * j - 1), 1:ncoords]
        nabla.K[1:ncoords, 2] <- -2 * phi[i]^2 * exp(-phi[i]^2 * dists.2[1:ncoords, j]^2) * dists.3[2 * j, 1:ncoords]
        nabla.K[1:ncoords, 3] <- -2 * phi[i]^2 * (1 - 2 * phi[i]^2 * dists.3[(2 * j - 1), 1:ncoords]^2) * exp(-phi[i]^2 * dists.2[1:ncoords, j]^2)
        nabla.K[1:ncoords, 4] <- 4 * phi[i]^2 * exp(-phi[i]^2 * dists.2[1:ncoords, j]^2) * dists.3[(2 * j - 1), 1:ncoords] * dists.3[2 * j, 1:ncoords]
        nabla.K[1:ncoords, 5] <- -2 * phi[i]^2 * (1 - 2 * phi[i]^2 * dists.3[2 * j, 1:ncoords]^2) * exp(-phi[i]^2 * dists.2[1:ncoords, j]^2)
        V0[1, 1] <- 2 * phi[i]^2
        V0[2, 2] <- 2 * phi[i]^2
        V0[3, 3] <- 12 * phi[i]^4
        V0[3, 5] <- 4 * phi[i]^4
        V0[4, 4] <- 4 * phi[i]^4
        V0[5, 3] <- 4 * phi[i]^4
        V0[5, 5] <- 12 * phi[i]^4

        nabla.K.t[1:2, 1:ncoords] <- -t(nabla.K[1:ncoords, 1:2])
        nabla.K.t[3:5, 1:ncoords] <- t(nabla.K[1:ncoords, 3:5])

        tmp[1:5, 1:ncoords] <- nabla.K.t[1:5, 1:ncoords] %*% inverse(Sigma[1:ncoords, 1:ncoords])
        Sigma.grad[1:5, 1:5] <- sigma2[i] * (V0[1:5, 1:5] - tmp[1:5, 1:ncoords] %*% nabla.K[1:ncoords, 1:5])
        mu.grad[1:5, 1] <- tmp[1:5, 1:ncoords] %*% z[i, 1:ncoords]


        result[i, (5 * j - 4):(5 * j)] <- rmnorm_chol(1,
                                                      mean = mu.grad[1:5, 1],
                                                      cholesky = chol(Sigma.grad[1:5, 1:5]),
                                                      prec_param = FALSE)
      }
      if(i == 500) cat("Iteration:", "\t", i, "...", "\t")
      if((i > 500) & (i%%500 == 0)) cat(i, "...","\t")
    }
    return(result)
  })
