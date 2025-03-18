#' Posterior samples for wombling measures from Matern with \eqn{\nu=3/2}
#'
#' @param coords coordinates
#' @param curve curve coordinates
#' @param dists distance matrix
#' @param tvec vector of t's
#' @param umat matrix of u's
#' @param z posterior samples of \eqn{Z(s)}
#' @param phi posterior samples of \eqn{\phi}
#' @param sigma2 posterior samples of \eqn{\sigma^2}
#' @keywords wombling_matern1
#' @import nimble
#' @export
wombling_matern1 <- nimble::nimbleFunction(
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
    result <- matrix(nrow = nmcmc, ncol = ncurve, init = FALSE)

    Sigma <- matrix(nrow = ncoords, ncol = ncoords, init = FALSE)
    gamma <- matrix(nrow = ncoords, ncol = 2, init = FALSE)
    tmp <- matrix(nrow = 1, ncol = ncoords, init = FALSE)

    for(i in 1:nmcmc){
      Sigma[1:ncoords, 1:ncoords] <- materncov1(dists = dists[1:ncoords, 1:ncoords],
                                                phi = phi[i],
                                                sigma2 = 1,
                                                tau2 = 0)
      for(j in 1:ncurve){
        # cross-covariance
        gamma[1:ncoords, 1:2] <- gamma1.mcov1(coords = coords[1:ncoords, 1:2],
                                              t = tvec[j],
                                              u = umat[j, 1:2],
                                              s0 = curve[j, 1:2],
                                              phi = phi[i])
        # closed form
        k11 <- 2 * sqrt(3) * phi[i] * tvec[j] * gamma_int(x = tvec[j],
                                                          a = 1,
                                                          b = 1/sqrt(3)/phi[i])

        tmp[1, 1:ncoords] <- gamma[1:ncoords, 2] %*% inverse(Sigma[1:ncoords, 1:ncoords])
        # inprod ensures they are scalars
        mean.w <- inprod(tmp[1, 1:ncoords], z[i, 1:ncoords])
        var.w <- sigma2[i] * (k11 - inprod(tmp[1, 1:ncoords], gamma[1:ncoords, 1]))

        # NIMBLE can only sample from N(0,1)
        result[i, j] <- sqrt(var.w) * rnorm(1) + mean.w
      }
      if(i == 100) cat("Iteration:", "\t", i, "...", "\t")
      if((i > 100) & (i%%100 == 0)) cat(i , "...","\t")
    }
    return(result)
  })
