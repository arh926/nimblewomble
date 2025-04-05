#' Matern Covariance with \eqn{\nu = 5/2}
#'
#' Computes the Matern covariance matrix with fractal parameter \eqn{\nu = 5/2}.
#' @param dists distance matrix
#' @param phi spatial range
#' @param sigma2 spatial variance
#' @param tau2 nugget variance
#' @keywords materncov2
#' @importFrom nimble nimbleFunction nimDim nimMatrix
#' @export
materncov2 <- nimble::nimbleFunction(
  run = function(dists = double(2),
                 phi = double(0),
                 sigma2 = double(0),
                 tau2 = double(0)){

    returnType(double(2))
    n <- dim(dists)[1]

    result <- matrix(nrow = n, ncol = n, init = FALSE)

    # result[1:n, 1:n] <- sigma2 *  (1 + phi * sqrt(5) * dists[1:n, 1:n] + phi^2 * 5 * dists[1:n, 1:n]^2/3) * exp(-phi * sqrt(5) * dists[1:n, 1:n]) + tau2 * diag(n)
    for(i in 1:n)
      for(j in 1:n){
        result[i, j] <- sigma2 *  (1 + phi * sqrt(5) * dists[i, j] + phi^2 * 5 * dists[i, j]^2/3) * exp(-phi * sqrt(5) * dists[i, j])
        if(i == j)  result[i, j] <- result[i, j] + tau2
        else result[i, j] <- result[i, j]
      }

    return(result)
  })
