#' Posterior samples of spatial effects and intercept for all kernels in the
#' presence of covariates
#'
#' For internal use only.
#'
#' @param y response
#' @param X covariates (supply as a matrix without intercept)
#' @param beta  posterior samples of \eqn{\beta} (supply as a matrix)
#' @keywords zXbeta
#' @importFrom nimble nimbleFunction nimMatrix
#' @examples
#' \dontrun{
#' #####################
#' # Internal use only #
#' #####################
#' # Example usage inside of nimblewomble::zbeta_samples(...)
#' zXb = nimble::compileNimble(zXbeta)
#' zb.samples = zXb(y = y, X = X, beta = beta)
#' }
#' @author Aritra Halder <aritra.halder@drexel.edu>, \cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @export
zXbeta <- nimble::nimbleFunction(
  run = function(y = double(1),
                 X = double(2),
                 beta = double(2)){
    returnType(double(2))
    n <- length(y)
    p <- dim(X)[2]
    nmcmc <- dim(beta)[1]

    result <- matrix(ncol = n + 1, nrow = nmcmc, init = FALSE)

    for(i in 1:nmcmc){
      result[i, 1:n] <- y[1:n] - X[1:n, 1:p] %*% beta[i, 1:p]
      result[i, (n + 1)] <- mean(result[i, 1:n])
      result[i, 1:n] <- result[i, 1:n] - result[i, (n + 1)]
    }
    return(result)
  })

