#' Posterior samples of spatial effects and intercept for Matern with \eqn{nu=3/2}
#'
#' For internal use only.
#' @param y response
#' @param X covariates (supply as a matrix without intercept)
#' @param coords coordinates
#' @param model matrix of posterior samples of \eqn{\phi}, \eqn{\sigma^2} and \eqn{\tau^2}
#' @param kernel choice of kernel; must be one of "matern1", "matern2", "gaussian"
#' @keywords zbeta_samples
#' @import nimble
#' @export
zbeta_samples <- function(coords = NULL,
                          y = NULL, X = NULL,
                          model = NULL,
                          kernel = c("matern1", "matern2", "gaussian")){
  N = length(y)
  cnames = colnames(model)

  if(sum(X[,1] == 1) == N) stop("Please supply X without the first column as 1! If no covariates leave X as NULL!")
  if(is.null(X)){
    dists = as.matrix(dist(coords))

    phi = model[, "phi"]
    sigma2 = model[, "sigma2"]
    tau2 = model[, "tau2"]

    if(kernel == "matern1"){
      zbM1 = nimble::compileNimble(zbeta_matern1)
      zb.samples = zbM1(y = y, dists = dists,
                        phi = phi,
                        sigma2 = sigma2,
                        tau2 = tau2)
    }else if(kernel == "matern2"){
      zbM2 = nimble::compileNimble(zbeta_matern2)
      zb.samples = zbM2(y = y, dists = dists,
                        phi = phi,
                        sigma2 = sigma2,
                        tau2 = tau2)
    }else{
      zbG = nimble::compileNimble(zbeta_gaussian)
      zb.samples = zbG(y = y, dists = dists,
                       phi = phi,
                       sigma2 = sigma2,
                       tau2 = tau2)
    }
    model = cbind(zb.samples[, (N + 1)], model, zb.samples[, 1:N])
    colnames(model) = c("beta[0]",
                        cnames,
                        paste0("z[", 1:N, "]"))
    return(model)
  }else{
    p = ncol(X)
    beta = model[, paste0("beta[", 1:p, "]")]

    zXb = nimble::compileNimble(zXbeta)
    zb.samples = zXb(y = y, X = X, beta = beta)

    model = cbind(zb.samples[, (N + 1)], model, zb.samples[, 1:N])
    colnames(model) = c("beta[0]",
                        cnames,
                        paste0("z[", 1:N, "]"))
    return(model)
  }
}
