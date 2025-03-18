#' Posterior samples of spatial effects and intercept for Matern with \eqn{nu=3/2}
#'
#' For internal use only.
#' @param y response
#' @param coords coordinates
#' @param model matrix of posterior samples of \eqn{\phi}, \eqn{\sigma^2} and \eqn{\tau^2}
#' @param kernel choice of kernel; must be one of "matern1", "matern2", "gaussian"
#' @keywords zbeta_samples
#' @import nimble
#' @export
zbeta_samples <- function(coords = NULL,
                          y = NULL,
                          model = NULL,
                          kernel = c("matern1", "matern2", "gaussian")){
  N = length(y)
  dists = as.matrix(dist(coords))
  cnames = colnames(model)

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
  model = cbind(model, zb.samples[, 1:N], zb.samples[, (N + 1)])
  colnames(model) = c(cnames,
                      paste("z", 1:N, sep = ""),
                      "b0")
  return(model)
}
