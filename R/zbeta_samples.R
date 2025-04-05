#' Posterior samples of spatial effects and intercept for Matern with \eqn{nu=3/2}
#'
#' For internal use only.
#' @param y response
#' @param X covariates (supply as a matrix without intercept)
#' @param coords coordinates
#' @param model matrix of posterior samples of \eqn{\phi}, \eqn{\sigma^2} and \eqn{\tau^2}
#' @param kernel choice of kernel; must be one of "matern1", "matern2", "gaussian"
#' @keywords zbeta_samples
#' @importFrom nimble compileNimble
#' @examples
#' \dontrun{
#' set.seed(1)
#' # Generated Simulated Data
#' N = 1e2
#' tau = 1
#' coords = matrix(runif(2 * N, -10, 10), ncol = 2); colnames(coords) = c("x", "y")
#' y = rnorm(N, mean = 20 * sin(sqrt(coords[, 1]^2  + coords[, 2]^2)), sd = tau)
#'
#' # Posterior samples for theta
#' mc_sp = gp_fit(coords = coords, y = y, kernel = "matern2")
#' # Posterior samples for Z(s) and beta
#' model = zbeta_samples(y = y, coords = coords,
#'                       model = mc_sp$mcmc,
#'                       kernel = "matern2")
#' estimates = t(round(apply(model, 2, quantile, probs = c(0.5, 0.025, 0.975)), 3))
#' yfit = estimates[paste("z", 1:N, sep = ""), "50%"] + estimates["b0", "50%"]
#' ylow = estimates[paste("z", 1:N, sep = ""), "2.5%"] + estimates["b0", "2.5%"]
#' yhigh = estimates[paste("z", 1:N, sep = ""), "97.5%"] + estimates["b0", "97.5%"]
#' fit_frame = cbind(true = y, est = yfit, `2.5%` = ylow, `97.5%` = yhigh)
#' fit_frame$sig = significance(data_frame = data.frame(fit_frame[,-1]))
#'
#' # Plot
#' sp_ggplot(data_frame = data.frame(coords, z = yfit, sig = fit_frame$sig))
#' }
#' @author Aritra Halder <aritra.halder@drexel.edu>, \cr Sudipto Banerjee <sudipto@ucla.edu>
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
      zbM1 = compileNimble(zbeta_matern1)
      zb.samples = zbM1(y = y, dists = dists,
                        phi = phi,
                        sigma2 = sigma2,
                        tau2 = tau2)
    }else if(kernel == "matern2"){
      zbM2 = compileNimble(zbeta_matern2)
      zb.samples = zbM2(y = y, dists = dists,
                        phi = phi,
                        sigma2 = sigma2,
                        tau2 = tau2)
    }else{
      zbG = compileNimble(zbeta_gaussian)
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
