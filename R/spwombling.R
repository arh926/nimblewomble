#' Posterior samples for wombling measures
#'
#' @param coords coordinates
#' @param curve coordinates of the curve for wombling
#' @param model posterior samples of \eqn{Z(s)}, \eqn{\phi}, \eqn{\sigma^2}
#' @param kernel choice of kernel; must be one of "matern1", "matern2", "gaussian"
#' @keywords sprates
#' @import nimble coda
#' @export
spwombling <- function(coords = NULL,
                       curve = NULL,
                       model = NULL,
                       kernel =  c("matern1", "matern2", "gaussian")){
  N = nrow(coords)
  distM = as.matrix(dist(coords))
  tvec = sapply(1:(nrow(curve) - 1), function(x) sqrt(sum((curve[(x + 1),] - curve[x,])^2)))
  umat = as.matrix(t(sapply(1:(nrow(curve) - 1), function(x) (curve[(x + 1),] - curve[x,])))/tvec)

  z = model[, paste("z", 1:N, sep = "")]
  phi = model[, "phi"]
  sigma2 = model[, "sigma2"]
  nmcmc = length(phi)
  if(kernel == "matern1"){
    WM1 = nimble::compileNimble(wombling_matern1)
    wmeasure = WM1(coords = coords,
                    curve = curve,
                    dists = distM,
                    tvec = tvec,
                    umat = umat,
                    z = z,
                    phi = phi,
                    sigma2 = sigma2)
  }else if(kernel == "matern2"){
    WM2 = nimble::compileNimble(wombling_matern2)
    wmeasure = WM2(coords = coords,
                    curve = curve,
                    dists = distM,
                    tvec = tvec,
                    umat = umat,
                    z = z,
                    phi = phi,
                    sigma2 = sigma2)
  }else{
    WG = nimble::compileNimble(wombling_gaussian)
    wmeasure = WG(coords = coords,
                    curve = curve,
                    dists = distM,
                    tvec = tvec,
                    umat = umat,
                    z = z,
                    phi = phi,
                    sigma2 = sigma2)
  }

  nbatch = 100

  split_samples = split(1:nmcmc, ceiling(seq_along(1:nmcmc)/(nmcmc/nbatch)))

  if(kernel == "matern1"){

    batch.wmeasure.mcmc = do.call(rbind, lapply(split_samples, function(x) apply(wmeasure[x,], 2, median)))
    estimate.wm = data.frame(round(t(apply(batch.wmeasure.mcmc, 2, quantile, probs = c(0.5, 0.025, 0.975))), 3))
    colnames(estimate.wm) = c("50%", "2.5%", "97.5%")
    estimate.wm$sig = significance(estimate.wm)

  }else{
    batch.wmeasure.mcmc = do.call(rbind, lapply(split_samples, function(x) apply(wmeasure[x,], 2, median)))
    estimate.wm.1 = data.frame(round(t(apply(batch.wmeasure.mcmc[, seq(1, 2 * (nrow(curve) - 1), by = 2)], 2, quantile, probs = c(0.5, 0.025, 0.975))), 3))
    estimate.wm.2 = data.frame(round(t(apply(batch.wmeasure.mcmc[, seq(2, 2 * (nrow(curve) - 1), by = 2)], 2, quantile, probs = c(0.5, 0.025, 0.975))), 3))
    colnames(estimate.wm.1) = colnames(estimate.wm.2) = c("50%", "2.5%", "97.5%")
    estimate.wm.1$sig = significance(estimate.wm.1)
    estimate.wm.2$sig = significance(estimate.wm.2)
  }


  if(kernel == "matern1"){
    return(list(wm.mcmc = batch.wmeasure.mcmc,
                estimate.wm = estimate.wm))
  }else{
    return(list(wm.mcmc = batch.wmeasure.mcmc,
                estimate.wm.1 = estimate.wm.1,
                estimate.wm.2 = estimate.wm.2))
  }
}
