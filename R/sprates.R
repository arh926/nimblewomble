#' Posterior samples for rates of change
#'
#' @param coords coordinates
#' @param grid grid for sampling the rates of change
#' @param model posterior samples of \eqn{Z(s)}, \eqn{\phi}, \eqn{\sigma^2}
#' @param kernel choice of kernel; must be one of "matern1", "matern2", "gaussian"
#' @keywords sprates
#' @import nimble coda
#' @export
sprates <- function(coords = NULL,
                    grid = NULL,
                    model = NULL,
                    kernel =  c("matern1", "matern2", "gaussian")){
  N = nrow(coords)
  distM = as.matrix(dist(coords))
  dist.2 = sapply(1:nrow(grid), function(y) apply(coords, 1, function(x) sqrt(sum((x - grid[y,])^2))))
  dist.3 = do.call(rbind, sapply(1:nrow(grid), function(y) apply(coords, 1, function(x) (x - grid[y,])), simplify = FALSE))

  z = model[, paste("z", 1:N, sep = "")]
  phi = model[, "phi"]
  sigma2 = model[, "sigma2"]
  nmcmc = length(phi)
  if(kernel == "matern1"){
    GM1 = nimble::compileNimble(gradients_matern1)
    sprates = GM1(dists.1 = distM,
                  dists.2 = dist.2,
                  dists.3 = dist.3,
                  z = z,
                  phi = phi,
                  sigma2 = sigma2)
  }else if(kernel == "matern2"){
    CM2 = nimble::compileNimble(curvatures_matern2)
    sprates = CM2(dists.1 = distM,
                  dists.2 = dist.2,
                  dists.3 = dist.3,
                  z = z,
                  phi = phi,
                  sigma2 = sigma2)
  }else{
    CG = nimble::compileNimble(curvatures_gaussian)
    sprates = CG(dists.1 = distM,
                 dists.2 = dist.2,
                 dists.3 = dist.3,
                 z = z,
                 phi = phi,
                 sigma2 = sigma2)
  }
  nbatch = 100

  split_samples = split(1:nmcmc, ceiling(seq_along(1:nmcmc)/(nmcmc/nbatch)))
  if(kernel == "matern1"){
    sx.mcmc = sprates[, seq(1, (2 * nrow(grid)), by = 2)]
    sy.mcmc = sprates[, seq(2, (2 * nrow(grid)), by = 2)]


    batch.sx.mcmc = do.call(rbind, lapply(split_samples, function(x) apply(sx.mcmc[x,], 2, median)))
    batch.sy.mcmc = do.call(rbind, lapply(split_samples, function(x) apply(sy.mcmc[x,], 2, median)))

    estimate.sx = data.frame(round(t(apply(batch.sx.mcmc, 2, quantile, probs = c(0.5, 0.025, 0.975))), 3))
    estimate.sy = data.frame(round(t(apply(batch.sy.mcmc, 2, quantile, probs = c(0.5, 0.025, 0.975))), 3))

    colnames(estimate.sx) = colnames(estimate.sy) = c("50%", "2.5%", "97.5%")
    estimate.sx$sig = significance(data_frame = estimate.sx)
    estimate.sy$sig = significance(data_frame = estimate.sy)

  }else{
    sx.mcmc = sprates[, seq(1, (5 * nrow(grid)), by = 5)]
    sy.mcmc = sprates[, seq(2, (5 * nrow(grid)), by = 5)]
    sxx.mcmc = sprates[, seq(3, (5 * nrow(grid)), by = 5)]
    sxy.mcmc = sprates[, seq(4, (5 * nrow(grid)), by = 5)]
    syy.mcmc = sprates[, seq(5, (5 * nrow(grid)), by = 5)]


    batch.sx.mcmc = do.call(rbind, lapply(split_samples, function(x) apply(sx.mcmc[x,], 2, median)))
    batch.sy.mcmc = do.call(rbind, lapply(split_samples, function(x) apply(sy.mcmc[x,], 2, median)))
    batch.sxx.mcmc = do.call(rbind, lapply(split_samples, function(x) apply(sxx.mcmc[x,], 2, median)))
    batch.sxy.mcmc = do.call(rbind, lapply(split_samples, function(x) apply(sxy.mcmc[x,], 2, median)))
    batch.syy.mcmc = do.call(rbind, lapply(split_samples, function(x) apply(syy.mcmc[x,], 2, median)))


    estimate.sx = data.frame(round(t(apply(batch.sx.mcmc, 2, quantile, probs = c(0.5, 0.025, 0.975))), 3))
    estimate.sy = data.frame(round(t(apply(batch.sy.mcmc, 2, quantile, probs = c(0.5, 0.025, 0.975))), 3))
    estimate.sxx = data.frame(round(t(apply(batch.sxx.mcmc, 2, quantile, probs = c(0.5, 0.025, 0.975))), 3))
    estimate.sxy = data.frame(round(t(apply(batch.sxy.mcmc, 2, quantile, probs = c(0.5, 0.025, 0.975))), 3))
    estimate.syy = data.frame(round(t(apply(batch.syy.mcmc, 2, quantile, probs = c(0.5, 0.025, 0.975))), 3))

    colnames(estimate.sx) = colnames(estimate.sy) = colnames(estimate.sxx) = colnames(estimate.sxy) = colnames(estimate.syy) = c("50%", "2.5%", "97.5%")
    estimate.sx$sig = significance(data_frame = estimate.sx)
    estimate.sy$sig = significance(data_frame = estimate.sy)
    estimate.sxx$sig = significance(data_frame = estimate.sxx)
    estimate.sxy$sig = significance(data_frame = estimate.sxy)
    estimate.syy$sig = significance(data_frame = estimate.syy)
  }


  if(kernel == "matern1"){
    return(list(sx.mcmc = coda::as.mcmc(batch.sx.mcmc),
                sy.mcmc = coda::as.mcmc(batch.sy.mcmc),
                estimate.sx = estimate.sx,
                estimate.sy = estimate.sy))
  }else{
    return(list(sx.mcmc = coda::as.mcmc(batch.sx.mcmc),
                sy.mcmc = coda::as.mcmc(batch.sy.mcmc),
                sxx.mcmc = coda::as.mcmc(batch.sxx.mcmc),
                sxy.mcmc = coda::as.mcmc(batch.sxy.mcmc),
                syy.mcmc = coda::as.mcmc(batch.syy.mcmc),
                estimate.sx = estimate.sx,
                estimate.sy = estimate.sy,
                estimate.sxx = estimate.sxx,
                estimate.sxy = estimate.sxy,
                estimate.syy = estimate.syy))
  }
}
