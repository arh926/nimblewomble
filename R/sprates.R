#' Posterior samples for rates of change
#'
#' @param coords coordinates
#' @param grid grid for sampling the rates of change
#' @param model posterior samples of \eqn{Z(s)}, \eqn{\phi}, \eqn{\sigma^2}
#' @param kernel choice of kernel; must be one of "matern1", "matern2", "gaussian"
#' @keywords sprates
#' @importFrom nimble compileNimble
#' @importFrom coda as.mcmc
#' @importFrom stats quantile
#' @examples
#' \dontrun{
#' require(patchwork)
#'
#' set.seed(1)
#' # Generated Simulated Data
#' N = 1e2
#' tau = 1
#' coords = matrix(runif(2 * N, -10, 10), ncol = 2); colnames(coords) = c("x", "y")
#' y = rnorm(N, mean = 20 * sin(sqrt(coords[, 1]^2  + coords[, 2]^2)), sd = tau)
#'
#' # Create equally spaced grid of points
#' xsplit = ysplit = seq(-10, 10, by = 1)[-c(1, 21)]
#' grid = as.matrix(expand.grid(xsplit, ysplit), ncol = 2)
#' colnames(grid) = c("x", "y")
#'
#' ####################################
#' # Process for True Rates of Change #
#' ####################################
#' # Gradient along x
#' true_sx = round(20 * cos(sqrt(grid[,1]^2 + grid[,2]^2)) *
#'                 grid[,1]/sqrt(grid[,1]^2 + grid[,2]^2), 3)
#' # Gradient along y
#' true_sy = round(20 * cos(sqrt(grid[,1]^2 + grid[,2]^2)) *
#'                 grid[,2]/sqrt(grid[,1]^2 + grid[,2]^2), 3)
#' # Curvature along x
#' true_sxx = round(20 * cos(sqrt(grid[,1]^2 + grid[,2]^2))/
#'                   sqrt(grid[,1]^2 + grid[,2]^2) -
#'                  20 * cos(sqrt(grid[,1]^2 + grid[,2]^2)) *
#'                  grid[,1]^2/(grid[,1]^2 + grid[,2]^2)^(3/2) -
#'                  20 * sin(sqrt(grid[,1]^2 + grid[,2]^2)) *
#'                  grid[,1]^2/(grid[,1]^2 + grid[,2]^2), 3)
#' # Mixed Curvature
#' true_sxy = round(-20 * (cos(sqrt(grid[,1]^2 + grid[,2]^2)) -
#'                  sin(sqrt(grid[,1]^2 + grid[,2]^2))) * grid[,1]
#'                   * grid[,2]/(grid[,1]^2 + grid[,2]^2), 3)
#' # Curvature along y
#' true_syy = round(20 * cos(sqrt(grid[,1]^2 + grid[,2]^2))/
#'                   sqrt(grid[,1]^2 + grid[,2]^2) -
#'                  20 * cos(sqrt(grid[,1]^2 + grid[,2]^2)) *
#'                  grid[,2]^2/(grid[,1]^2 + grid[,2]^2)^(3/2) -
#'                  20 * sin(sqrt(grid[,1]^2 + grid[,2]^2)) *
#'                  grid[,2]^2/(grid[,1]^2 + grid[,2]^2), 3)
#' # Create the plots
#' p1 = sp_ggplot(data_frame = data.frame(coords, z = y))
#' p2 = sp_ggplot(data_frame = data.frame(grid[-which(is.nan(true_sx)),],
#'                z = true_sx[-which(is.nan(true_sx))]))
#' p3 = sp_ggplot(data_frame = data.frame(grid[-which(is.nan(true_sy)),],
#'                z = true_sy[-which(is.nan(true_sy))]))
#' p4 = sp_ggplot(data_frame = data.frame(grid[-which(is.nan(true_sxx)),],
#'                z = true_sxx[-which(is.nan(true_sxx))]))
#' p5 = sp_ggplot(data_frame = data.frame(grid[-which(is.nan(true_sxy)),],
#'                z = true_sxy[-which(is.nan(true_sxy))]))
#' p6 = sp_ggplot(data_frame = data.frame(grid[-which(is.nan(true_syy)),],
#'                z = true_syy[-which(is.nan(true_syy))]))
#'
#' ((p1 + p2 + p3)/(p4 + p5 + p6))
#'
#' ##########################
#' # Fit a Gaussian Process #
#' ##########################
#' # Posterior samples for theta
#' mc_sp = gp_fit(coords = coords, y = y, kernel = "matern2")
#' # Posterior samples for Z(s) and beta
#' model = zbeta_samples(y = y, coords = coords,
#'                       model = mc_sp$mcmc,
#'                       kernel = "matern2")
#' ###################
#' # Rates of Change #
#' ###################
#' gradients = sprates(grid = grid,
#'                     coords = coords,
#'                     model = model,
#'                     kernel = "matern2")
#' p8 = sp_ggplot(data_frame = data.frame(grid,
#'                 z = gradients$estimate.sx[,"50%"],
#'                sig = gradients$estimate.sx$sig))
#' p9 = sp_ggplot(data_frame = data.frame(grid,
#'                z = gradients$estimate.sy[,"50%"],
#'                sig = gradients$estimate.sy$sig))
#' p10 = sp_ggplot(data_frame = data.frame(grid,
#'                  z = gradients$estimate.sxx[,"50%"],
#'                 sig = gradients$estimate.sxx$sig))
#' p11 = sp_ggplot(data_frame = data.frame(grid,
#'                  z = gradients$estimate.sxy[,"50%"],
#'                  sig = gradients$estimate.sxy$sig))
#' p12 = sp_ggplot(data_frame = data.frame(grid,
#'                 z = gradients$estimate.syy[,"50%"],
#'                 sig = gradients$estimate.syy$sig))
#'
#' # compare with true estimates
#' ((p7 + p8 + p9)/(p10 + p11 + p12))
#' }
#' @author Aritra Halder <aritra.halder@drexel.edu>, \cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @export
sprates <- function(coords = NULL,
                    grid = NULL,
                    model = NULL,
                    kernel =  c("matern1", "matern2", "gaussian")){
  N = nrow(coords)
  distM = as.matrix(dist(coords))
  dist.2 = sapply(1:nrow(grid), function(y) apply(coords, 1, function(x) sqrt(sum((x - grid[y,])^2))))
  dist.3 = do.call(rbind, sapply(1:nrow(grid), function(y) apply(coords, 1, function(x) (x - grid[y,])), simplify = FALSE))

  z = model[, paste0("z[", 1:N, "]")]
  phi = model[, "phi"]
  sigma2 = model[, "sigma2"]
  nmcmc = length(phi)
  if(kernel == "matern1"){
    GM1 = compileNimble(gradients_matern1)
    sprates = GM1(dists.1 = distM,
                  dists.2 = dist.2,
                  dists.3 = dist.3,
                  z = z,
                  phi = phi,
                  sigma2 = sigma2)
  }else if(kernel == "matern2"){
    CM2 = compileNimble(curvatures_matern2)
    sprates = CM2(dists.1 = distM,
                  dists.2 = dist.2,
                  dists.3 = dist.3,
                  z = z,
                  phi = phi,
                  sigma2 = sigma2)
  }else{
    CG = compileNimble(curvatures_gaussian)
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
    return(list(sx.mcmc = as.mcmc(batch.sx.mcmc),
                sy.mcmc = as.mcmc(batch.sy.mcmc),
                estimate.sx = estimate.sx,
                estimate.sy = estimate.sy))
  }else{
    return(list(sx.mcmc = as.mcmc(batch.sx.mcmc),
                sy.mcmc = as.mcmc(batch.sy.mcmc),
                sxx.mcmc = as.mcmc(batch.sxx.mcmc),
                sxy.mcmc = as.mcmc(batch.sxy.mcmc),
                syy.mcmc = as.mcmc(batch.syy.mcmc),
                estimate.sx = estimate.sx,
                estimate.sy = estimate.sy,
                estimate.sxx = estimate.sxx,
                estimate.sxy = estimate.sxy,
                estimate.syy = estimate.syy))
  }
}
