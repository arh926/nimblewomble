#' Posterior samples for wombling measures
#'
#' @param coords coordinates
#' @param curve coordinates of the curve for wombling
#' @param model posterior samples of \eqn{Z(s)}, \eqn{\phi}, \eqn{\sigma^2}
#' @param kernel choice of kernel; must be one of "matern1", "matern2", "gaussian"
#' @keywords sprates
#' @importFrom nimble compileNimble
#' @importFrom coda as.mcmc
#' @importFrom stats median quantile
#' @examples
#' \dontrun{
#' require(nimble)
#' require(nimblewomble)
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
#' ##########################
#' # Fit a Gaussian Process #
#' ##########################
#' # Posterior samples for theta
#' mc_sp = gp_fit(coords = coords, y = y, kernel = "matern2")
#' # Posterior samples for Z(s) and beta
#' model = zbeta_samples(y = y, coords = coords,
#'                       model = mc_sp$mcmc,
#'                       kernel = "matern2")
#' ############
#' # Wombling #
#' ############
#' # Pick any curve (contour) of your choice
#' # curve = your contour
#' tvec = sapply(1:(nrow(curve) - 1), function(x) sqrt(sum((curve[(x + 1),] - curve[x,])^2)))
#' umat = as.matrix(t(sapply(1:(nrow(curve) - 1), function(x) (curve[(x + 1),] - curve[x,])))/tvec)
#'
#' wm = spwombling(coords = coords,
#'                 curve = curve,
#'                 model = model,
#'                 kernel = "matern2")
#'
#' # Total wombling measure for gradient
#' colSums(wm$estimate.wm.1[,-4]); colSums(wm$estimate.wm.1[,-4])/sum(tvec)
#' # Total wombling measure for curvature
#' colSums(wm$estimate.wm.2[,-4]); colSums(wm$estimate.wm.2[,-4])/sum(tvec)
#'
#'
#' # Color code points based on significance
#' col.pts.1 = sapply(wm$estimate.wm.1$sig, function(x){
#'   if(x == 1) return("green")
#'   else if(x == -1) return("cyan")
#'   else return(NA)
#'   })
#'
#'   col.pts.2 = sapply(wm$estimate.wm.2$sig, function(x){
#'   if(x == 1) return("green")
#'   else if(x == -1) return("cyan")
#'   else return(NA)
#'   })
#'
#' p13 = sp_ggplot(data_frame = data.frame(coords, y))
#' p14 = p13 + geom_path(curve, mapping = aes(x, y), linewidth = 2)
#' p15 = p13 + geom_path(curve, mapping = aes(x, y), linewidth = 2) +
#' geom_path(curve, mapping = aes(x, y),
#'           colour = c(col.pts.1, NA), linewidth = 1, na.rm = TRUE)
#' p16 = p13 + geom_path(curve, mapping = aes(x, y), linewidth = 2) +
#'             geom_path(curve, mapping = aes(x, y),
#'             colour = c(col.pts.2, NA), linewidth = 1, na.rm = TRUE)
#'
#' p14 + (p15/p16)
#'
#' ###############
#' # True Values #
#' ###############
#' truth = matrix(0, nrow = nrow(curve) - 1, ncol = 2)
#' rule = seq(0, 1, by = 0.01)
#'
#' for(i in 1:(nrow(curve) - 1)){
#'   u.perp = c(umat[i, 2], - umat[i, 1])
#'   s0 = curve[i,]
#'
#' truth.lsegment = sapply(rule * tvec[i], function(x){
#'   s.t = s0 + x * umat[i,]
#'   true_sx = 20 * cos(sqrt(s.t[1]^2 + s.t[2]^2)) * s.t[1]/
#'             sqrt(s.t[1]^2 + s.t[2]^2)
#'   true_sy = 20 * cos(sqrt(s.t[1]^2 + s.t[2]^2)) * s.t[2]/
#'             sqrt(s.t[1]^2 + s.t[2]^2)
#'   true_sx * u.perp[1] + true_sy * u.perp[2]
#'   })
#'
#' truth[i, 1] = sum(truth.lsegment * (tvec[i]/101))
#'
#' truth.lsegment = sapply(rule * tvec[i], function(x){
#'   s.t = s0 + x * umat[i,]
#'   true_sxx = 20 * cos(sqrt(s.t[1]^2 + s.t[2]^2))/sqrt(s.t[1]^2 + s.t[2]^2) -
#'     20 * cos(sqrt(s.t[1]^2 + s.t[2]^2)) *
#'             s.t[1]^2/(s.t[1]^2 + s.t[2]^2)^(3/2) -
#'     20 * sin(sqrt(s.t[1]^2 + s.t[2]^2)) * s.t[1]^2/(s.t[1]^2 + s.t[2]^2)
#'   true_sxy = -20 * (cos(sqrt(s.t[1]^2 + s.t[2]^2)) -
#'                 sin(sqrt(s.t[1]^2 + s.t[2]^2))) *
#'                       s.t[1] * s.t[2]/(s.t[1]^2 + s.t[2]^2)
#'   true_syy = 20 * cos(sqrt(s.t[1]^2 + s.t[2]^2))/sqrt(s.t[1]^2 + s.t[2]^2) -
#'     20 * cos(sqrt(s.t[1]^2 + s.t[2]^2)) *
#'                s.t[2]^2/(s.t[1]^2 + s.t[2]^2)^(3/2) -
#'     20 * sin(sqrt(s.t[1]^2 + s.t[2]^2)) * s.t[2]^2/(s.t[1]^2 + s.t[2]^2)
#'   true_sxx * u.perp[1]^2 + 2 * true_sxy * u.perp[1] * u.perp[2] +
#'                 true_syy * u.perp[2]^2
#'   })
#'   truth[i, 2] = sum(truth.lsegment * (tvec[i]/101))
#' }
#' true.total = colSums(truth); true.total
#' true.avg.total = true.total/sum(tvec); true.avg.total
#' }
#' @author Aritra Halder <aritra.halder@drexel.edu>, \cr
#'  Sudipto Banerjee <sudipto@ucla.edu>
#' @export
spwombling <- function(coords = NULL,
                       curve = NULL,
                       model = NULL,
                       kernel =  c("matern1", "matern2", "gaussian")){
  N = nrow(coords)
  distM = as.matrix(dist(coords))
  tvec = sapply(1:(nrow(curve) - 1), function(x) sqrt(sum((curve[(x + 1),] - curve[x,])^2)))
  umat = as.matrix(t(sapply(1:(nrow(curve) - 1), function(x) (curve[(x + 1),] - curve[x,])))/tvec)

  z = model[, paste0("z[", 1:N, "]")]
  phi = model[, "phi"]
  sigma2 = model[, "sigma2"]
  nmcmc = length(phi)
  if(kernel == "matern1"){
    WM1 = compileNimble(wombling_matern1)
    wmeasure = WM1(coords = coords,
                    curve = curve,
                    dists = distM,
                    tvec = tvec,
                    umat = umat,
                    z = z,
                    phi = phi,
                    sigma2 = sigma2)
  }else if(kernel == "matern2"){
    WM2 = compileNimble(wombling_matern2)
    wmeasure = WM2(coords = coords,
                    curve = curve,
                    dists = distM,
                    tvec = tvec,
                    umat = umat,
                    z = z,
                    phi = phi,
                    sigma2 = sigma2)
  }else{
    WG = compileNimble(wombling_gaussian)
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
    return(list(wm.mcmc = as.mcmc(batch.wmeasure.mcmc),
                estimate.wm = estimate.wm))
  }else{
    return(list(wm.mcmc = as.mcmc(batch.wmeasure.mcmc),
                estimate.wm.1 = estimate.wm.1,
                estimate.wm.2 = estimate.wm.2))
  }
}
