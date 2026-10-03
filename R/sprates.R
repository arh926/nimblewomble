#' Posterior samples for rates of change
#'
#' @param coords coordinates
#' @param grid grid for sampling the rates of change
#' @param model posterior samples of \eqn{Z(s)}, \eqn{\phi}, \eqn{\sigma^2}
#' @param kernel choice of kernel; must be one of 'matern1', 'matern2', 'gaussian'
#' @returns A list containing MCMC samples for gradients and curvatures and the associated estimates and 95% confidence intervals
#' @keywords sprates
#' @importFrom nimble compileNimble
#' @importFrom coda as.mcmc
#' @importFrom stats quantile
#' @importFrom matrixStats colMedians colQuantiles
#' @examples
#' \donttest{
#' library(nimble)
#' library(nimblewomble)
#'
#' set.seed(1)
#' # Generated Simulated Data
#' N = 1e2
#' tau = 1
#' coords = matrix(runif(2 * N, -10, 10), ncol = 2); colnames(coords) = c('x', 'y')
#' y = rnorm(N, mean = 20 * sin(sqrt(coords[, 1]^2  + coords[, 2]^2)), sd = tau)
#'
#' # Create equally spaced grid of points
#' xsplit = ysplit = seq(-10, 10, by = 1)[-c(1, 21)]
#' grid = as.matrix(expand.grid(xsplit, ysplit), ncol = 2)
#' colnames(grid) = c('x', 'y')
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
#' ##########################
#' # Fit a Gaussian Process #
#' ##########################
#' # Posterior samples for theta
#' mc_sp = gp_fit(coords = coords, y = y, kernel = 'matern2')
#' # Posterior samples for Z(s) and beta
#' model = zbeta_samples(y = y, coords = coords,
#'                       model = mc_sp$mcmc,
#'                       kernel = 'matern2')
#'
#' ###################
#' # Rates of Change #
#' ###################
#' gradients = sprates(grid = grid,
#'                     coords = coords,
#'                     model = model,
#'                     kernel = 'matern2')
#' p8 = sp_ggplot(data_frame = data.frame(grid,
#'                 z = gradients$estimate.sx[,'50%'],
#'                sig = gradients$estimate.sx$sig))
#' p9 = sp_ggplot(data_frame = data.frame(grid,
#'                z = gradients$estimate.sy[,'50%'],
#'                sig = gradients$estimate.sy$sig))
#' p10 = sp_ggplot(data_frame = data.frame(grid,
#'                  z = gradients$estimate.sxx[,'50%'],
#'                 sig = gradients$estimate.sxx$sig))
#' p11 = sp_ggplot(data_frame = data.frame(grid,
#'                  z = gradients$estimate.sxy[,'50%'],
#'                  sig = gradients$estimate.sxy$sig))
#' p12 = sp_ggplot(data_frame = data.frame(grid,
#'                 z = gradients$estimate.syy[,'50%'],
#'                 sig = gradients$estimate.syy$sig))
#' }
#' @author Aritra Halder <aritra.halder@drexel.edu>, \cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @export
sprates <- function(
  coords = NULL,
  grid = NULL,
  model = NULL,
  kernel = c("matern1", "matern2", "gaussian")
) {
  kernel <- match.arg(kernel)

  N <- nrow(coords)
  G <- nrow(grid)
  distM <- as.matrix(dist(coords))

  # Distances and signed component displacements from each grid point to every
  # observation — fed to the compiled derivative kernel below.
  dist.2 <- sapply(seq_len(G), function(j) {
    apply(coords, 1, function(x) sqrt(sum((x - grid[j, ])^2)))
  })
  dist.3 <- do.call(
    rbind,
    sapply(
      seq_len(G),
      function(j) {
        apply(coords, 1, function(x) {
          (x -
            grid[j, ])
        })
      },
      simplify = FALSE
    )
  )

  z <- model[, paste0("z[", 1:N, "]")]
  phi <- model[, "phi"]
  sigma2 <- model[, "sigma2"]
  nmcmc <- length(phi)

  # Kernel-specific compile + evaluation. compileNimble regenerates C++ every call
  # (nimble does not cache), so this ~30-60s step dominates wall-clock.
  compiled <- switch(
    kernel,
    matern1 = compileNimble(gradients_matern1),
    matern2 = compileNimble(curvatures_matern2),
    gaussian = compileNimble(curvatures_gaussian)
  )
  raw <- compiled(
    dists.1 = distM,
    dists.2 = dist.2,
    dists.3 = dist.3,
    z = z,
    phi = phi,
    sigma2 = sigma2
  )

  # matern1 returns gradients only; matern2/gaussian add three curvatures.  Components
  # are interleaved column-wise in `raw` — one array() reshape gets us a per-component
  # view without copying.
  comp_names <- if (kernel == "matern1") {
    c("sx", "sy")
  } else {
    c("sx", "sy", "sxx", "sxy", "syy")
  }
  D <- length(comp_names)
  arr <- array(raw, dim = c(nmcmc, D, G))

  # Batch medians reduce MCMC autocorrelation before quantiling.
  nbatch <- 100
  split_samples <- split(
    seq_len(nmcmc),
    ceiling(seq_len(nmcmc) / (nmcmc / nbatch))
  )

  # Per-component pipeline: batched draws -> quantile summary + significance.  colMedians
  # / colQuantiles from matrixStats are ~10-20x faster than apply(., 2, median|quantile)
  # and avoid the apply()->list->rbind detour.
  summarize <- function(mat) {
    batches <- matrix(0, nrow = nbatch, ncol = ncol(mat))
    for (b in seq_along(split_samples)) {
      batches[b, ] <- matrixStats::colMedians(
        mat[split_samples[[b]], , drop = FALSE],
        na.rm = TRUE
      )
    }
    q <- round(
      matrixStats::colQuantiles(batches, probs = c(0.5, 0.025, 0.975)),
      3
    )
    df <- as.data.frame(q)
    colnames(df) <- c("50%", "2.5%", "97.5%")
    df$sig <- significance(data_frame = df)
    list(mcmc = as.mcmc(batches), estimate = df)
  }

  results <- lapply(seq_len(D), function(d) summarize(arr[, d, ]))
  names(results) <- comp_names

  # Preserve the original return shape: sx.mcmc, sy.mcmc, ..., estimate.sx, ...
  c(
    setNames(lapply(results, `[[`, "mcmc"), paste0(comp_names, ".mcmc")),
    setNames(lapply(results, `[[`, "estimate"), paste0("estimate.", comp_names))
  )
}
