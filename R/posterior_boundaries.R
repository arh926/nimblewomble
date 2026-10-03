#' Posterior samples of a level-set boundary on a grid
#'
#' Draws posterior samples of the fitted surface via [posterior_field()]
#' and extracts the level-`c_level` contour from each draw. The set of
#' contours *is* the posterior on the boundary location — sharp where
#' the field is well-identified, diffuse where it is not.
#'
#' The input grid can be either a full rectangular grid (all cells
#' sampled) or a rectangular grid restricted to in-polygon cells via
#' `in_poly` — the latter mirrors the pattern where out-of-polygon
#' cells stay `NA` so `contourLines` terminates cleanly at the study
#' boundary rather than extrapolating.
#'
#' @param coords   n × 2 matrix of observed coordinates.
#' @param grid_uc  data.frame or matrix of all rectangular grid coords
#'   (`x`, `y`) — including cells outside the study polygon.
#' @param X_obs    n × p covariate matrix at observed coords.
#' @param X_new    (# in-polygon cells) × p covariate matrix at the
#'   in-polygon subset of `grid_uc`.
#' @param model    posterior samples matrix from [zbeta_samples()].
#' @param kernel   one of `'matern1'`, `'matern2'`, `'gaussian'`.
#' @param c_level  contour level to extract (on the ŷ scale, e.g.
#'   `log(500)` for the Dutch soil-quality zinc threshold).
#' @param M        number of posterior draws.
#' @param in_poly  optional logical vector of length `nrow(grid_uc)`
#'   marking cells to sample; `NA` elsewhere. Defaults to all `TRUE`.
#' @returns A list with:
#'   \describe{
#'     \item{`paths`}{data.frame `(x, y, draw, seg)` — every contour
#'       vertex from every posterior draw. Feed to [wm_spaghetti_layer()].}
#'     \item{`median_paths`}{data.frame `(x, y, seg)` — contour of the
#'       cell-wise median field.}
#'     \item{`Z_grid`}{`M × nrow(grid_uc)` matrix of joint posterior
#'       draws (in-polygon cells filled; out-of-polygon `NA`). Reuse
#'       for exceedance-probability maps: `colMeans(Z_grid > c_level)`.}
#'     \item{`gx`, `gy`}{unique sorted x / y coords of the rectangular grid.}
#'     \item{`c_level`, `M`}{echoed inputs.}
#'   }
#' @keywords posterior_boundaries
#' @importFrom grDevices contourLines
#' @importFrom stats median
#' @examples
#' \dontrun{
#' pb <- posterior_boundaries(coords, grid_uc,
#'                             X_obs = X, X_new = X_in,
#'                             model = zbeta, kernel = 'matern1',
#'                             c_level = log(500), M = 50,
#'                             in_poly = !is.na(dist_uc))
#' # spaghetti + median line as a ggplot layer:
#' p2 + wm_spaghetti_layer(pb)
#' # per-cell exceedance probability:
#' colMeans(pb$Z_grid > pb$c_level)
#' }
#' @author Reed Benkendorf
#' @export
posterior_boundaries <- function(
  coords,
  grid_uc,
  X_obs,
  X_new,
  model,
  kernel = c("matern1", "matern2", "gaussian"),
  c_level,
  M = 50,
  in_poly = NULL
) {
  kernel <- match.arg(kernel)
  grid_uc <- as.data.frame(grid_uc)
  if (!all(c("x", "y") %in% names(grid_uc))) {
    stop("grid_uc must have columns 'x' and 'y'.")
  }
  if (is.null(in_poly)) {
    in_poly <- rep(TRUE, nrow(grid_uc))
  }

  grid_in <- as.matrix(grid_uc[in_poly, c("x", "y")])
  Z_in <- posterior_field(
    coords,
    grid_in,
    X_obs = X_obs,
    X_new = X_new,
    model = model,
    kernel = kernel,
    M = M
  )
  Z_grid <- matrix(NA_real_, nrow = nrow(Z_in), ncol = nrow(grid_uc))
  Z_grid[, in_poly] <- Z_in

  # contourLines expects z[i, j] = f(gx[i], gy[j]); expand.grid varies x fastest, so a
  # column-major reshape lines up without transposition.
  gx <- sort(unique(grid_uc$x))
  gy <- sort(unique(grid_uc$y))
  contour_df <- function(vec_row, i) {
    z_mat <- matrix(vec_row, nrow = length(gx), ncol = length(gy))
    cl <- contourLines(gx, gy, z_mat, levels = c_level)
    if (!length(cl)) {
      return(NULL)
    }
    do.call(
      rbind,
      lapply(seq_along(cl), function(k) {
        data.frame(x = cl[[k]]$x, y = cl[[k]]$y, draw = i, seg = k)
      })
    )
  }

  paths <- do.call(
    rbind,
    lapply(seq_len(nrow(Z_grid)), function(i) contour_df(Z_grid[i, ], i))
  )

  med_field <- apply(Z_grid, 2, median)
  med_mat <- matrix(med_field, nrow = length(gx), ncol = length(gy))
  med_cl <- contourLines(gx, gy, med_mat, levels = c_level)
  median_paths <- do.call(
    rbind,
    lapply(seq_along(med_cl), function(k) {
      data.frame(x = med_cl[[k]]$x, y = med_cl[[k]]$y, seg = k)
    })
  )

  list(
    paths = paths,
    median_paths = median_paths,
    Z_grid = Z_grid,
    gx = gx,
    gy = gy,
    c_level = c_level,
    M = M
  )
}
