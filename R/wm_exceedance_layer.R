#' ggplot layer: pointwise exceedance probability + boundary credible zone
#'
#' Returns a list of ggplot2 layers that render, for every grid cell, the
#' posterior probability that the fitted surface exceeds `pb$c_level`
#' (i.e. `colMeans(pb$Z_grid > pb$c_level)`). Overlaid are the `P = 0.5`
#' contour (the median boundary) and the `ci` contours bracketing a
#' *spatial* credible interval on the boundary itself. Both line types are
#' shown as legend keys.
#'
#' The grid is rebuilt from `pb$gx` / `pb$gy` in `expand.grid` order (x
#' varying fastest), matching the column order of `pb$Z_grid`.
#'
#' @param pb            output of [posterior_boundaries()].
#' @param ci            length-2 numeric vector of probability levels for
#'   the credible-zone contours. Default `c(0.05, 0.95)`.
#' @param fill_colours  length-3 character vector of `low`, `mid`, `high`
#'   colours for the diverging fill (midpoint at `P = 0.5`). Ignored when
#'   `fill_scale` is supplied. Default `c('steelblue', 'white', 'firebrick')`.
#' @param fill_scale    optional ggplot2 continuous fill scale (e.g.
#'   `ggplot2::scale_fill_viridis_c()`) to use in place of the default
#'   diverging gradient. If it has no `name`, `fill_name` is used.
#' @param fill_labels   labels for the default fill scale: a function of the
#'   breaks, e.g. `scales::label_number(drop0trailing = TRUE)`. Ignored when
#'   `fill_scale` is supplied. The default drops trailing zeros
#'   (`0, 0.25, 0.5, 0.75, 1`).
#' @param fill_name     fill legend title. Defaults to
#'   `P(z(s) > <c_level>)`.
#' @param median_colour colour of the `P = 0.5` contour. Default `'black'`.
#' @param ci_colour     colour of the credible-zone contours. Default
#'   `'black'`.
#' @param linewidth     contour line width. Default `0.4`.
#' @param median_linetype,ci_linetype line types for the median and
#'   credible-zone contours. Defaults `'solid'` and `'dashed'`.
#' @param interpolate   passed to [ggplot2::geom_raster()]. Default `TRUE`.
#' @param show_legend   logical; draw the median / CI line legend? Default
#'   `TRUE`.
#' @returns A list of ggplot2 layers and scales, addable to any ggplot
#'   (including an empty `ggplot()`) with `+`.
#' @keywords wm_exceedance_layer
#' @importFrom ggplot2 aes geom_raster geom_contour scale_fill_gradient2
#' @importFrom ggplot2 scale_linetype_manual guide_legend guide_colourbar
#' @importFrom stats setNames
#' @examples
#' \dontrun{
#' ggplot2::ggplot() + wm_exceedance_layer(pb) +
#'   ggplot2::labs(title = 'Boundary exceedance probability')
#'
#' # custom palette and line colours
#' ggplot2::ggplot() +
#'   wm_exceedance_layer(pb,
#'     fill_scale = ggplot2::scale_fill_viridis_c(limits = c(0, 1)),
#'     median_colour = 'white', ci_colour = 'grey80')
#' }
#' @author Reed Benkendorf
#' @export
wm_exceedance_layer <- function(
  pb,
  ci = c(0.05, 0.95),
  fill_colours = c("steelblue", "white", "firebrick"),
  fill_scale = NULL,
  fill_labels = function(b) formatC(b, format = "g"),
  fill_name = sprintf("P(z(s) > %.2f)", pb$c_level),
  median_colour = "black",
  ci_colour = "black",
  linewidth = 0.4,
  median_linetype = "solid",
  ci_linetype = "dashed",
  interpolate = TRUE,
  show_legend = TRUE
) {
  if (length(ci) != 2) {
    stop("`ci` must be a length-2 vector of probabilities.")
  }
  excee_df <- expand.grid(x = pb$gx, y = pb$gy)
  excee_df$p <- colMeans(pb$Z_grid > pb$c_level)

  med_lab <- "Median"
  ci_lab <- sprintf("%g%% CI", round(100 * abs(diff(ci)), 1))
  labs_ <- c(med_lab, ci_lab)

  if (is.null(fill_scale)) {
    fill_scale <- scale_fill_gradient2(
      low = fill_colours[1],
      mid = fill_colours[2],
      high = fill_colours[3],
      midpoint = 0.5,
      limits = c(0, 1),
      na.value = NA,
      labels = fill_labels,
      name = fill_name,
      guide = guide_colourbar(order = 1)
    )
  } else if (inherits(fill_scale$name, "waiver")) {
    fill_scale$name <- fill_name
  }

  list(
    geom_raster(
      data = excee_df,
      mapping = aes(.data$x, .data$y, fill = .data$p),
      interpolate = interpolate,
      na.rm = TRUE,
      inherit.aes = FALSE
    ),
    fill_scale,
    geom_contour(
      data = excee_df,
      mapping = aes(.data$x, .data$y, z = .data$p, linetype = med_lab),
      breaks = 0.5,
      colour = median_colour,
      linewidth = linewidth,
      na.rm = TRUE,
      inherit.aes = FALSE
    ),
    geom_contour(
      data = excee_df,
      mapping = aes(.data$x, .data$y, z = .data$p, linetype = ci_lab),
      breaks = ci,
      colour = ci_colour,
      linewidth = linewidth,
      na.rm = TRUE,
      inherit.aes = FALSE
    ),
    scale_linetype_manual(
      name = NULL,
      values = stats::setNames(c(median_linetype, ci_linetype), labs_),
      breaks = labs_,
      guide = if (show_legend) {
        # order = 99 sorts this below the fill legend when stacked.
        guide_legend(
          order = 99,
          override.aes = list(colour = c(median_colour, ci_colour))
        )
      } else {
        "none"
      }
    )
  )
}
