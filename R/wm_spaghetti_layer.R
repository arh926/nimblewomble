#' ggplot layer: boundary spaghetti + median line
#'
#' Returns a list of ggplot2 layers that render the posterior on a
#' level-set boundary as a spaghetti overlay: one thin translucent line
#' per posterior draw (from `pb$paths`), plus a thicker opaque line for
#' the median-field contour (from `pb$median_paths`).
#'
#' @param pb            output of [posterior_boundaries()].
#' @param draw_colour   colour for the per-draw spaghetti. Default
#'   `'white'` (reads well on a viridis backdrop).
#' @param draw_alpha    per-draw line alpha. Default `0.5`.
#' @param draw_size     per-draw line width. Default `0.45`.
#' @param median_colour colour for the median-field boundary line.
#'   Default `'black'`.
#' @param median_size   median line width. Default `0.9`.
#' @param draw_legend_colour colour of the 'Posterior Draws' legend key.
#'   Kept separate from `draw_colour` because the default white draws would
#'   be invisible against a legend background. Default `'grey50'`.
#' @param show_legend   logical; draw a legend with 'Median' and
#'   'Posterior Draws' keys? Default `TRUE`. The legend uses the colour
#'   scale, so it will replace any colour scale already on the plot.
#' @returns A list of two `geom_path` layers and a colour scale, addable to
#'   any ggplot with `+`.
#' @keywords wm_spaghetti_layer
#' @importFrom ggplot2 aes geom_path scale_colour_manual guide_legend
#' @examples
#' \dontrun{
#' p2 + wm_spaghetti_layer(pb) +
#'   ggplot2::labs(title = 'Boundary spaghetti')
#' }
#' @author Reed Benkendorf
#' @export
wm_spaghetti_layer <- function(
  pb,
  draw_colour = "white",
  draw_alpha = 0.5,
  draw_size = 0.45,
  median_colour = "black",
  median_size = 0.9,
  draw_legend_colour = "grey50",
  show_legend = TRUE
) {
  labs_ <- c("Median", "Posterior Draws")
  list(
    geom_path(
      data = pb$paths,
      mapping = aes(
        .data$x,
        .data$y,
        group = interaction(.data$draw, .data$seg),
        colour = "Posterior Draws"
      ),
      alpha = draw_alpha,
      linewidth = draw_size,
      inherit.aes = FALSE
    ),
    geom_path(
      data = pb$median_paths,
      mapping = aes(.data$x, .data$y, group = .data$seg, colour = "Median"),
      linewidth = median_size,
      inherit.aes = FALSE
    ),
    scale_colour_manual(
      name = NULL,
      values = c("Median" = median_colour, "Posterior Draws" = draw_colour),
      breaks = labs_,
      guide = if (show_legend) {
        # order = 99 sorts this below any fill legend when stacked.
        guide_legend(
          order = 99,
          override.aes = list(
            colour = c(median_colour, draw_legend_colour),
            linewidth = c(median_size, draw_size),
            alpha = 1
          )
        )
      } else {
        "none"
      }
    )
  )
}
