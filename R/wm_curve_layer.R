#' ggplot layer: colour a wombling curve by CI significance
#'
#' Returns a list of ggplot2 layers (two `geom_path` passes + a
#' `scale_colour_manual`) that draw the boundary in two passes: a thick
#' black baseline, then a narrower coloured top layer whose colour encodes
#' per-segment CI significance (positive = green, negative = cyan,
#' non-significant segments are transparent and dropped from the legend).
#'
#' `geom_path` colours each segment by its *starting* point's aesthetic,
#' so `sig` (length `nrow(curve) - 1`, one per segment) is padded with an
#' `NA` to line up with `curve`'s rows.
#'
#' Handles both single- and multi-boundary cases with the same call:
#' pass a matrix / data.frame + vector for one boundary, or matching lists
#' for several (each element becomes an independently grouped `geom_path`).
#'
#' @param curve  Either an `n × 2` matrix / data.frame of coordinates
#'   (`x`, `y`), or a list of such objects — one per connected boundary
#'   segment (as produced by `split(as.data.frame(geom(cv)), ...)` on the
#'   output of `terra::as.contour()`).
#' @param sig    CI significance flags, values in `c(-1, 0, 1, NA)`. Either a
#'   length-`(nrow(curve) - 1)` vector (matching a single curve), or a
#'   list of such vectors (matching `curve`'s list). Values are those
#'   returned by [significance()] and stored in
#'   `spwombling()$estimate.wm$sig`.
#' @param values named character vector of colours for `'positive'` and
#'   `'negative'` segments.
#' @param legend_name legend title.
#' @returns A list of ggplot2 layers that can be added to a ggplot with `+`.
#' @keywords wm_curve_layer
#' @importFrom ggplot2 aes geom_path scale_colour_manual guide_legend
#' @examples
#' \dontrun{
#' # Single boundary:
#' p2 + wm_curve_layer(curve, wm$estimate.wm$sig)
#'
#' # Multiple boundary segments, each with its own spwombling() result:
#' p2 + wm_curve_layer(sub_segments,
#'                     lapply(wm_list, function(w) w$estimate.wm$sig))
#' }
#' @author Reed Benkendorf
#' @export
wm_curve_layer <- function(
  curve,
  sig,
  values = c(positive = "#00A651", negative = "#00B4D8"),
  legend_name = "CI significance"
) {
  # Normalise inputs to parallel lists so the multi- and single-boundary branches share
  # one code path.
  if (!is.list(curve) || is.data.frame(curve)) {
    curve <- list(curve)
    sig <- list(sig)
  }
  if (length(curve) != length(sig)) {
    stop("`curve` and `sig` must have matching lengths.")
  }

  curve_df <- do.call(
    rbind,
    lapply(seq_along(curve), function(i) {
      df <- as.data.frame(curve[[i]])
      s <- sig[[i]]
      # Pad sig with NA so length matches nrow(df) (geom_path colours each segment by its
      # starting point's aesthetic).
      df$sig <- factor(
        c(
          ifelse(
            s == 1,
            "positive",
            ifelse(s == -1, "negative", NA_character_)
          ),
          NA
        ),
        levels = c("positive", "negative")
      )
      df$boundary <- i
      df
    })
  )

  list(
    geom_path(
      data = curve_df,
      aes(.data$x, .data$y, group = .data$boundary),
      linewidth = 2,
      inherit.aes = FALSE
    ),
    geom_path(
      data = curve_df,
      aes(.data$x, .data$y, colour = .data$sig, group = .data$boundary),
      linewidth = 1,
      na.rm = TRUE,
      inherit.aes = FALSE
    ),
    scale_colour_manual(
      name = legend_name,
      values = values,
      na.translate = FALSE,
      # order = 99 sorts this below any fill legend when stacked.
      guide = guide_legend(order = 99)
    )
  )
}
