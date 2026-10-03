#' Spatial Plot Function
#'
#' @param data_frame data frame consisting of coordinates and data
#' @param sp logical parameter indicating whether to make a spatial plot
#' @param shape if sp = TRUE shape file should be provided (should be an sf object)
#' @param legend.key.height height of legend (defaults to .7)
#' @param legend.key.width width of legend (defaults to .4)
#' @param text.size size of legend text (defaults to 10)
#' @param point.size size of points to be plotted (defaults to 0.7)
#' @param clr.pt color of point to be plotted (defaults to black)
#' @param bound.box bounding box for spatial maps (leave as NULL if not known)
#' @param extend logical parameter indicating whether to extend the interpolation (defaults to TRUE)
#' @param title title of the plot (defaults to NULL)
#' @param palette (optional) color palette
#' @returns A ggplot object.
#' @keywords sp_plot
#' @importFrom methods as
#' @importFrom sf st_coordinates st_transform st_as_sf
#' @importFrom sp proj4string proj4string<- over as.image.SpatialGridDataFrame
#' @importFrom MBA mba.surf
#' @importFrom terra project rast
#' @importFrom ggplot2 ggplot theme_bw aes after_stat labs scale_fill_distiller scale_colour_manual geom_sf geom_point theme element_line element_text unit element_rect geom_raster ggtitle xlim ylim
#' @importFrom ggplot2 guide_legend guide_colourbar
#' @importFrom metR geom_contour2
#' @importFrom ggspatial layer_spatial
#' @examples
#' library(nimblewomble)
#'
#' set.seed(1)
#' # Generated Simulated Data
#' N = 1e2
#' tau = 1
#' coords = matrix(runif(2 * N, -10, 10), ncol = 2)
#' colnames(coords) = c('x', 'y')
#' y = rnorm(N, mean = 20 * sin(sqrt(coords[, 1]^2  + coords[, 2]^2)), sd = tau)
#'
#' sp_ggplot(data_frame = data.frame(coords, z = y))
#' @author Aritra Halder <aritra.halder@drexel.edu>, \cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @export
sp_ggplot <- function(
  data_frame = NULL,
  sp = FALSE,
  shape = NULL,
  legend.key.height = 0.7,
  legend.key.width = 0.4,
  text.size = 10,
  point.size = 0.7,
  clr.pt = "black",
  palette = "Spectral",
  extend = TRUE,
  title = NULL,
  bound.box = NULL
) {
  if (sp && is.null(shape)) {
    stop("Please provide a shape file!")
  }

  # When a 4th column is present it's the CI significance flag from `significance()`
  # (values in {-1, 0, 1}). Map to a factor so ggplot can emit a legend for the point
  # overlay; sig == 0 stays NA and is dropped from the plot via na.translate = FALSE
  # below.
  ci_path <- (ncol(data_frame) == 4)
  if (ci_path) {
    sig_lbl <- factor(
      ifelse(
        data_frame[, 4] == 1,
        "positive",
        ifelse(data_frame[, 4] == -1, "negative", NA_character_)
      ),
      levels = c("positive", "negative")
    )
  } else {
    col.pt <- rep("black", nrow(data_frame))
    color.pt <- rep(clr.pt, nrow(data_frame))
  }

  cnames <- colnames(data_frame)[1:3]
  coords <- data.frame(data_frame[, c(1, 2)])
  if (ci_path) {
    coords$sig <- sig_lbl
  }

  # pch=21 (bordered circle) with `color` mapped to significance and a fixed white fill
  # gives a visible ring against the coloured surface. Fill can't be mapped here because
  # ggplot2 permits only one fill scale per plot and scale_fill_distiller() above owns it
  # (avoids needing ggnewscale).
  build_pt_layer <- function(pts_df) {
    if (ci_path) {
      list(
        geom_point(
          data = pts_df,
          mapping = aes(
            x = .data[[cnames[1]]],
            y = .data[[cnames[2]]],
            colour = sig
          ),
          fill = "white",
          size = point.size,
          stroke = 0.5,
          pch = 21,
          na.rm = TRUE
        ),
        scale_colour_manual(
          name = "CI significance",
          values = c(positive = "#00A651", negative = "#00B4D8"),
          na.translate = FALSE,
          guide = guide_legend(order = 99)
        )
      )
    } else {
      geom_point(
        data = pts_df,
        mapping = aes(x = .data[[cnames[1]]], y = .data[[cnames[2]]]),
        size = point.size,
        color = color.pt,
        fill = col.pt,
        stroke = 0.5,
        pch = 21,
        na.rm = TRUE
      )
    }
  }

  surf <- mba.surf(
    data_frame[, 1:3],
    no.X = 300,
    no.Y = 300,
    h = 5,
    m = 2,
    extend = extend,
    sp = sp
  )$xyz.est

  if (sp) {
    if (!is.null(bound.box)) {
      box_proj <- st_coordinates(st_transform(
        st_as_sf(data.frame(t(bound.box)), coords = c("x", "y"), crs = 4326),
        crs = 5070
      ))
    }

    coords_proj <- st_coordinates(st_transform(
      st_as_sf(coords, coords = cnames[1:2], crs = 4326),
      crs = 5070
    ))
    colnames(coords_proj) <- cnames[1:2]
    coords_proj <- as.data.frame(coords_proj)
    if (ci_path) {
      coords_proj$sig <- sig_lbl
    }

    shape.sp <- as(shape, "Spatial")
    proj4string(surf) <- proj4string(shape.sp)
    surf.tmp <- try(surf[!is.na(over(surf, shape.sp))[, 1], ], silent = TRUE)
    surf <- as.image.SpatialGridDataFrame(surf.tmp)

    surf.rast <- project(rast(surf.tmp), "EPSG:5070")
    df_surf.rast <- as.data.frame(surf.rast, xy = TRUE)
    colnames(df_surf.rast) <- cnames

    gplot <- ggplot() +
      theme_bw() +
      layer_spatial(
        data = surf.rast,
        mapping = aes(fill = "after_stat(band1)")
      ) +
      labs(x = cnames[1], y = cnames[2], fill = "") +
      scale_fill_distiller(
        palette = palette,
        label = function(x) sprintf("%.2f", x),
        na.value = NA,
        direction = -1,
        guide = guide_colourbar(order = 1)
      ) +
      geom_contour2(
        data = df_surf.rast,
        mapping = aes(
          x = .data[[cnames[1]]],
          y = .data[[cnames[2]]],
          z = .data[[cnames[3]]],
          label = round(after_stat(level), 2)
        ),
        linewidth = 0.1,
        label_size = 2.5
      ) +
      geom_sf(data = shape, fill = "transparent") +
      build_pt_layer(coords_proj) +
      theme(
        axis.line = element_line(linewidth = 0.2),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 10),
        legend.key.height = unit(1.5, "cm"),
        legend.key.width = unit(0.6, "cm"),
        legend.text = element_text(size = 20),
        plot.margin = unit(rep(0.15, 4), "cm"),
        panel.background = element_rect(colour = "black", linewidth = 1)
      )
  } else {
    gg.grid <- expand.grid(surf$x, surf$y) # st_coordinates(surf.sf)
    df.gg <- cbind(gg.grid, as.vector(surf$z))
    colnames(df.gg) <- cnames

    gplot <- ggplot() +
      geom_raster(
        data = df.gg,
        mapping = aes(
          x = .data[[cnames[1]]],
          y = .data[[cnames[2]]],
          fill = .data[[cnames[3]]]
        )
      ) +
      labs(x = cnames[1], y = cnames[2], fill = "") +
      scale_fill_distiller(
        palette = palette,
        label = function(x) sprintf("%.2f", x),
        na.value = NA,
        direction = -1,
        guide = guide_colourbar(order = 1)
      ) +
      geom_contour2(
        data = df.gg,
        mapping = aes(
          x = .data[[cnames[1]]],
          y = .data[[cnames[2]]],
          z = .data[[cnames[3]]],
          label = round(after_stat(level), 2)
        ),
        linewidth = 0.1,
        label_size = 2.5
      ) +
      build_pt_layer(coords) +
      theme_bw() +
      theme(
        axis.line = element_line(linewidth = 0.2),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 10),
        legend.key.height = unit(legend.key.height, "cm"),
        legend.key.width = unit(legend.key.width, "cm"),
        legend.text = element_text(size = text.size),
        plot.margin = unit(rep(0.15, 4), "cm"),
        panel.background = element_rect(colour = "black", linewidth = 1)
      )
  }
  if (!is.null(title)) {
    gplot <- gplot + ggtitle(title)
  }
  if (!is.null(bound.box)) {
    gplot <- gplot + xlim(box_proj[, 1]) + ylim(box_proj[, 2])
  }
  return(gplot)
}
