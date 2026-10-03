#' nimblewomble: Bayesian Wombling using nimble
#'
#' Fit Gaussian process models and perform boundary analysis (wombling)
#' using the nimble Bayesian modeling environment.
#'
#' @keywords internal
#' @importFrom stats setNames
"_PACKAGE"

# Silences R CMD check 'no visible binding' notes: .data — rlang pronoun used inside
# aes() (rlang comes via ggplot2).  sig — bare column name referenced in aes() in sp_plot
# / wm_curve_layer.  level — bare column name referenced in aes() in sp_plot.
utils::globalVariables(c(".data", "sig", "level"))
