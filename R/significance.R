#' Determines significance for posterior estimates
#'
#' For internal use only.
#'
#' @param data_frame data frame consisting of median, lower and upper confidence interval for estimates
#' @return A data frame consisting median, lower and upper confidence interval for estimates and significance (0 or 1). For internal use only.
#' @keywords significance
#' @examples
#' \dontrun{
#' #####################
#' # Internal use only #
#' #####################
#' # Example usage inside of nimblewomble::spwombling(...)
#' estimate.wm$sig = significance(estimate.wm)
#' }
#' @author Aritra Halder <aritra.halder@drexel.edu>, \cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @export
significance <- function(data_frame = NULL) {
  apply(data_frame, 1, function(x) {
    # NAs propagate up from wombling measures whose posterior variance rounded slightly
    # negative (sqrt(-eps) -> NaN). Return NA so downstream summaries / plotting can skip
    # these rather than error.
    if (anyNA(x)) {
      return(NA_integer_)
    }
    if (x[3] < 0) {
      return(-1)
    } else if (x[2] > 0) {
      return(1)
    } else {
      return(0)
    }
  })
}
