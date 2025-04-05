#' Determines significance for posterior estimates
#'
#' For internal use only.
#'
#' @param data_frame matrix consisting of median, lower and upper Confidence Interval
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
significance <- function(data_frame = NULL){
  apply(data_frame, 1, function(x){
    if(x[3] < 0) return(-1)
    else if(x[2] > 0 ) return(1)
    else return(0)
  })
}
