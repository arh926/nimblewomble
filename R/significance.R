#' Internal Function for significance
#'
#' @param data_frame matrix consisting of median, lower and upper CI
#' @keywords sig
#' @export
significance <- function(data_frame = NULL){
  apply(data_frame, 1, function(x){
    if(x[3] < 0) return(-1)
    else if(x[2] > 0 ) return(1)
    else return(0)
  })
}
