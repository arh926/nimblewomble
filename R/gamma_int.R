#' Incomplete Gamma Function
#'
#' For internal use only. Use integration as a limit of a sum to numerically
#'  compute the incomplete gamma integral
#'
#' @param x gamma quantile
#' @param a shape parameter
#' @param b scale parameter
#' @keywords gamma_int
#' @importFrom nimble nimSeq
#' @examples
#' \dontrun{
#' #####################
#' # Internal use only #
#' #####################
#' # Example usage in nimblewomble::wombling_matern1(...) or,
#' # nimblewomble::wombling_matern2(...)
#' gamma_int(x = 1, a = 1, b = 1/sqrt(3))
#' }
#' @author Aritra Halder <aritra.halder@drexel.edu>, \cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @export
gamma_int <- nimble::nimbleFunction(
  run = function(x = double(0),
                 a = double(0),
                 b = double(0)){
    returnType(double(0))

    rule = seq(0, 1, by = 0.01)

    result = sum((rule[1:100] * x)^(a - 1) * exp(-rule[1:100] * x/b) * x/100)
    return(result)
    })

