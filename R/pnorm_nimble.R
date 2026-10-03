#' Computes the Cumulative Distribution Function (CDF) for the standard Gaussian
#' probability distribution
#'
#' For internal use only.
#'
#' @param x standard Gaussian quantile
#' @returns A numeric value of the CDF. For internal use only.
#' @keywords pnorm_nimble
#' @importFrom nimble nimbleFunction nimSeq
#' @examples
#' \donttest{
#' #####################
#' # Internal use only #
#' #####################
#' # Example usage inside of nimblewomble::wombling_gaussian()
#' library(nimble)
#' library(nimblewomble)
#'
#' cPnorm_nimble = compileNimble(pnorm_nimble)
#' cPnorm_nimble(1)
#' }
#' @author Aritra Halder <aritra.halder@drexel.edu>, \cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @export
pnorm_nimble <- nimble::nimbleFunction(run = function(x = double(0)) {
  returnType(double(0))

  rule <- seq(0, 1, length.out = 100)

  result <- sum(1 / sqrt(2 * pi) * exp(-(rule[1:100] * x)^2 / 2) * x / 100) +
    0.5
  return(result)
})
