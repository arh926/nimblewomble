#' Computes the Cumulative Distribution Function for the standard Gaussian
#' probability distribution
#'
#' For internal use only.
#'
#' @param x standard Gaussian quantile
#' @keywords pnorm_nimble
#' @importFrom nimble nimbleFunction nimSeq
#' @examples
#' \dontrun{
#' #####################
#' # Internal use only #
#' #####################
#' # Example usage inside of nimblewomble::wombling_gaussian()
#' pnorm_nimble(sqrt(2) * phi[i] * tvec[j]) - 1)
#' }
#' @author Aritra Halder <aritra.halder@drexel.edu>, \cr
#' Sudipto Banerjee <sudipto@ucla.edu>
#' @export
pnorm_nimble <- nimble::nimbleFunction(
  run = function(x = double(0)){
    returnType(double(0))

    rule <- seq(0, 1, length.out = 100)

    result <- sum(1/sqrt(2 * pi) * exp(-(rule[1:100] * x)^2/2) * x/100) + 0.5
    return(result)
  })
