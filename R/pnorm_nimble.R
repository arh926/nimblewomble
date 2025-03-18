#' Gaussian CDF
#'
#' internal use only
#' @param x standard Gaussian quantile
#' @keywords pnorm_nimble
#' @import nimble
#' @export
pnorm_nimble <- nimble::nimbleFunction(
  run = function(x = double(0)){
    returnType(double(0))

    rule <- seq(0, 1, length.out = 100)

    result <- sum(1/sqrt(2 * pi) * exp(-(rule[1:100] * x)^2/2) * x/100) + 0.5
    return(result)
  })
# cPnorm_nimble = compileNimble(pnorm_nimble)
