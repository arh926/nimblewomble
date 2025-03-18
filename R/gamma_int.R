#' Incomplete Gamma Function
#'
#' internal only
#' @param x gamma quantile
#' @param a shape parameter
#' @param b scale parameter
#' @keywords gamma_int
#' @import nimble
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
# cGammint = compileNimble(gamma_int)

