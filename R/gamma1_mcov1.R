#' Internal function for computing cross-covariance of wombling measures for Matern \eqn{\nu=3/2}
#'
#' @param coords coordinates
#' @param t value of t
#' @param u vector of u
#' @param s0 starting point on curve \eqn{s_0}
#' @param phi posterior sample of \eqn{\phi}
#' @keywords gamma1.mcov1
#' @importFrom nimble nimbleFunction nimDim nimMatrix nimSeq
#' @export
gamma1.mcov1 <- nimble::nimbleFunction(
  run = function(coords = double(2),
                 t = double(0),
                 u = double(1),
                 s0 = double(1),
                 phi = double(0)){
    returnType(double(2))

    n <- dim(coords)[1]

    result <- matrix(nrow = n, ncol = 2, init = FALSE)

    # 100 points for sum
    rule <- seq(0, 1, by = 0.01)
    u.perp <- matrix(nrow = 2, ncol = 1, init = FALSE)
    delta <- matrix(nrow = 100, ncol = 2, init = FALSE)

    u.perp[1, 1] <- u[2]
    u.perp[2, 1] <- -u[1]

    for(i in 1:n){
      delta[1:100, 1] <- -t * u[1] * rule[1:100] + (coords[i, 1] - s0[1])
      delta[1:100, 2] <- -t * u[2] * rule[1:100] + (coords[i, 2] - s0[2])

      # 1-Dimensional Quadrature: Integral as a limit of a sum
      result[i, 1] <- sum(-3 * phi^2 * exp(-sqrt(3) * phi * sqrt(delta[1:100, 1]^2 + delta[1:100, 2]^2)) * (delta[1:100, 1:2] %*% u.perp[1:2, 1])) * (t/101)
      result[i, 2] <- -result[i, 1]
    }
    return(result)
  })


