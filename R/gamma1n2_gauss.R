#' Internal function for computing cross-covariance of wombling measures for the squared exponential kernel
#'
#' @param coords coordinates
#' @param t value of t
#' @param u vector of u
#' @param s0 starting point on curve \eqn{s_0}
#' @param phi posterior sample of \eqn{\phi}
#' @keywords gamma1n2.gauss
#' @import nimble
#' @export
gamma1n2.gauss <- nimble::nimbleFunction(
  run = function(coords = double(2),
                 t = double(0),
                 u = double(1),
                 s0 = double(1),
                 phi = double(0)){
    returnType(double(2))

    n <- dim(coords)[1]

    result <- matrix(nrow = n, ncol = 4, init = FALSE)

    # 100 points for sum
    u.perp <- matrix(nrow = 2, ncol = 1, init = FALSE)
    delta <- matrix(nrow = 2, ncol = 1, init = FALSE)
    delta.t <- matrix(nrow = 2, ncol = 1, init = FALSE)

    u.perp[1, 1] <- u[2]
    u.perp[2, 1] <- -u[1]

    for(i in 1:n){
      delta[1, 1] <- (coords[i, 1] - s0[1])
      delta[2, 1] <- (coords[i, 2] - s0[2])

      # 1-Dimensional Quadrature: Integral as a limit of a sum
      result[i, 1] <- -2 * sqrt(pi) * phi * inprod(u.perp, delta) * exp(-phi^2 * inprod(u.perp, delta)^2) * (pnorm_nimble(sqrt(2) * phi * (t + inprod(u, delta))) - pnorm_nimble(sqrt(2) * phi * inprod(u, delta)))
      result[i, 2] <- result[i, 1] * (1 - 2 * phi^2 * inprod(u.perp, delta)^2)
      result[i, 3] <- -result[i, 1]
      result[i, 4] <- result[i, 2]

    }
    return(result)
  })
