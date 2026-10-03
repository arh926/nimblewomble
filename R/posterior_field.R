#' Joint posterior samples of the fitted surface on a new grid
#'
#' Given posterior samples of GP hyperparameters and the latent field at
#' observed coords (as returned by [zbeta_samples()]), draws `M` joint
#' samples of ŷ(s) = β0 + X_new · β + z(s) at a new grid via the GP
#' conditional:
#'
#'   z(s_new) | z(s_obs), φ, σ² ~ N(K_no K_oo⁻¹ z_obs,
#'                                   K_nn − K_no K_oo⁻¹ K_no').
#'
#' @param coords   n × 2 matrix of observed coordinates.
#' @param grid_new g × 2 matrix of prediction coordinates.
#' @param X_obs    n × p covariate matrix at observed coords (no intercept).
#' @param X_new    g × p covariate matrix at grid coords (no intercept).
#' @param model    posterior samples matrix from [zbeta_samples()]; must
#'   include columns `z[i]`, `phi`, `sigma2`, `beta[0]`, `beta[1..p]`.
#' @param kernel   one of `'matern1'`, `'matern2'`, `'gaussian'`.
#' @param M        number of posterior draws to return.
#' @param jitter_rel small relative jitter added to covariance diagonals
#'   for numerical stability. Default `1e-6`.
#' @returns An `M × g` matrix of joint posterior draws of ŷ on the grid.
#' @keywords posterior_field
#' @importFrom stats dist rnorm
#' @examples
#' \dontrun{
#' Z <- posterior_field(coords, grid_new, X_obs = X, X_new = X_in,
#'                      model = zbeta, kernel = 'matern1', M = 50)
#' }
#' @author Reed Benkendorf
#' @export
posterior_field <- function(
  coords,
  grid_new,
  X_obs,
  X_new,
  model,
  kernel = c("matern1", "matern2", "gaussian"),
  M = 50,
  jitter_rel = 1e-06
) {
  kernel <- match.arg(kernel)
  N <- nrow(coords)
  G <- nrow(grid_new)
  P <- ncol(X_obs)

  z_all <- model[, paste0("z[", seq_len(N), "]")]
  phi_all <- model[, "phi"]
  s2_all <- model[, "sigma2"]
  b0_all <- model[, "beta[0]"]
  beta_all <- model[, paste0("beta[", seq_len(P), "]"), drop = FALSE]
  draws <- round(seq(1, nrow(model), length.out = M))

  # Correlation-only kernel (variance and jitter are added below).
  K_corr <- switch(
    kernel,
    matern1 = function(D, phi) {
      (1 + sqrt(3) * phi * D) *
        exp(
          -sqrt(3) *
            phi *
            D
        )
    },
    matern2 = function(D, phi) {
      (1 + sqrt(5) * phi * D + 5 / 3 * phi^2 * D^2) *
        exp(-sqrt(5) * phi * D)
    },
    gaussian = function(D, phi) exp(-phi^2 * D^2)
  )

  D_oo <- as.matrix(dist(coords))
  D_nn <- as.matrix(dist(grid_new))
  # cross-distances via ||g||² + ||c||² − 2 g·c; pmax guards float roundoff.
  D_no <- sqrt(pmax(
    outer(rowSums(grid_new^2), rowSums(coords^2), "+") -
      2 *
        grid_new %*%
          t(coords),
    0
  ))

  out <- matrix(NA_real_, nrow = M, ncol = G)
  for (i in seq_along(draws)) {
    d <- draws[i]
    phi <- phi_all[d]
    s2 <- s2_all[d]
    jit <- jitter_rel * s2

    K_oo <- s2 * K_corr(D_oo, phi) + diag(jit, N)
    K_no <- s2 * K_corr(D_no, phi)
    K_nn <- s2 * K_corr(D_nn, phi) + diag(jit, G)

    L_oo <- chol(K_oo)
    alpha <- backsolve(L_oo, backsolve(L_oo, z_all[d, ], transpose = TRUE))
    mu_z <- as.numeric(K_no %*% alpha)
    lin <- b0_all[d] + as.numeric(X_new %*% beta_all[d, ])
    W <- backsolve(L_oo, t(K_no), transpose = TRUE)
    Sig_n <- K_nn - crossprod(W)

    L_n <- chol(Sig_n + diag(jit, G))
    out[i, ] <- lin + mu_z + as.numeric(t(L_n) %*% rnorm(G))
  }
  out
}
