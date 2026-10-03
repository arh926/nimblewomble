#' Fit a Gaussian process
#'
#' Fits a Gaussian process with the choice of three kernels. Uses `nimble` to generate posterior samples.
#'
#' @param coords spatial coordinats (supply as a matrix)
#' @param y response
#' @param X covariates (supply as a matrix without the intercept)
#' @param kernel choice of kernel; must be one of 'matern1', 'matern2', 'gaussian'
#' @param niter number of iterations
#' @param nburn burn-in
#' @returns A list of MCMC samples containing the covariance parameters and the parameter estimates with associated 95% credible intervals.
#' @keywords gp_fit
#' @importFrom nimble nimbleModel compileNimble configureMCMC buildMCMC runMCMC
#' @importFrom coda as.mcmc
#' @importFrom stats dist quantile
#' @examples
#' \donttest{
#' library(nimble)
#' library(nimblewomble)
#'
#' set.seed(1)
#' # Generated Simulated Data
#' N = 1e2
#' tau = 1
#' coords = matrix(runif(2 * N, -10, 10), ncol = 2)
#' colnames(coords) = c('x', 'y')
#' y = rnorm(N, mean = 20 * sin(sqrt(coords[, 1]^2  + coords[, 2]^2)), sd = tau)
#' # Posterior samples for theta
#' mc_sp = gp_fit(coords = coords, y = y, kernel = 'matern2')
#' mc_sp$estimates
#' }
#' @author Aritra Halder <aritra.halder@drexel.edu>, \cr
#'  Sudipto Banerjee <sudipto@ucla.edu>
#' @export
gp_fit <- function(
  coords = NULL,
  y = NULL,
  X = NULL,
  kernel = c("matern1", "matern2", "gaussian"),
  niter = NULL,
  nburn = NULL
) {
  kernel <- match.arg(kernel)
  N <- length(y)
  if (is.null(niter)) {
    warning("Defaulting to 1e4 iterations with 5e3 burn-in.")
    niter <- 10000
    nburn <- niter / 2
  }
  distM <- as.matrix(dist(coords))

  # Kernel dispatch: `kernel_fn` for inits$cov at call time; `kernel_sym` is spliced into
  # the nimble model via bquote so we don't repeat the model body three times.
  kernel_fn <- switch(
    kernel,
    matern1 = materncov1,
    matern2 = materncov2,
    gaussian = gaussian
  )
  kernel_sym <- as.name(switch(
    kernel,
    matern1 = "materncov1",
    matern2 = "materncov2",
    gaussian = "gaussian"
  ))

  inits <- list(sigma2 = 1, phi = 1, tau2 = 0.1)
  inits$cov <- kernel_fn(
    dists = distM,
    sigma2 = inits$sigma2,
    phi = inits$phi,
    tau2 = inits$tau2
  )

  if (is.null(X)) {
    # No covariates: y ~ MVN(0, K). No beta node.
    constants <- list(N = N, dists = distM, zeros = rep(0, N))
    data <- list(y = y)
    gp_model <- bquote({
      mu[1:N] <- zeros[1:N]
      cov[1:N, 1:N] <- .(kernel_sym)(dists[1:N, 1:N], phi, sigma2, tau2)

      sigma2 ~ dinvgamma(shape = 2, rate = 1)
      phi ~ dunif(0, 10)
      tau2 ~ dinvgamma(shape = 2, rate = 0.1)

      y[1:N] ~ dmnorm(mu[1:N], cov = cov[1:N, 1:N])
    })
  } else {
    if (sum(X[, 1] == 1) == N) {
      stop(
        "Please supply X without the first column as 1! ",
        "If no covariates leave X as NULL!"
      )
    }

    p <- ncol(X)
    inits$beta <- rep(1, p)

    if (p == 1) {
      # dmnorm with p=1 fails because the length-1 mean vector collapses to a scalar
      # during model construction (`mu_b[1:1]` -> scalar).  Use scalar dnorm and keep the
      # `beta[1]` name so downstream code that references `beta[1]` (zbeta_samples,
      # spwombling helpers) still works.
      constants <- list(N = N, X = X, dists = distM)
      data <- list(y = y)
      gp_model <- bquote({
        cov[1:N, 1:N] <- .(kernel_sym)(dists[1:N, 1:N], phi, sigma2, tau2)

        beta[1] ~ dnorm(0, sd = 100)
        sigma2 ~ dinvgamma(shape = 2, rate = 1)
        phi ~ dunif(0, 10)
        tau2 ~ dinvgamma(shape = 2, rate = 0.1)

        Xb[1:N] <- X[1:N, 1] * beta[1]
        y[1:N] ~ dmnorm(Xb[1:N], cov = cov[1:N, 1:N])
      })
    } else {
      constants <- list(
        N = N,
        p = p,
        X = X,
        dists = distM,
        mu_b = rep(0, p),
        cov_b = 10000 *
          diag(p)
      )
      data <- list(y = y)
      gp_model <- bquote({
        cov[1:N, 1:N] <- .(kernel_sym)(dists[1:N, 1:N], phi, sigma2, tau2)

        beta[1:p] ~ dmnorm(mu_b[1:p], cov = cov_b[1:p, 1:p])
        sigma2 ~ dinvgamma(shape = 2, rate = 1)
        phi ~ dunif(0, 10)
        tau2 ~ dinvgamma(shape = 2, rate = 0.1)

        Xb[1:N] <- X[1:N, 1:p] %*% beta[1:p]
        y[1:N] ~ dmnorm(Xb[1:N], cov = cov[1:N, 1:N])
      })
    }
  }

  model <- nimbleModel(
    gp_model,
    constants = constants,
    data = data,
    inits = inits
  )
  cModel <- compileNimble(model)
  conf <- configureMCMC(model)
  MCMC <- buildMCMC(conf)
  cMCMC <- compileNimble(MCMC, project = cModel)

  samples <- runMCMC(cMCMC, niter = niter, nburnin = nburn)
  samples <- coda::as.mcmc(samples)

  q <- function(nm) {
    round(quantile(samples[, nm], probs = c(0.5, 0.025, 0.975)), 3)
  }

  estimates <- rbind(tau2 = q("tau2"), sigma2 = q("sigma2"), phi = q("phi"))
  if (!is.null(X)) {
    beta_cols <- paste0("beta[", seq_len(ncol(X)), "]")
    beta_est <- t(apply(
      samples[, beta_cols, drop = FALSE],
      2,
      quantile,
      probs = c(0.5, 0.025, 0.975)
    ))
    rownames(beta_est) <- beta_cols
    estimates <- rbind(estimates, round(beta_est, 3))
  }

  list(mcmc = samples, estimates = estimates)
}
