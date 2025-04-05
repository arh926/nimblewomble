#' Fit a Gaussian process
#'
#' Fits a Gaussian process with the choice of three kernels. Uses `nimble` to generate posterior samples.
#'
#' @param coords spatial coordinats (supply as a matrix)
#' @param y response
#' @param X covariates (supply as a matrix without the intercept)
#' @param kernel choice of kernel; must be one of "matern1", "matern2", "gaussian"
#' @param niter number of iterations
#' @param nburn burn-in
#' @keywords gp_fit
#' @importFrom nimble nimbleModel compileNimble configureMCMC buildMCMC runMCMC
#' @importFrom coda as.mcmc
#' @importFrom stats dist quantile
#' @examples
#' \dontrun{
#' require(nimble)
#' require(nimblewomble)
#'
#' set.seed(1)
#' # Generated Simulated Data
#' N = 1e2
#' tau = 1
#' coords = matrix(runif(2 * N, -10, 10), ncol = 2)
#' colnames(coords) = c("x", "y")
#' y = rnorm(N, mean = 20 * sin(sqrt(coords[, 1]^2  + coords[, 2]^2)), sd = tau)
#' # Posterior samples for theta
#' mc_sp = gp_fit(coords = coords, y = y, kernel = "matern2")
#' mc_sp$estimates
#' }
#' @author Aritra Halder <aritra.halder@drexel.edu>, \cr
#'  Sudipto Banerjee <sudipto@ucla.edu>
#' @export
gp_fit <- function(coords = NULL,
                   y = NULL, X = NULL,
                   kernel = c("matern1", "matern2", "gaussian"),
                   niter = NULL, nburn = NULL){
  N = length(y)
  if(is.null(niter)){
    warning("Defaulting to 1e4 iterations with 5e3 burn-in.")
    niter = 1e4
    nburn = niter/2
  }
  distM = as.matrix(dist(coords))

  if(sum(X[,1] == 1) == N) stop("Please supply X without the first column as 1! If no covariates leave X as NULL!")
  if(is.null(X)){
    p = ncol(X)

    constants = list(N = N,
                     dists = distM,
                     zeros = rep(0, N))
    data = list(y = y)
    inits = list(sigma2 = 1, phi = 1, tau2 = 0.1)

    if(kernel ==  "matern1") inits$cov = materncov1(dists = distM, sigma2 =  inits$sigma2, phi = inits$phi, tau2 = inits$tau2)
    else if(kernel ==  "matern2") inits$cov = materncov2(dists = distM, sigma2 =  inits$sigma2, phi = inits$phi, tau2 = inits$tau2)
    else inits$cov = gaussian(dists = distM, sigma2 =  inits$sigma2, phi = inits$phi, tau2 = inits$tau2)

    if(kernel ==  "matern1"){
      gp_model = quote({
        # Initialization #
        mu[1:N] <- zeros[1:N]
        cov[1:N, 1:N] <- materncov1(dists[1:N, 1:N], phi, sigma2, tau2)

        # Priors #
        sigma2 ~ dinvgamma(shape = 2, rate = 1)
        phi ~ dunif(0, 10)
        tau2 ~ dinvgamma(shape = 2, rate = 0.1)

        # Likelihood #
        y[1:N] ~ dmnorm(mu[1:N], cov = cov[1:N, 1:N])
      })

    }else if(kernel ==  "matern2"){

      gp_model = quote({
        # Initialization #
        mu[1:N] <- zeros[1:N]
        cov[1:N, 1:N] <- materncov2(dists[1:N, 1:N], phi, sigma2, tau2)

        # Priors #
        sigma2 ~ dinvgamma(shape = 2, rate = 1)
        phi ~ dunif(0, 10)
        tau2 ~ dinvgamma(shape = 2, rate = 0.1)

        # Likelihood #
        y[1:N] ~ dmnorm(mu[1:N], cov = cov[1:N, 1:N])
      })

    }else{

      gp_model = quote({
        # Initialization #
        mu[1:N] <- zeros[1:N]
        cov[1:N, 1:N] <- gaussian(dists[1:N, 1:N], phi, sigma2, tau2)

        # Priors #
        sigma2 ~ dinvgamma(shape = 2, rate = 1)
        phi ~ dunif(0, 10)
        tau2 ~ dinvgamma(shape = 2, rate = 0.1)

        # Likelihood #
        y[1:N] ~ dmnorm(mu[1:N], cov = cov[1:N, 1:N])
      })

    }

    model = nimbleModel(gp_model, constants = constants, data = data, inits = inits)
    cModel = compileNimble(model)

    conf = configureMCMC(model)
    MCMC = buildMCMC(conf)
    cMCMC = compileNimble(MCMC, project = cModel)

    samples = runMCMC(cMCMC, niter = niter, nburnin = nburn)
    samples = coda::as.mcmc(samples)

    estimates = rbind(round(quantile(samples[, "tau2"], probs = c(0.5, 0.025, 0.975)), 3),
                      round(quantile(samples[, "sigma2"], probs = c(0.5, 0.025, 0.975)), 3),
                      round(quantile(samples[, "phi"], probs = c(0.5, 0.025, 0.975)), 3))
    rownames(estimates) = c("tau2", "sigma2", "phi")
    return(list(mcmc = samples,
                estimates = estimates))
  }else{
    p = ncol(X)

    constants = list(N = N, p = p, X = X,
                     dists = distM,
                     zeros = rep(0, p),
                     Ip = diag(p))
    data = list(y = y)
    inits = list(sigma2 = 1, phi = 1, tau2 = 0.1, beta = rep(1, p))

    if(kernel ==  "matern1") inits$cov = materncov1(dists = distM, sigma2 =  inits$sigma2, phi = inits$phi, tau2 = inits$tau2)
    else if(kernel ==  "matern2") inits$cov = materncov2(dists = distM, sigma2 =  inits$sigma2, phi = inits$phi, tau2 = inits$tau2)
    else inits$cov = gaussian(dists = distM, sigma2 =  inits$sigma2, phi = inits$phi, tau2 = inits$tau2)

    if(kernel ==  "matern1"){
      gp_model = quote({
        # Initialization #
        mu_b[1:p] <- zeros[1:p]
        cov_b[1:p, 1:p] <- 1e4 * Ip[1:p, 1:p]
        cov[1:N, 1:N] <- materncov1(dists[1:N, 1:N], phi, sigma2, tau2)

        # Priors #
        beta[1:p] ~ dmnorm(mu_b[1:p], cov = cov_b[1:p, 1:p])
        sigma2 ~ dinvgamma(shape = 2, rate = 1)
        phi ~ dunif(0, 10)
        tau2 ~ dinvgamma(shape = 2, rate = 0.1)

        # Likelihood #
        Xb[1:N] <- X[1:N, 1:p] %*% beta[1:p]
        y[1:N] ~ dmnorm(Xb[1:N], cov = cov[1:N, 1:N])
      })

    }else if(kernel ==  "matern2"){

      gp_model = quote({
        # Initialization #
        mu_b[1:p] <- zeros[1:p]
        cov_b[1:p, 1:p] <- 1e4 * Ip[1:p, 1:p]
        cov[1:N, 1:N] <- materncov2(dists[1:N, 1:N], phi, sigma2, tau2)

        # Priors #
        beta[1:p] ~ dmnorm(mu_b[1:p], cov = cov_b[1:p, 1:p])
        sigma2 ~ dinvgamma(shape = 2, rate = 1)
        phi ~ dunif(0, 10)
        tau2 ~ dinvgamma(shape = 2, rate = 0.1)

        # Likelihood #
        Xb[1:N] <- X[1:N, 1:p] %*% beta[1:p]
        y[1:N] ~ dmnorm(Xb[1:N], cov = cov[1:N, 1:N])
      })

    }else{

      gp_model = quote({
        # Initialization #
        mu_b[1:p] <- zeros[1:p]
        cov_b[1:p, 1:p] <- 1e4 * Ip[1:p, 1:p]
        cov[1:N, 1:N] <- gaussian(dists[1:N, 1:N], phi, sigma2, tau2)

        # Priors #
        beta[1:p] ~ dmnorm(mu_b[1:p], cov = cov_b[1:p, 1:p])
        sigma2 ~ dinvgamma(shape = 2, rate = 1)
        phi ~ dunif(0, 10)
        tau2 ~ dinvgamma(shape = 2, rate = 0.1)

        # Likelihood #
        Xb[1:N] <- X[1:N, 1:p] %*% beta[1:p]
        y[1:N] ~ dmnorm(Xb[1:N], cov = cov[1:N, 1:N])
      })

    }

    model = nimbleModel(gp_model, constants = constants, data = data, inits = inits)
    cModel = compileNimble(model)

    conf = configureMCMC(model)
    MCMC = buildMCMC(conf)
    cMCMC = compileNimble(MCMC, project = cModel)

    samples = runMCMC(cMCMC, niter = niter, nburnin = nburn)
    samples = as.mcmc(samples)

    estimates = rbind(round(quantile(samples[, "tau2"], probs = c(0.5, 0.025, 0.975)), 3),
                      round(quantile(samples[, "sigma2"], probs = c(0.5, 0.025, 0.975)), 3),
                      round(quantile(samples[, "phi"], probs = c(0.5, 0.025, 0.975)), 3),
                      round(t(apply(samples[, paste0("beta[", 1:p,"]")], 2, quantile, probs = c(0.5, 0.025, 0.975))), 3))
    rownames(estimates) = c("tau2", "sigma2", "phi", paste0("beta[", 1:p,"]"))
    return(list(mcmc = samples,
                estimates = estimates))
  }
}
