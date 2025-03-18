#' Fit a Gaussian process
#'
#' Fits a Gaussian process with the choice of three kernels. Uses `nimble` to generate posterior samples.
#'
#' @param coords spatial coordinats (supply as a matrix)
#' @param y response
#' @param kernel choice of kernel; must be one of "matern1", "matern2", "gaussian"
#' @param niter number of iterations
#' @param nburn burn-in
#' @keywords gp_fit
#' @import nimble coda stats
#' @export
gp_fit <- function(coords = NULL,
                   y = NULL,
                   kernel = c("matern1", "matern2", "gaussian"),
                   niter = NULL, nburn = NULL){
  N = length(y)
  if(is.null(niter)){
    warning("Defaulting to 1e4 iterations with 5e3 burn-in.")
    niter = 1e4
    nburn = niter/2
  }
  distM = as.matrix(dist(coords))

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
      # Priors #
      sigma2 ~ dinvgamma(shape = 1, rate = 1)
      phi ~ dunif(0, 10)
      tau2 ~ dinvgamma(shape = 2, rate = 1)

      # Initialization #
      mu[1:N] <- zeros[1:N]
      cov[1:N, 1:N] <- materncov1(dists[1:N, 1:N], phi, sigma2, tau2)

      # Likelihood #
      y[1:N] ~ dmnorm(mu[1:N], cov = cov[1:N, 1:N])
    })

  }else if(kernel ==  "matern2"){

    gp_model = quote({
      # Priors #
      sigma2 ~ dinvgamma(shape = 1, rate = 1)
      phi ~ dunif(0, 10)
      tau2 ~ dinvgamma(shape = 2, rate = 1)

      # Initialization #
      mu[1:N] <- zeros[1:N]
      cov[1:N, 1:N] <- materncov2(dists[1:N, 1:N], phi, sigma2, tau2)

      # Likelihood #
      y[1:N] ~ dmnorm(mu[1:N], cov = cov[1:N, 1:N])
    })

  }else{

    gp_model = quote({
      # Priors #
      sigma2 ~ dinvgamma(shape = 1, rate = 1)
      phi ~ dunif(0, 10)
      tau2 ~ dinvgamma(shape = 2, rate = 1)

      # Initialization #
      mu[1:N] <- zeros[1:N]
      cov[1:N, 1:N] <- gaussian(dists[1:N, 1:N], phi, sigma2, tau2)

      # Likelihood #
      y[1:N] ~ dmnorm(mu[1:N], cov = cov[1:N, 1:N])
    })

  }

  model = nimble::nimbleModel(gp_model, constants = constants, data = data, inits = inits)
  cModel = nimble::compileNimble(model)

  conf = nimble::configureMCMC(model)
  MCMC = nimble::buildMCMC(conf)
  cMCMC = nimble::compileNimble(MCMC, project = cModel)

  samples = nimble::runMCMC(cMCMC, niter = niter, nburnin = nburn)
  samples = coda::as.mcmc(samples)

  estimates = rbind(round(quantile(samples[, "tau2"], probs = c(0.5, 0.025, 0.975)), 3),
                    round(quantile(samples[, "sigma2"], probs = c(0.5, 0.025, 0.975)), 3),
                    round(quantile(samples[, "phi"], probs = c(0.5, 0.025, 0.975)), 3))
  rownames(estimates) = c("tau2", "sigma2", "phi")
  return(list(mcmc = samples,
              estimates = estimates))
}
