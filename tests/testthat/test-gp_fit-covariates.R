# Regression + shape tests for gp_fit's covariate handling. These tests all call gp_fit,
# which triggers nimble MCMC and compileNimble (~30-60s cold).  Marked skip_on_cran to
# keep CRAN checks fast.  NOTE: We explicitly attach nimble here — its model-building
# path looks up internal helpers (e.g., calc_dmnormAltParams) in a scope that fails when
# nimble is only imported (not attached), which is the default when tests only load
# nimblewomble. Symptom without attach: 'R function 'calc_dmnormAltParams' does not
# exist' during nimbleModel().

skip_on_cran()
# nimble is attached via helper-nimble.R

# Tiny problem: N = 20, niter = 200. Keeps each test under ~1 min once the nimble
# toolchain is warm.
set.seed(1)
N <- 20
coords <- matrix(runif(2 * N, -5, 5), ncol = 2)
colnames(coords) <- c("x", "y")
y <- rnorm(N, mean = 2 * coords[, 1], sd = 0.5)

test_that("gp_fit runs with X = NULL (no covariates)", {
  mc <- gp_fit(
    coords = coords,
    y = y,
    kernel = "matern1",
    niter = 200,
    nburn = 100
  )
  expect_type(mc, "list")
  expect_named(mc, c("mcmc", "estimates"))
  expect_true(all(c("phi", "sigma2", "tau2") %in% colnames(mc$mcmc)))
  expect_true(all(c("phi", "sigma2", "tau2") %in% rownames(mc$estimates)))
})

test_that("gp_fit runs with p = 1 covariate (dnorm dispatch)", {
  # Previously the p=1 case blew up because `mu_b[1:1]` collapsed to a scalar inside
  # dmnorm. gp_fit now dispatches to scalar `dnorm(0, sd=100)` when p == 1 and keeps the
  # `beta[1]` name so downstream code (zbeta_samples, spwombling) that references it
  # still works.
  X <- matrix(coords[, 1], ncol = 1)
  mc <- gp_fit(
    coords = coords,
    y = y,
    X = X,
    kernel = "matern1",
    niter = 200,
    nburn = 100
  )
  expect_true("beta[1]" %in% colnames(mc$mcmc))
  expect_true("beta[1]" %in% rownames(mc$estimates))
})

test_that("gp_fit runs with p = 2 covariates", {
  X <- cbind(coords[, 1], coords[, 1]^2)
  mc <- gp_fit(
    coords = coords,
    y = y,
    X = X,
    kernel = "matern1",
    niter = 200,
    nburn = 100
  )
  expect_true(all(c("beta[1]", "beta[2]") %in% colnames(mc$mcmc)))
})

test_that("gp_fit rejects X with a leading intercept column", {
  X_bad <- cbind(1, coords[, 1])
  expect_error(
    gp_fit(
      coords = coords,
      y = y,
      X = X_bad,
      kernel = "matern1",
      niter = 200,
      nburn = 100
    ),
    "without the first column as 1"
  )
})
