# Kernel-value tests for the three covariance functions. These are pure numerical checks
# — no MCMC — so they're fast and safe on CRAN.  Each kernel is called uncompiled (as a
# nimbleFunction directly), which is slower per call but avoids the compileNimble
# round-trip; we only need a handful of evaluations per test.

# Small distance matrix for repeated use.
D <- matrix(c(0, 1, 2, 1, 0, 1.5, 2, 1.5, 0), 3, 3, byrow = TRUE)

test_that("materncov1 matches closed-form Matern nu=3/2", {
  phi <- 0.7
  sigma2 <- 2
  tau2 <- 0.1
  got <- materncov1(dists = D, phi = phi, sigma2 = sigma2, tau2 = tau2)
  expected_off <- sigma2 * (1 + sqrt(3) * phi * D) * exp(-sqrt(3) * phi * D)
  expected <- expected_off + tau2 * diag(3)

  expect_equal(got, expected, tolerance = 1e-10)
  expect_equal(diag(got), rep(sigma2 + tau2, 3), tolerance = 1e-10)
  expect_true(isSymmetric(got))
})

test_that("materncov2 matches closed-form Matern nu=5/2", {
  phi <- 0.5
  sigma2 <- 1.5
  tau2 <- 0.05
  got <- materncov2(dists = D, phi = phi, sigma2 = sigma2, tau2 = tau2)
  expected_off <- sigma2 *
    (1 + sqrt(5) * phi * D + 5 / 3 * phi^2 * D^2) *
    exp(-sqrt(5) * phi * D)
  expected <- expected_off + tau2 * diag(3)

  expect_equal(got, expected, tolerance = 1e-10)
  expect_equal(diag(got), rep(sigma2 + tau2, 3), tolerance = 1e-10)
  expect_true(isSymmetric(got))
})

test_that("gaussian kernel matches closed-form squared-exponential", {
  phi <- 0.3
  sigma2 <- 1
  tau2 <- 0
  got <- gaussian(dists = D, phi = phi, sigma2 = sigma2, tau2 = tau2)
  expected <- sigma2 * exp(-phi^2 * D^2) + tau2 * diag(3)

  expect_equal(got, expected, tolerance = 1e-10)
  expect_true(isSymmetric(got))
  # tau2 = 0 case: diagonal should be exactly sigma2.
  expect_equal(diag(got), rep(sigma2, 3), tolerance = 1e-10)
})

test_that("kernels return zero-distance = sigma2 + tau2 on the diagonal", {
  D1 <- matrix(0, 2, 2)
  sigma2 <- 3
  tau2 <- 0.5
  for (kfun in list(materncov1, materncov2, gaussian)) {
    K <- kfun(dists = D1, phi = 1, sigma2 = sigma2, tau2 = tau2)
    expect_equal(diag(K), rep(sigma2 + tau2, 2), tolerance = 1e-10)
  }
})
