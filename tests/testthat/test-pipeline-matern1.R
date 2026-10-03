# End-to-end integration test for the matern1 pipeline: gp_fit -> zbeta_samples ->
# sprates -> spwombling High-value single test: exercises hundreds of lines across the
# exported user-facing functions plus the internal nimbleFunctions (materncov1,
# zbeta_matern1, gradients_matern1, wombling_matern1). A regression in any of these blows
# up this test, and the shape/invariant checks catch silent geometry bugs (asymmetric
# Sigma.grad, sign flips) that the # PATCH comments in the sources document as historic
# pitfalls.

skip_on_cran()

# ---------- shared tiny fit ------------------------------------------------
set.seed(1)
N <- 20
coords <- matrix(runif(2 * N, -5, 5), ncol = 2)
colnames(coords) <- c("x", "y")
y <- rnorm(N, mean = 2 * coords[, 1], sd = 0.5)

# niter=200 -> 100 post-burn == nbatch=100 -> batch_size=1. This used to trip
# spwombling's `apply(vec, 2, median)` because the vector-collapse dropped dims; the
# refactored version uses `drop = FALSE` and matrixStats::colMedians, so this
# configuration exercises the fix.
mc_sp <- gp_fit(
  coords = coords,
  y = y,
  kernel = "matern1",
  niter = 200,
  nburn = 100
)
zbeta <- zbeta_samples(
  y = y,
  coords = coords,
  model = mc_sp$mcmc,
  kernel = "matern1"
)

# ---------- zbeta_samples --------------------------------------------------

test_that("zbeta_samples appends beta[0] and z[1..N] columns", {
  expect_true("beta[0]" %in% colnames(zbeta))
  expect_true(all(paste0("z[", 1:N, "]") %in% colnames(zbeta)))
  # Hyperparameters are still present (not overwritten).
  expect_true(all(c("phi", "sigma2", "tau2") %in% colnames(zbeta)))
  expect_equal(nrow(zbeta), nrow(mc_sp$mcmc))
})

# ---------- sprates --------------------------------------------------------

test_that("sprates on matern1 returns gradients only (2 components)", {
  xs <- seq(-4, 4, length.out = 4)
  grid <- as.matrix(expand.grid(xs, xs))
  colnames(grid) <- c("x", "y")

  gr <- sprates(coords = coords, grid = grid, model = zbeta, kernel = "matern1")

  # matern1 = gradients only, no curvatures. Return shape.
  expect_named(gr, c("sx.mcmc", "sy.mcmc", "estimate.sx", "estimate.sy"))
  for (nm in c("estimate.sx", "estimate.sy")) {
    est <- gr[[nm]]
    expect_s3_class(est, "data.frame")
    expect_equal(colnames(est), c("50%", "2.5%", "97.5%", "sig"))
    expect_equal(nrow(est), nrow(grid))
    # Quantile ordering — a broken chol/asymmetric Sigma.grad path can return nonsense
    # that violates this.
    expect_true(all(est[["2.5%"]] <= est[["50%"]]))
    expect_true(all(est[["50%"]] <= est[["97.5%"]]))
    # No NAs in the point estimates.
    expect_false(any(is.na(est[["50%"]])))
    # sig is +/-1/0.
    expect_true(all(est$sig %in% c(-1, 0, 1)))
  }
})

# ---------- spwombling -----------------------------------------------------

test_that("spwombling on matern1 returns a single wombling measure", {
  # Simple curve: 8-point path across the domain.
  t <- seq(0, 1, length.out = 8)
  curve <- cbind(-3 + 6 * t, -1 + 2 * t)

  wm <- spwombling(
    coords = coords,
    curve = curve,
    model = zbeta,
    kernel = "matern1"
  )

  # matern1 wombling: one measure per curve segment, so 1 mcmc block + 1 estimate frame.
  expect_named(wm, c("wm.mcmc", "estimate.wm"))
  est <- wm$estimate.wm
  expect_s3_class(est, "data.frame")
  expect_equal(colnames(est), c("50%", "2.5%", "97.5%", "sig"))
  expect_equal(nrow(est), nrow(curve) - 1)
  expect_true(all(est[["2.5%"]] <= est[["50%"]]))
  expect_true(all(est[["50%"]] <= est[["97.5%"]]))
  expect_false(any(is.na(est[["50%"]])))
  expect_true(all(est$sig %in% c(-1, 0, 1)))
})
