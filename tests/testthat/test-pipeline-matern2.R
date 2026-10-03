# End-to-end integration test for the matern2 pipeline. Distinct from matern1 because
# sprates/spwombling take separate branches: - sprates(matern2) returns 5 components: sx,
# sy, sxx, sxy, syy (adds curvatures) - spwombling(matern2) returns 2 measures instead of
# 1 and uses different internal nimbleFunctions (curvatures_matern2, wombling_matern2).
# This covers those branches; matern1 covers the other side.

skip_on_cran()

set.seed(2)
N <- 20
coords <- matrix(runif(2 * N, -5, 5), ncol = 2)
colnames(coords) <- c("x", "y")
y <- rnorm(N, mean = coords[, 1]^2 * 0.1, sd = 0.5)

# See test-pipeline-matern1.R for why niter = 200 exercises the batch-collapse fix in
# spwombling.
mc_sp <- gp_fit(
  coords = coords,
  y = y,
  kernel = "matern2",
  niter = 200,
  nburn = 100
)
zbeta <- zbeta_samples(
  y = y,
  coords = coords,
  model = mc_sp$mcmc,
  kernel = "matern2"
)

test_that("sprates on matern2 returns gradients + 3 curvature components", {
  xs <- seq(-4, 4, length.out = 4)
  grid <- as.matrix(expand.grid(xs, xs))
  colnames(grid) <- c("x", "y")

  gr <- sprates(coords = coords, grid = grid, model = zbeta, kernel = "matern2")

  # matern2 exposes the full 2 gradients + 3 curvature components.
  expect_named(
    gr,
    c(
      "sx.mcmc",
      "sy.mcmc",
      "sxx.mcmc",
      "sxy.mcmc",
      "syy.mcmc",
      "estimate.sx",
      "estimate.sy",
      "estimate.sxx",
      "estimate.sxy",
      "estimate.syy"
    )
  )

  for (nm in c(
    "estimate.sx",
    "estimate.sy",
    "estimate.sxx",
    "estimate.sxy",
    "estimate.syy"
  )) {
    est <- gr[[nm]]
    expect_s3_class(est, "data.frame")
    expect_equal(colnames(est), c("50%", "2.5%", "97.5%", "sig"), info = nm)
    expect_equal(nrow(est), nrow(grid), info = nm)
    expect_true(all(est[["2.5%"]] <= est[["50%"]]), info = nm)
    expect_true(all(est[["50%"]] <= est[["97.5%"]]), info = nm)
    expect_false(any(is.na(est[["50%"]])), info = nm)
    expect_true(all(est$sig %in% c(-1, 0, 1)), info = nm)
  }
})

test_that("spwombling on matern2 returns two wombling measures", {
  t <- seq(0, 1, length.out = 8)
  curve <- cbind(-3 + 6 * t, -1 + 2 * t)

  wm <- spwombling(
    coords = coords,
    curve = curve,
    model = zbeta,
    kernel = "matern2"
  )

  # matern2 wombling exposes both first-order (gradient) and second-order (curvature)
  # measures along the curve — two estimate frames.
  expect_named(wm, c("wm.mcmc", "estimate.wm.1", "estimate.wm.2"))
  for (nm in c("estimate.wm.1", "estimate.wm.2")) {
    est <- wm[[nm]]
    expect_s3_class(est, "data.frame")
    expect_equal(colnames(est), c("50%", "2.5%", "97.5%", "sig"), info = nm)
    expect_equal(nrow(est), nrow(curve) - 1, info = nm)
    expect_true(all(est[["2.5%"]] <= est[["50%"]]), info = nm)
    expect_true(all(est[["50%"]] <= est[["97.5%"]]), info = nm)
    expect_false(any(is.na(est[["50%"]])), info = nm)
    expect_true(all(est$sig %in% c(-1, 0, 1)), info = nm)
  }
})
