# End-to-end integration test for the gaussian (squared-exponential) pipeline.
# Structurally mirrors matern2 (5-component sprates, 2-measure spwombling) but exercises
# the DIFFERENT internal nimbleFunctions: curvatures_gaussian, wombling_gaussian,
# zbeta_gaussian and the pnorm_nimble helper used inside wombling_gaussian's k[1, 1].

skip_on_cran()

set.seed(3)
N <- 20
coords <- matrix(runif(2 * N, -5, 5), ncol = 2)
colnames(coords) <- c("x", "y")
y <- rnorm(N, mean = 0.3 * coords[, 1] - 0.2 * coords[, 2], sd = 0.5)

# See test-pipeline-matern1.R for why niter = 200 exercises the batch-collapse fix in
# spwombling.
mc_sp <- gp_fit(
  coords = coords,
  y = y,
  kernel = "gaussian",
  niter = 200,
  nburn = 100
)
zbeta <- zbeta_samples(
  y = y,
  coords = coords,
  model = mc_sp$mcmc,
  kernel = "gaussian"
)

test_that("sprates on gaussian returns gradients + 3 curvature components", {
  xs <- seq(-4, 4, length.out = 4)
  grid <- as.matrix(expand.grid(xs, xs))
  colnames(grid) <- c("x", "y")

  gr <- sprates(
    coords = coords,
    grid = grid,
    model = zbeta,
    kernel = "gaussian"
  )

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

test_that("spwombling on gaussian returns two wombling measures", {
  t <- seq(0, 1, length.out = 8)
  curve <- cbind(-3 + 6 * t, -1 + 2 * t)

  wm <- spwombling(
    coords = coords,
    curve = curve,
    model = zbeta,
    kernel = "gaussian"
  )

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
