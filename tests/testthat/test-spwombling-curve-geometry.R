# The spwombling() refactor replaced sapply-based tvec/umat with diff(curve) +
# rowSums/rownorm. This test reproduces both the old and new forms on the same input to
# lock in numerical equivalence — if a future edit breaks either, this catches it before
# an MCMC round-trip.

test_that("diff-based tvec/umat matches the sapply-based reference", {
  set.seed(1)
  # Randomish 2-D curve, 25 vertices.
  curve <- cbind(cumsum(rnorm(25, 0, 0.4)), cumsum(rnorm(25, 0, 0.4)))

  # New (in-package) form.
  delta_new <- diff(curve)
  tvec_new <- sqrt(rowSums(delta_new^2))
  umat_new <- delta_new / tvec_new

  # Old sapply-based form (verbatim from pre-refactor spwombling).
  tvec_old <- sapply(1:(nrow(curve) - 1), function(x) {
    sqrt(sum((curve[(x + 1), ] - curve[x, ])^2))
  })
  umat_old <- as.matrix(
    t(sapply(1:(nrow(curve) - 1), function(x) {
      (curve[(x + 1), ] - curve[x, ])
    })) /
      tvec_old
  )

  expect_equal(tvec_new, tvec_old, tolerance = 1e-12)
  expect_equal(unname(umat_new), unname(umat_old), tolerance = 1e-12)
})

test_that("unit tangents from spwombling geometry have unit norm", {
  set.seed(2)
  curve <- cbind(cumsum(rnorm(10)), cumsum(rnorm(10)))
  delta <- diff(curve)
  tvec <- sqrt(rowSums(delta^2))
  umat <- delta / tvec
  expect_equal(sqrt(rowSums(umat^2)), rep(1, nrow(umat)), tolerance = 1e-12)
})
