# pnorm_nimble should match stats::pnorm within numerical tolerance.  The nimble
# implementation is a series expansion used internally by wombling_gaussian, so it needs
# to agree with pnorm at the points typical of that context (roughly [-4, 4], where the
# CDF is non-trivial).

test_that("pnorm_nimble agrees with stats::pnorm in the useful range", {
  # Tail accuracy at |x| >= 3 is ~1e-3 absolute (series approximation); accuracy in the
  # body is much tighter. wombling_gaussian only ever hits |x| < ~3 in practice. Using
  # absolute-diff checks because testthat's tolerance = X is relative, which behaves
  # badly near zero.
  xs_body <- c(-2, -1, -0.5, 0, 0.5, 1, 2)
  got_body <- vapply(xs_body, function(x) pnorm_nimble(x), numeric(1))
  expect_lt(max(abs(got_body - stats::pnorm(xs_body))), 0.005)

  xs_tail <- c(-3, 3)
  got_tail <- vapply(xs_tail, function(x) pnorm_nimble(x), numeric(1))
  expect_lt(max(abs(got_tail - stats::pnorm(xs_tail))), 0.002)
})

test_that("pnorm_nimble is monotone", {
  xs <- seq(-3, 3, length.out = 25)
  ys <- vapply(xs, function(x) pnorm_nimble(x), numeric(1))
  expect_true(all(diff(ys) > 0))
})

test_that("pnorm_nimble at 0 is ~0.5", {
  expect_equal(pnorm_nimble(0), 0.5, tolerance = 1e-04)
})
