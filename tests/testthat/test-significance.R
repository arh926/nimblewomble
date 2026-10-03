# significance() returns +1 / -1 / 0 depending on whether a CI band crosses 0.  The
# function positionally reads x[2] (lower) and x[3] (upper) — these tests lock in that
# contract so a future column-order regression breaks visibly.

test_that("significance returns 1 when CI is strictly positive", {
  df <- data.frame(
    `50%` = 1.5,
    `2.5%` = 0.5,
    `97.5%` = 2.5,
    check.names = FALSE
  )
  expect_equal(significance(df), 1)
})

test_that("significance returns -1 when CI is strictly negative", {
  df <- data.frame(
    `50%` = -1.5,
    `2.5%` = -2.5,
    `97.5%` = -0.5,
    check.names = FALSE
  )
  expect_equal(significance(df), -1)
})

test_that("significance returns 0 when CI straddles zero", {
  df <- data.frame(`50%` = 0.1, `2.5%` = -1, `97.5%` = 1, check.names = FALSE)
  expect_equal(significance(df), 0)
})

test_that("significance is vectorized across rows", {
  df <- data.frame(
    `50%` = c(1.5, -1.5, 0.1),
    `2.5%` = c(0.5, -2.5, -1),
    `97.5%` = c(2.5, -0.5, 1),
    check.names = FALSE
  )
  expect_equal(significance(df), c(1, -1, 0))
})

test_that("boundary case: upper CI == 0 is treated as non-significant", {
  # Rule: x[3] < 0 for -1, x[2] > 0 for +1, else 0. Exact zero on boundary → 0. Locks in
  # current behavior; change with intent.
  df <- data.frame(`50%` = -0.5, `2.5%` = -1, `97.5%` = 0, check.names = FALSE)
  expect_equal(significance(df), 0)
})
