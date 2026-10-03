# gamma_int(x, a, b) = int_0^x s^(a-1) * exp(-s/b) ds Closed forms: a = 1: b * (1 -
# exp(-x/b)) a = 2: b^2 * (1 - exp(-x/b)) - b*x*exp(-x/b) The nimble implementation is a
# Riemann-sum approximation, so we allow a looser tolerance than the closed form warrants
# — this test locks in the approximation's calibration rather than requiring bit-exact
# equality.

# gamma_int is a Riemann-sum approximation; observed error is up to ~3% of the integral
# value, growing with x/b. Absolute tolerance of 0.02 locks in that scale — a regression
# to O(0.1) or a sign flip will still fail.
CLOSED_FORM_TOL <- 0.02

test_that("gamma_int matches closed form for a = 1", {
  for (params in list(
    list(x = 1, b = 1),
    list(x = 2, b = 1),
    list(x = 0.5, b = 2),
    list(x = 3, b = 0.5)
  )) {
    expected <- params$b * (1 - exp(-params$x / params$b))
    got <- gamma_int(x = params$x, a = 1, b = params$b)
    expect_lt(
      abs(got - expected),
      CLOSED_FORM_TOL,
      label = sprintf(
        "|gamma_int(a=1, x=%g, b=%g) - expected|",
        params$x,
        params$b
      )
    )
  }
})

test_that("gamma_int matches closed form for a = 2", {
  for (params in list(
    list(x = 1, b = 1),
    list(x = 2, b = 1),
    list(x = 0.5, b = 2),
    list(x = 3, b = 0.5)
  )) {
    expected <- params$b^2 *
      (1 - exp(-params$x / params$b)) -
      params$b * params$x * exp(-params$x / params$b)
    got <- gamma_int(x = params$x, a = 2, b = params$b)
    expect_lt(
      abs(got - expected),
      CLOSED_FORM_TOL,
      label = sprintf(
        "|gamma_int(a=2, x=%g, b=%g) - expected|",
        params$x,
        params$b
      )
    )
  }
})

test_that("gamma_int is monotone in x for fixed a, b", {
  b <- 1
  for (a in c(1, 2)) {
    xs <- seq(0.1, 5, by = 0.5)
    ys <- vapply(xs, function(x) gamma_int(x = x, a = a, b = b), numeric(1))
    expect_true(
      all(diff(ys) > 0),
      info = sprintf("gamma_int not monotone for a=%d", a)
    )
  }
})
