# Frozen reference values for the deterministic (non-mrfDepth) methods.
# These lock in the exact numerical output so that future changes cannot
# silently alter results for the classical estimators. If a change is
# intended to modify these numbers, update the reference values here
# deliberately -- do not loosen the tolerance.

test_that("raiseReg sequential predictions are frozen", {
  fit <- raiseReg(mpg ~ ., mtcars, threshold = 5, method = "sequential")
  expect_equal(unname(round(predict(fit, mtcars[5:8, ]), 5)),
               c(17.22376, 20.42377, 13.91919, 22.55782), tolerance = 1e-4)
})

test_that("raiseReg simultaneous coefficients are frozen", {
  fit <- raiseReg(mpg ~ ., mtcars, threshold = 5, method = "simultaneous")
  expect_equal(unname(round(coef(fit)["wt"], 7)), -2.4951978, tolerance = 1e-5)
  expect_equal(unname(round(fit$cn["after"], 4)), 9.7621, tolerance = 1e-3)
})

test_that("ridgeReg predictions are frozen", {
  fit <- ridgeReg(mpg ~ ., mtcars, k = 2)
  expect_equal(unname(round(predict(fit, mtcars[5:8, ]), 5)),
               c(17.12237, 20.37197, 14.20233, 23.34057), tolerance = 1e-4)
})

test_that("raiseReg diagnostics are frozen", {
  fit <- raiseReg(mpg ~ ., mtcars, threshold = 5, method = "sequential")
  expect_equal(unname(round(cooks.distance(fit)["Merc 230"], 7)), 0.3792206,
               tolerance = 1e-5)
  expect_equal(unname(round(hatvalues(fit)["Merc 230"], 7)), 0.7422870,
               tolerance = 1e-5)
})
