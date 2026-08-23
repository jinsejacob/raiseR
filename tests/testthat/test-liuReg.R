test_that("liuReg dmm/dopt match the reference liureg package (verified externally)", {
  fit <- liuReg(mpg ~ ., data = mtcars, d = "dmm", clip = FALSE)
  expect_equal(fit$d, -2.121227, tolerance = 1e-4)

  fit_opt <- liuReg(mpg ~ ., data = mtcars, d = "dopt", clip = FALSE)
  expect_equal(fit_opt$d, -0.4315992, tolerance = 1e-4)
})

test_that("liuReg d = 1 reduces to OLS on the centred design", {
  fit <- liuReg(mpg ~ disp + hp + wt, data = mtcars, d = 1)
  lm_fit <- lm(mpg ~ disp + hp + wt, data = mtcars)
  expect_equal(unname(coef(fit)), unname(coef(lm_fit)), tolerance = 1e-8)
})

test_that("liuReg clip = TRUE snaps out-of-range d to [0, 1]", {
  fit <- liuReg(mpg ~ ., data = mtcars, clip = TRUE)
  expect_true(fit$d >= 0 && fit$d <= 1)
})

test_that("liuReg manual d is used as-is", {
  fit <- liuReg(mpg ~ disp + hp + wt, data = mtcars, d = 0.5)
  expect_equal(fit$d, 0.5)
})

test_that("liuReg predict() works with dot formulas and newdata", {
  fit <- liuReg(mpg ~ ., data = mtcars, d = "dmm", clip = FALSE)
  expect_silent(p <- predict(fit, newdata = mtcars[1:3, ]))
  expect_length(p, 3)
})

test_that("robLiu MM produces finite coefficients and predictions", {
  fit <- robLiu(mpg ~ disp + hp + wt, data = mtcars)
  expect_true(all(is.finite(coef(fit))))
  expect_silent(p <- predict(fit, newdata = mtcars[1:2, ]))
  expect_length(p, 2)
})

test_that("robLiu SDO produces finite coefficients for dmm, dopt, and manual d", {
  skip_if_not_installed("mrfDepth")
  fit_dmm <- robLiu(mpg ~ disp + hp + wt, data = mtcars, type = "SDO",
                     seed = 1)
  expect_true(all(is.finite(coef(fit_dmm))))

  fit_dopt <- robLiu(mpg ~ disp + hp + wt, data = mtcars, type = "SDO",
                      d = "dopt", seed = 1)
  expect_true(all(is.finite(coef(fit_dopt))))

  fit_manual <- robLiu(mpg ~ disp + hp + wt, data = mtcars, type = "SDO",
                        d = 0.5, seed = 1)
  expect_equal(fit_manual$d, 0.5)
})
