test_that("raiseReg diagnostics match an equivalent lm() fit exactly", {
  fit <- raiseReg(mpg ~ ., data = mtcars)
  keep <- setdiff(setdiff(colnames(mtcars), "mpg"), fit$dropped)
  lm_fit <- lm(as.formula(paste("mpg ~", paste(keep, collapse = " + "))), data = mtcars)

  expect_equal(unname(hatvalues(fit)), unname(hatvalues(lm_fit)), tolerance = 1e-8)
  expect_equal(unname(cooks.distance(fit)), unname(cooks.distance(lm_fit)), tolerance = 1e-8)
  expect_equal(unname(covRatio(fit)), unname(covRatio(lm_fit)), tolerance = 1e-8)

  db_fit <- dfbetas(fit); db_lm <- dfbetas(lm_fit)
  expect_equal(unname(db_fit[, colnames(db_lm)]), unname(db_lm), tolerance = 1e-6)
})

test_that("predict() works with dot formulas for every constructor", {
  skip_if_not_installed("mrfDepth")
  fits <- list(
    raiseReg(mpg ~ ., data = mtcars),
    ridgeReg(mpg ~ ., data = mtcars),
    robRaise(mpg ~ ., data = mtcars),
    robRidge(mpg ~ ., data = mtcars),
    liuReg(mpg ~ ., data = mtcars),
    robLiu(mpg ~ ., data = mtcars)
  )
  for (fit in fits) {
    expect_silent(p <- predict(fit, newdata = mtcars[1:2, ]))
    expect_length(p, 2)
  }
})

test_that("bptest and ncvTest work on raiseReg objects (skipped if packages absent)", {
  skip_if_not_installed("lmtest")
  skip_if_not_installed("car")
  fit <- raiseReg(mpg ~ ., data = mtcars)
  expect_silent(bp <- lmtest::bptest(fit))
  expect_silent(nc <- car::ncvTest(fit))
  expect_true(is.numeric(bp$statistic))
  expect_true(is.numeric(nc$ChiSquare))
})
