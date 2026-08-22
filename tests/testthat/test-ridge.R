test_that("ridgeReg(k = 0) reduces exactly to lm()", {
  set.seed(1)
  n <- 200
  x1 <- rnorm(n); x2 <- 0.97 * x1 + rnorm(n, sd = 0.05); x3 <- rnorm(n)
  dat <- data.frame(y = 3 + 2 * x1 + 1.5 * x2 - x3 + rnorm(n), x1, x2, x3)
  fit_lm <- lm(y ~ x1 + x2 + x3, data = dat)
  fit0 <- ridgeReg(y ~ x1 + x2 + x3, data = dat, k = 0)
  expect_equal(unname(coef(fit0)), unname(coef(fit_lm)), tolerance = 1e-8)
})

test_that("ridgeReg picks a positive biasing parameter under collinearity and shrinks coefficients", {
  set.seed(1)
  n <- 200
  x1 <- rnorm(n); x2 <- 0.97 * x1 + rnorm(n, sd = 0.05); x3 <- rnorm(n)
  dat <- data.frame(y = 3 + 2 * x1 + 1.5 * x2 - x3 + rnorm(n), x1, x2, x3)
  fit <- ridgeReg(y ~ x1 + x2 + x3, data = dat)
  expect_true(fit$k > 0)
  fit_lm <- lm(y ~ x1 + x2 + x3, data = dat)
  ols_norm <- sqrt(sum(coef(fit_lm)[-1]^2))
  ridge_norm <- sqrt(sum(coef(fit)[-1]^2))
  expect_true(ridge_norm < ols_norm)
})

test_that("robRidge downweights outliers and produces finite output", {
  set.seed(1)
  n <- 150
  x1 <- rnorm(n); x2 <- rnorm(n)
  dat <- data.frame(y = 1 + 2 * x1 - x2 + rnorm(n), x1, x2)
  dat$y[1:5] <- dat$y[1:5] + 15
  fit <- robRidge(y ~ x1 + x2, data = dat)
  expect_true(all(is.finite(coef(fit))))
  expect_true(sum(fit$weights == 0) >= 3)
})

test_that("ridgeReg/robRidge R-squared is the squared correlation of actual and fitted values", {
  set.seed(1)
  n <- 150
  x1 <- rnorm(n); x2 <- rnorm(n)
  dat <- data.frame(y = 1 + 2 * x1 - x2 + rnorm(n), x1, x2)
  fit_rg <- ridgeReg(y ~ x1 + x2, data = dat)
  expect_equal(fit_rg$r.squared, cor(fitted(fit_rg), dat$y)^2, tolerance = 1e-10)

  dat$y[1:5] <- dat$y[1:5] + 15
  fit_rrg <- robRidge(y ~ x1 + x2, data = dat)
  expect_equal(fit_rrg$r.squared, cor(fitted(fit_rrg), dat$y)^2, tolerance = 1e-10)
})

test_that("ridgeReg/robRidge report no formal significance test (no t/p-value columns)", {
  set.seed(1)
  n <- 100
  x1 <- rnorm(n); x2 <- rnorm(n)
  dat <- data.frame(y = 1 + x1 - x2 + rnorm(n), x1, x2)
  fit_rg <- ridgeReg(y ~ x1 + x2, data = dat)
  fit_rrg <- robRidge(y ~ x1 + x2, data = dat)
  expect_false(any(c("t value", "Pr(>|t|)") %in% colnames(fit_rg$coef_table)))
  expect_false(any(c("t value", "Pr(>|t|)") %in% colnames(fit_rrg$coef_table)))
})

test_that("predict() on ridgeReg/robRidge matches fitted() on the training data", {
  set.seed(1)
  n <- 100
  x1 <- rnorm(n); x2 <- rnorm(n)
  dat <- data.frame(y = 1 + x1 - x2 + rnorm(n), x1, x2)
  fit_rg <- ridgeReg(y ~ x1 + x2, data = dat)
  fit_rrg <- robRidge(y ~ x1 + x2, data = dat)
  expect_equal(unname(predict(fit_rg, newdata = dat)), unname(fitted(fit_rg)), tolerance = 1e-8)
  expect_equal(unname(predict(fit_rrg, newdata = dat)), unname(fitted(fit_rrg)), tolerance = 1e-8)
})
