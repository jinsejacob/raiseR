test_that("vif matches car::vif-style computation and is fooled by outliers", {
  skip_if_not_installed("mrfDepth")
  set.seed(1)
  n <- 200
  x1 <- rnorm(n); x2 <- 0.97 * x1 + rnorm(n, sd = 0.05); x3 <- rnorm(n)
  y <- 3 + 2 * x1 + 1.5 * x2 - x3 + rnorm(n)
  dat <- data.frame(y, x1, x2, x3)

  v <- vif(y ~ x1 + x2 + x3, data = dat)
  expect_s3_class(v, "vif")
  expect_true(all(v >= 1))
  expect_true(unname(v["x1"]) > 100)  # strong collinearity correctly flagged

  # classic textbook cross-check: VIF_i == 1 / (1 - R^2 from auxiliary regression)
  aux <- lm(x1 ~ x2 + x3, data = dat)
  r2 <- summary(aux)$r.squared
  expect_equal(unname(v["x1"]), 1 / (1 - r2), tolerance = 1e-8)

  dat_out <- dat
  dat_out$y[1:6] <- dat_out$y[1:6] + rnorm(6, 15, 3)
  dat_out$x1[7:10] <- dat_out$x1[7:10] + 8
  v_out <- vif(y ~ x1 + x2 + x3, data = dat_out)
  rv_out <- rvif(~ x1 + x2 + x3, data = dat_out, seed = 1)
  expect_true(max(v_out) < 10)         # classical VIF is fooled
  expect_true(max(rv_out$rvif) > 10)   # robust VIF is not
})

test_that("vif() and rvif() report a Condition Number", {
  skip_if_not_installed("mrfDepth")
  set.seed(1)
  n <- 200
  x1 <- rnorm(n); x2 <- 0.9 * x1 + rnorm(n, sd = 0.1); x3 <- rnorm(n)
  dat <- data.frame(y = rnorm(n), x1, x2, x3)
  v <- vif(y ~ x1 + x2 + x3, data = dat)
  expect_true(is.numeric(attr(v, "cn")))
  expect_true(attr(v, "cn") > 1)

  rv <- rvif(~ x1 + x2 + x3, data = dat, seed = 1)
  expect_true(is.numeric(rv$cn))
  expect_true(rv$cn > 1)
})
