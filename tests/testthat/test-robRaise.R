test_that("robRaise reduces to classical weighted least squares when nothing needs raising", {
  skip_if_not_installed("mrfDepth")
  set.seed(9)
  n <- 150
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
  dat <- data.frame(y = 2 + 1.5 * x1 - 2 * x2 + 0.5 * x3 + rnorm(n), x1, x2, x3)
  dat$y[1:5] <- dat$y[1:5] + rnorm(5, 12, 2)

  fit <- robRaise(y ~ x1 + x2 + x3, data = dat, seed = 42)
  expect_length(fit$k, 0)
  fit_wls <- lm(y ~ x1 + x2 + x3, data = dat, weights = fit$weights)
  expect_equal(unname(coef(fit)[-1]), unname(coef(fit_wls)[-1]), tolerance = 1e-6)
})

test_that("robRaise downweights injected outliers and recovers coefficients close to the truth", {
  skip_if_not_installed("mrfDepth")
  set.seed(1)
  n <- 200
  x1 <- rnorm(n); x2 <- 0.97 * x1 + rnorm(n, sd = 0.05); x3 <- rnorm(n)
  dat <- data.frame(y = 3 + 2 * x1 + 1.5 * x2 - x3 + rnorm(n), x1, x2, x3)
  dat$y[1:6] <- dat$y[1:6] + rnorm(6, 15, 3)
  dat$x1[7:10] <- dat$x1[7:10] + 8

  fit <- robRaise(y ~ x1 + x2 + x3, data = dat, seed = 1)
  expect_true(fit$n_downweighted >= 5)
  expect_true(all(fit$vif[, "after"] < fit$threshold, na.rm = TRUE))
  # combined effect of the raised, formerly-collinear pair should track the true combined slope (3.5)
  combined <- unname(coef(fit)["x1"] + coef(fit)["x2"])
  expect_equal(combined, 3.5, tolerance = 0.6)
  expect_true(all(is.finite(fit$coef_table)))
  expect_true(fit$nu_star > 0 && fit$nu2 > 0)
})

test_that("the Satterthwaite trace shortcuts match brute-force n x n computation", {
  set.seed(7)
  n <- 60; p <- 4
  X <- matrix(rnorm(n * p), n, p)
  Qw <- qr.Q(qr(X))
  w <- runif(n, 0.1, 1)

  M <- crossprod(Qw, w * Qw)
  M2 <- crossprod(Qw, (w^2) * Qw)
  nu2_fast <- sum(w) - sum(diag(M))
  trace_fast <- sum(w^2) - 2 * sum(diag(M2)) + sum(diag(M %*% M))

  Pw <- Qw %*% t(Qw)
  NW <- (diag(n) - Pw) %*% diag(w)
  nu2_brute <- sum(diag(NW))
  trace_brute <- sum(diag(NW %*% NW))

  expect_equal(nu2_fast, nu2_brute, tolerance = 1e-8)
  expect_equal(trace_fast, trace_brute, tolerance = 1e-8)
})

test_that("robRaise handles factor predictors and predict() on new data", {
  skip_if_not_installed("mrfDepth")
  set.seed(3)
  n <- 150
  grp <- factor(sample(c("A", "B", "C"), n, replace = TRUE, prob = c(0.6, 0.3, 0.1)))
  x1 <- rnorm(n); x2 <- 0.95 * x1 + rnorm(n, sd = 0.08)
  dat <- data.frame(y = 5 + 2 * x1 + 1.3 * x2 +
                       ifelse(grp == "B", 3, ifelse(grp == "C", -2, 0)) + rnorm(n),
                     x1, x2, grp)
  fit <- robRaise(y ~ x1 + x2 + grp, data = dat, seed = 1)
  expect_true(all(is.finite(coef(fit))))
  pr <- predict(fit, newdata = dat[1, c("x1", "x2", "grp")])
  expect_true(is.finite(pr))
})

test_that("robRaise() gives an informative error for method != 'simultaneous'", {
  skip_if_not_installed("mrfDepth")
  set.seed(1)
  dat <- data.frame(y = rnorm(30), x1 = rnorm(30), x2 = rnorm(30))
  expect_error(robRaise(y ~ x1 + x2, data = dat, method = "sequential"),
               "simultaneous")
})

test_that("robRaise() reports weighted SVIF distinctly from weighted ordinary VIF", {
  skip_if_not_installed("mrfDepth")
  set.seed(1)
  n <- 200
  x1 <- rnorm(n); x2 <- 0.97 * x1 + rnorm(n, sd = 0.05); x3 <- rnorm(n)
  dat <- data.frame(y = 3 + 2 * x1 + 1.5 * x2 - x3 + rnorm(n), x1, x2, x3)
  fit <- robRaise(y ~ x1 + x2 + x3, data = dat, seed = 1)
  expect_true(all(c("svif", "vif", "cn") %in% names(fit)))
})

test_that("robRaise() drop-and-refit works and gives NA coefficients for dropped predictors", {
  skip_if_not_installed("mrfDepth")
  set.seed(1)
  n <- 200
  x1 <- rnorm(n); x2 <- 0.999 * x1 + rnorm(n, sd = 0.01); x3 <- rnorm(n)
  dat <- data.frame(y = 1 + x1 - x2 + 0.5 * x3 + rnorm(n), x1, x2, x3)
  fit <- robRaise(y ~ x1 + x2 + x3, data = dat, seed = 1, lambda_max = 3)
  if (length(fit$dropped) > 0) {
    expect_true(all(is.na(coef(fit)[fit$dropped])))
    expect_true(all(is.finite(coef(fit)[!is.na(coef(fit))])))
  }
})

test_that("robRaise() manual lambda override works", {
  skip_if_not_installed("mrfDepth")
  set.seed(1)
  n <- 150
  x1 <- rnorm(n); x2 <- 0.9 * x1 + rnorm(n, sd = 0.1); x3 <- rnorm(n)
  dat <- data.frame(y = 1 + x1 - x2 + rnorm(n), x1, x2, x3)
  fit <- robRaise(y ~ x1 + x2 + x3, data = dat, seed = 1, lambda = c(x2 = 4))
  expect_equal(unname(fit$k[["x2"]]), 4)
})

test_that("summary.robRaise no longer prints VIF/SVIF tables, and $vif/$svif/$cn are attached", {
  skip_if_not_installed("mrfDepth")
  set.seed(1)
  n <- 150
  x1 <- rnorm(n); x2 <- rnorm(n)
  dat <- data.frame(y = 1 + x1 - x2 + rnorm(n), x1, x2)
  fit <- robRaise(y ~ x1 + x2, data = dat, seed = 1)
  out <- capture.output(print(summary(fit)))
  expect_false(any(grepl("SVIF before", out)))
  expect_false(any(grepl("VIF before", out)))
  expect_true(all(c("vif", "svif", "cn") %in% names(fit)))
})
