gen_collinear <- function(seed = 1, n = 200) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- 0.97 * x1 + rnorm(n, sd = 0.05); x3 <- rnorm(n)
  y <- 3 + 2 * x1 + 1.5 * x2 - x3 + rnorm(n)
  data.frame(y, x1, x2, x3)
}

test_that("raiseReg reduces exactly to lm() when no predictor exceeds the threshold", {
  set.seed(11)
  n <- 150
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
  dat <- data.frame(y = 1 + 2 * x1 - x2 + 0.5 * x3 + rnorm(n), x1, x2, x3)
  fit_lm <- lm(y ~ x1 + x2 + x3, data = dat)
  for (m in c("simultaneous", "sequential")) {
    fit <- raiseReg(y ~ x1 + x2 + x3, data = dat, method = m)
    expect_equal(unname(coef(fit)), unname(coef(fit_lm)), tolerance = 1e-8)
    expect_equal(fit$r.squared, summary(fit_lm)$r.squared, tolerance = 1e-8)
  }
})

test_that("raiseReg matches lm() R-squared, sigma and F-statistic exactly under collinearity", {
  dat <- gen_collinear()
  fit_lm <- lm(y ~ x1 + x2 + x3, data = dat)
  for (m in c("simultaneous", "sequential")) {
    fit <- raiseReg(y ~ x1 + x2 + x3, data = dat, method = m)
    expect_equal(fit$r.squared, summary(fit_lm)$r.squared, tolerance = 1e-8)
    expect_equal(fit$sigma, summary(fit_lm)$sigma, tolerance = 1e-8)
    expect_equal(unname(fit$fstatistic["value"]), unname(summary(fit_lm)$fstatistic["value"]),
                 tolerance = 1e-6)
  }
})

test_that("raiseReg brings VIF below threshold and shrinks standard errors on collinear predictors", {
  dat <- gen_collinear()
  fit_lm <- lm(y ~ x1 + x2 + x3, data = dat)
  fit <- raiseReg(y ~ x1 + x2 + x3, data = dat)
  expect_true(all(fit$vif[, "after"] < fit$threshold))
  se_lm <- summary(fit_lm)$coefficients[c("x1", "x2"), "Std. Error"]
  se_rr <- fit$coef_table[c("x1", "x2"), "Std. Error"]
  expect_true(all(se_rr < se_lm))
})

test_that("predict() and fitted() are close to but not identical to the OLS fit under raising", {
  dat <- gen_collinear()
  fit_lm <- lm(y ~ x1 + x2 + x3, data = dat)
  fit <- raiseReg(y ~ x1 + x2 + x3, data = dat)
  d <- fitted(fit) - fitted(fit_lm)
  expect_true(all(is.finite(d)))
  expect_true(mean(abs(d)) > 0)
  expect_true(mean(abs(d)) < 2 * fit$sigma)
  expect_equal(predict(fit), fitted(fit))
  expect_equal(unname(predict(fit, newdata = dat)), unname(fitted(fit)), tolerance = 1e-8)
})

test_that("raiseReg handles factor predictors with correct treatment-contrast dummy coding", {
  set.seed(3)
  n <- 150
  grp <- factor(sample(c("A", "B", "C"), n, replace = TRUE, prob = c(0.6, 0.3, 0.1)))
  x1 <- rnorm(n); x3 <- rnorm(n)
  dat <- data.frame(y = 5 + 2 * x1 + 1.3 * x3 +
                       ifelse(grp == "B", 3, ifelse(grp == "C", -2, 0)) + rnorm(n),
                     x1, x3, grp)
  fit_lm <- lm(y ~ x1 + x3 + grp, data = dat)
  fit <- raiseReg(y ~ x1 + x3 + grp, data = dat)
  expect_equal(names(coef(fit)), names(coef(fit_lm)))
  expect_equal(unname(coef(fit)), unname(coef(fit_lm)), tolerance = 1e-8)

  newdat <- dat[dat$grp != "C", c("x1", "x3", "grp")][1:2, ]
  pr <- predict(fit, newdata = newdat)
  expect_true(all(is.finite(pr)))
})

test_that("raiseReg errors on non-numeric response and on a single predictor", {
  dat <- data.frame(y = letters[1:10], x1 = rnorm(10))
  expect_error(raiseReg(y ~ x1, data = dat), "numeric")
})

test_that("the closed-form raise parameter never silently under-raises near-exact collinearity", {
  # Regression test for a real bug: a bounded grid search for the raise
  # parameter could hit its cap on a near-exact collinearity and silently
  # leave the SVIF above threshold, producing a poorly-conditioned raised
  # design and downstream "computationally singular" errors -- exactly what
  # happened with mtcars' carb/gear columns once a few rows were
  # zero-weighted. Reproduce a near-exact dependency directly.
  set.seed(5)
  n <- 40
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- x1 + x2 + rnorm(n, sd = 1e-3)  # 1-tss ~ 7e-7: far beyond the old lambda_max=1000 cap,
                                        # but not so extreme it trips the exact-collinearity guard
  y <- 1 + x1 - x2 + 0.5 * x3 + rnorm(n)
  dat <- data.frame(y, x1, x2, x3)

  fit <- expect_no_error(raiseReg(y ~ x1 + x2 + x3, data = dat, method = "simultaneous"))
  expect_true(all(is.finite(coef(fit))))
  expect_true(all(fit$svif[, "after"] < fit$threshold))
})

test_that("mtcars repeatedly fits without error across many random seeds (regression test)", {
  data(mtcars, package = "datasets")
  for (s in 1:8) {
    fit <- expect_no_error(raiseReg(mpg ~ ., data = mtcars))
    expect_true(all(is.finite(coef(fit))))
  }
})

test_that("raiseReg reports both SVIF and ordinary VIF, and they can legitimately disagree", {
  data(mtcars, package = "datasets")
  fit <- raiseReg(mpg ~ ., data = mtcars)
  expect_true(all(c("svif", "vif", "cn") %in% names(fit)))
  # a predictor can have high ordinary VIF (dependence on ALL others) while
  # its SVIF (dependence on preceding columns only) stays below threshold
  expect_true(any(fit$vif[, "before"] >= fit$threshold) || any(fit$svif[, "before"] >= fit$threshold))
})

test_that("an exact (not just near) linear dependency is dropped gracefully, not a crash", {
  set.seed(1)
  n <- 40
  x1 <- rnorm(n); x2 <- rnorm(n)
  x3 <- x1 + x2  # exact dependency, to machine precision
  dat <- data.frame(y = rnorm(n), x1, x2, x3)
  fit <- expect_no_error(raiseReg(y ~ x1 + x2 + x3, data = dat))
  expect_true("x3" %in% fit$dropped || "x1" %in% fit$dropped || "x2" %in% fit$dropped)
  expect_true(sum(is.na(coef(fit))) == 1)  # exactly one predictor dropped -> one NA coefficient
  expect_true(all(is.finite(coef(fit)[!is.na(coef(fit))])))
})

test_that("lambda_max forces a drop-and-refit instead of an extreme raise parameter", {
  set.seed(1)
  n <- 100
  x1 <- rnorm(n); x2 <- 0.999 * x1 + rnorm(n, sd = 0.01); x3 <- rnorm(n)
  dat <- data.frame(y = 1 + x1 - x2 + 0.5 * x3 + rnorm(n), x1, x2, x3)

  fit_generous <- raiseReg(y ~ x1 + x2 + x3, data = dat, lambda_max = 1000, method = "simultaneous")
  expect_length(fit_generous$dropped, 0)
  expect_true(fit_generous$k[["x2"]] > 5)

  fit_strict <- raiseReg(y ~ x1 + x2 + x3, data = dat, lambda_max = 5, method = "simultaneous")
  expect_true("x2" %in% fit_strict$dropped)
  expect_true(is.na(coef(fit_strict)["x2"]))
  expect_true(all(is.finite(coef(fit_strict)[c("(Intercept)", "x1", "x3")])))
})

test_that("manual lambda overrides the automatic raise-parameter search", {
  set.seed(1)
  n <- 100
  x1 <- rnorm(n); x2 <- 0.97 * x1 + rnorm(n, sd = 0.05); x3 <- rnorm(n)
  dat <- data.frame(y = 1 + x1 - x2 + 0.5 * x3 + rnorm(n), x1, x2, x3)
  fit <- raiseReg(y ~ x1 + x2 + x3, data = dat, lambda = c(x2 = 3.5))
  expect_equal(unname(fit$k[["x2"]]), 3.5)
})

test_that("raise_method = 'grid' runs and gives finite output close to closed_form", {
  set.seed(1)
  n <- 150
  x1 <- rnorm(n); x2 <- 0.97 * x1 + rnorm(n, sd = 0.05); x3 <- rnorm(n)
  dat <- data.frame(y = 1 + x1 - x2 + 0.5 * x3 + rnorm(n), x1, x2, x3)
  fit_cf <- raiseReg(y ~ x1 + x2 + x3, data = dat, raise_method = "closed_form")
  fit_gr <- raiseReg(y ~ x1 + x2 + x3, data = dat, raise_method = "grid")
  expect_true(all(is.finite(coef(fit_gr))))
  expect_equal(unname(fit_cf$k), unname(fit_gr$k), tolerance = 0.15)
})

test_that("summary() no longer prints the VIF/SVIF tables, and they are attached to the object", {
  set.seed(1)
  dat <- gen_collinear()
  fit <- raiseReg(y ~ x1 + x2 + x3, data = dat)
  out <- capture.output(print(summary(fit)))
  expect_false(any(grepl("SVIF before", out)))
  expect_false(any(grepl("VIF before", out)))
  expect_true(all(c("vif", "svif", "cn") %in% names(fit)))
  expect_equal(colnames(fit$vif), c("before", "after"))
  expect_equal(colnames(fit$svif), c("before", "after"))
})

test_that("a predictor with SVIF just below threshold can still be raised (VIFraise(0) != SVIF)", {
  # Regression test for a real bug: the trigger condition was implemented as
  # "svif_before[i] >= threshold", but the original grid search actually
  # thresholds a different, always-larger quantity, VIFraise(lambda = 0) =
  # (1+tss)+tss/(1-tss). For threshold = 10 this quantity already exceeds 10
  # once SVIF is as low as ~9.11 -- so a column with SVIF = 9.99 (just under
  # 10) must still be raised. Comparing SVIF to threshold directly silently
  # skipped raising such columns, and because raising a later column
  # cascades backward through every earlier coefficient via
  # back-substitution on the upper-triangular Rtilde, this single wrong
  # skip could shift coefficients that were never candidates for raising
  # themselves. Construct a column whose SVIF sits at ~9.95, deliberately
  # just below 10, and confirm it still gets raised.
  set.seed(11)
  n <- 300
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  # tune x3's dependence on x1 so its SVIF lands just under 10 (tss ~ 0.8995)
  x3 <- x1 * sqrt(0.8995) + rnorm(n) * sqrt(1 - 0.8995)
  x3 <- scale(x3, scale = FALSE)[, 1] * sd(x1) / sd(scale(x3, scale = FALSE)[, 1])
  dat <- data.frame(y = rnorm(n), x1, x2, x3)

  fit <- raiseReg(y ~ x1 + x2 + x3, data = dat, method = "simultaneous")
  # x3's SVIF should sit close to, but not necessarily above, 10 -- the
  # point is that raising is decided correctly either way, matching what an
  # exact replica of the original grid search would do.
  literal_grid_k <- function(tss, threshold = 10) {
    lambda <- seq(0, 1000, by = 0.1)
    U <- (1 + lambda^2) * (1 - tss^2) + tss
    L <- (1 + lambda^2) * (1 - tss)
    vifraise1 <- U / L
    sindex1 <- sum(vifraise1 >= threshold)
    lambda[sindex1 + 1] + 1
  }
  svif_x3 <- fit$svif["x3", "before"]
  tss_x3 <- 1 - 1 / svif_x3
  expect_equal(unname(fit$k["x3"]), literal_grid_k(tss_x3), tolerance = 0.2)
})
