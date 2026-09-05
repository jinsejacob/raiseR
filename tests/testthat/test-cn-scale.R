test_that("cn() ordinary matches manual eigen-based computation", {
  X <- model.matrix(~ disp + hp + wt + drat, data = mtcars)[, -1]
  R <- cor(X)
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  expect_equal(cn(mpg ~ disp + hp + wt + drat, data = mtcars)$cn,
               sqrt(max(ev) / min(ev)), tolerance = 1e-8)
})

test_that("cn() robust runs and returns valid indices", {
  skip_if_not_installed("mrfDepth")
  out <- cn(mpg ~ disp + hp + wt + drat, data = mtcars, type = "R",
            seed = 1)
  expect_true(out$cn >= 1)
  expect_true(all(out$ci >= 1))
})

test_that("scaleDat() classical matches base::scale()", {
  expect_equal(unclass(scaleDat(mtcars[, 1:3])),
               unclass(base::scale(mtcars[, 1:3])), tolerance = 1e-10,
               ignore_attr = TRUE)
})

test_that("scaleDat() median uses median/MAD", {
  out <- scaleDat(mtcars[, 1:3], type = "median")
  expect_equal(attr(out, "scaled:center")[["mpg"]], median(mtcars$mpg))
  expect_equal(attr(out, "scaled:scale")[["mpg"]], mad(mtcars$mpg))
})

test_that("scaleDat() range produces values in [0, 1]", {
  out <- scaleDat(mtcars[, 1:3], type = "range")
  expect_true(all(out >= 0 & out <= 1))
})

test_that("scaleDat() weighted runs and returns weights attribute", {
  skip_if_not_installed("mrfDepth")
  out <- scaleDat(mtcars[, 1:3], type = "weighted", seed = 1)
  expect_true(!is.null(attr(out, "weights")))
  expect_length(attr(out, "weights"), nrow(mtcars))
})
