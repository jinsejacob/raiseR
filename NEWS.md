# raiseR 0.1.0

Initial release.

## Collinearity diagnostics

* `vif()` -- Variance Inflation Factor, with `type = "O"` (ordinary,
  correlation-matrix based) and `type = "R"` (robust, the RVIF of Jacob &
  Varadharajan, 2024). Accepts a formula, a fitted model object from this
  package (or an `lm`), or a numeric predictor matrix. `rvif()` is a
  convenience alias for `vif(type = "R")`.
* `cn()` -- Condition Number and full vector of Condition Indices, with
  ordinary (`type = "O"`) and robust weighted (`type = "R"`) versions.

## Regression estimators

* `raiseReg()` -- ordinary Raise Regression, with a sequential
  single-variable method (the default) and a simultaneous SVIF/QR method.
  Because raising leaves the column space of the design unchanged, its
  fitted values, residual standard error, R-squared and F-statistic are
  numerically identical to OLS, and its coefficient `t`-tests stay exactly
  valid.
* `robRaise()` -- Robust Raise Regression with exact finite-sample
  inference (sandwich covariance, effective degrees of freedom,
  Satterthwaite-Welch correction), downweighting outliers via Stahel-Donoho
  projection outlyingness and Tukey's biweight function.
* `ridgeReg()` / `robRidge()` -- ordinary and robust Ridge Regression.
  `robRidge()` offers `type = "MM"` (IRLS with an MM seed) and
  `type = "SDO"` (Stahel-Donoho weighting, then ordinary ridge).
* `liuReg()` / `robLiu()` -- ordinary and robust Liu Regression. The
  biasing parameter `d` defaults to the MSE-optimal `dopt`, with `dmm` and
  a manual value also available, and out-of-range estimates clipped to
  `[0, 1]` by default. `robLiu()` offers `type = "MM"` (following the
  MM-Liu estimator of Filzmoser & Kurnaz, with a closed-form biasing
  parameter) and `type = "SDO"`.

## Utilities and methods

* `scaleDat()` -- scale a dataset by one of four conventions: classical
  (mean/sd), robust weighted (Stahel-Donoho/Tukey biweight), median/MADN,
  or min-max to `[0, 1]`.
* `print()`, `summary()`, `coef()`, `fitted()`, `residuals()`,
  `predict(newdata = )` and `plot()` methods for every fitted object,
  modelled on the corresponding `lm` methods.
* For `raiseReg()` (the exact, unbiased fit), the standard influence
  diagnostics `hatvalues()`, `cooks.distance()`, `dfbetas()` and
  `covRatio()`, plus `lmtest::bptest()` and `car::ncvTest()`, are supported
  and return values numerically identical to an equivalent `lm()` fit.
  `robRaise()` supports `bptest()` and `ncvTest()`.

## Dependencies

* The robust methods (`robRaise()`, `robRidge(type = "SDO")`,
  `robLiu(type = "SDO")`, `vif(type = "R")`, `cn(type = "R")`,
  `scaleDat(type = "weighted")`) require the `mrfDepth` package, which
  provides the projection outlyingness measure used in the published
  methodology.
