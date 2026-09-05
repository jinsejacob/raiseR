# raiseR

Ordinary and Robust Raise, Ridge, and Liu Regression, with Condition Number
and Variance Inflation Factor diagnostics, for R.

`raiseR` implements the **Raise Regression**: an alternative to Ridge Regression
for combating multicollinearity that, unlike Ridge, leaves ordinary
least-squares inference intact. By construction, a raised design matrix
reproduces the OLS fitted values, residual standard error, R-squared and
F-statistic *exactly* -- raising only reallocates the fitted signal among
collinear predictors, stabilising their individual coefficients and standard
errors, without ever changing what the model actually explains.

The package also provides **robust** counterparts of every method for data
contaminated by outliers, a **Robust Variance Inflation Factor** and robust
**Condition Number** that resist being fooled by outliers the way their
classical versions can be, and ordinary and robust **Ridge** and **Liu**
regression.

## Installation

```r
# install.packages("remotes")
remotes::install_github("jinsejacob/raiseR")
```

`raiseR` depends on `mrfDepth` (for the projection outlyingness measure that
powers all of the robust methods) and `MASS` (for the MM-estimator used by
`robRidge(type = "MM")` and `robLiu(type = "MM")`). Both install
automatically from CRAN.

## Why raise, not ridge?

Ridge Regression shrinks coefficients by adding a penalty to the normal
equations, but this breaks the usual `t`- and F-testing machinery: exact,
finite-sample inference for a ridge estimate is not generally available.
Raise Regression instead re-expresses a collinear predictor as a linear
combination of itself and the part of it left unexplained by the other
predictors, which:

* reduces its Variance Inflation Factor below a chosen threshold (10 by
  default),
* **never changes the column space of the design matrix**, so the fitted
  values, residual sum of squares, R-squared and F-statistic are unchanged,
  and
* keeps the ordinary `t`-test machinery exactly valid for the raised
  coefficients.

## Ordinary VIF and Raise Regression

```r
library(raiseR)

set.seed(1)
n <- 200
x1 <- rnorm(n)
x2 <- 0.97 * x1 + rnorm(n, sd = 0.05)   # strongly collinear with x1
x3 <- rnorm(n)
y  <- 3 + 2 * x1 + 1.5 * x2 - x3 + rnorm(n)
dat <- data.frame(y, x1, x2, x3)

vif(y ~ x1 + x2 + x3, data = dat)       # ordinary VIF (type = "O")

fit <- raiseReg(y ~ x1 + x2 + x3, data = dat)   # sequential (default)
summary(fit)
```

R-squared, sigma and the F-statistic of the raised fit are numerically
identical to `lm(y ~ x1 + x2 + x3, data = dat)`; only the coefficient split
between the two collinear predictors, and the precision with which each is
estimated, changes. `method = "simultaneous"` fits the QR/SVIF-based
simultaneous raise strategy instead of the one-variable-at-a-time default.

## When outliers are also present

A handful of outliers can mask real collinearity from the classical VIF and
distort an ordinary raise fit. The robust methods downweight observations by
Tukey's biweight function applied to their Stahel-Donoho projection
outlyingness before doing anything else:

```r
dat_out <- dat
dat_out$y[1:6]   <- dat_out$y[1:6] + rnorm(6, 15, 3)   # y-outliers
dat_out$x1[7:10] <- dat_out$x1[7:10] + 8               # X-outliers

vif(y ~ x1 + x2 + x3, data = dat_out)                  # may be fooled
vif(y ~ x1 + x2 + x3, data = dat_out, type = "R", seed = 1)   # robust: not fooled

fit_r <- robRaise(y ~ x1 + x2 + x3, data = dat_out, seed = 1)
summary(fit_r)
```

## Ridge and Liu Regression

```r
ridgeReg(mpg ~ ., data = mtcars)               # ordinary ridge, auto k
robRidge(mpg ~ ., data = mtcars, type = "MM")  # robust ridge (MM)
robRidge(mpg ~ ., data = mtcars, type = "SDO", seed = 1)  # robust ridge (SDO)

liuReg(mpg ~ ., data = mtcars)                 # ordinary Liu, d = dopt
robLiu(mpg ~ ., data = mtcars, type = "MM")    # robust Liu (MM)
```

## Diagnostics and utilities

```r
cn(mpg ~ disp + hp + wt, data = mtcars)        # Condition Number + indices
cn(mpg ~ disp + hp + wt, data = mtcars, type = "R")   # robust version

scaleDat(mtcars, type = "median")              # median/MADN scaling
scaleDat(mtcars, type = "range")               # min-max to [0, 1]
```

## Fitted-object methods

Every fitted object supports `print()`, `summary()`, `coef()`, `fitted()`,
`residuals()`, `predict(newdata = )` and `plot()`. For `raiseReg()` -- the
exact, unbiased fit -- the standard influence diagnostics `hatvalues()`,
`cooks.distance()`, `dfbetas()` and `covRatio()`, plus `lmtest::bptest()`
and `car::ncvTest()`, are also available and return values numerically
identical to an equivalent `lm()` fit.

## A note on inference

For `raiseReg()`, the residual standard error, R-squared, F-statistic and
the standard errors underlying the coefficient table are computed from the
*raised* design and are therefore numerically identical to their OLS
counterparts -- this is the entire point of the raise regression. `fitted()`,
`residuals()` and `predict()`, which must apply to new data, apply the raise
coefficients to the *original, unraised* predictors instead, and so differ
very slightly from the training fit implied by the reported R-squared.

## References

Jacob, J. and Varadharajan, R. (2023). Raise Estimation: An Alternative
Approach in the Presence of Problematic Multicollinearity. *Mathematics and
Statistics*, 11(1), 51-64. <https://doi.org/10.13189/ms.2023.110106>

Jacob, J. and Varadharajan, R. (2022). Simultaneous raise regression: a
novel approach to combating collinearity in linear regression models.
*Quality & Quantity*, 57, 4365-4386.
<https://doi.org/10.1007/s11135-022-01557-9>

Jacob, J. and Varadharajan, R. (2024). Robust Variance Inflation Factor: A
Promising Approach for Collinearity Diagnostics in the Presence of Outliers.
*Sankhya B*, 86(2), 845-871.
<https://doi.org/10.1007/s13571-024-00342-y>

Jacob, J. (2025). Enhancing Linear Regression with Raise Techniques to
Effectively Tackle the Multicollinearity Problem. PhD thesis, SRM Institute
of Science and Technology, Kattankulathur, India.

Hoerl, A. E. and Kennard, R. W. (1970). Ridge regression: Biased estimation
for nonorthogonal problems. *Technometrics*, 12(1), 55-67.
<https://doi.org/10.1080/00401706.1970.10488634>

Liu, K. (1993). A new class of biased estimate in linear regression.
*Communications in Statistics - Theory and Methods*, 22(2), 393-402.
<https://doi.org/10.1080/03610929308831027>

Yohai, V. J. (1987). High breakdown-point and high efficiency robust
estimates for regression. *The Annals of Statistics*, 15(2), 642-656.
<https://doi.org/10.1214/aos/1176350366>

Filzmoser, P. and Kurnaz, F. S. (2018). A robust Liu regression estimator.
*Communications in Statistics - Simulation and Computation*, 47(2), 432-443.
<https://doi.org/10.1080/03610918.2016.1271889>

Kan, B., Alpu, O. and Yazici, B. (2013). Robust ridge and robust Liu
estimator for regression based on the LTS estimator. *Journal of Applied
Statistics*, 40(3), 644-655.
<https://doi.org/10.1080/02664763.2012.750285>

## License

GPL (>= 3)
