#' Ordinary Ridge Regression
#'
#' Fits Ridge Regression (Hoerl and Kennard, 1970) on a mean-centred design,
#' with the biasing parameter chosen as
#' \eqn{k = p\hat\sigma^2 / \hat\beta'\hat\beta}, where \eqn{\hat\beta} and
#' \eqn{\hat\sigma^2} come from the ordinary least squares fit.
#'
#' @param formula a two-sided formula.
#' @param data a data frame.
#' @param k optional fixed biasing parameter; if \code{NULL} (default) it is
#'   estimated from the data as described above.
#' @return an object of class \code{"ridgeReg"}.
#' @examples
#' fit <- ridgeReg(mpg ~ disp + hp + wt + drat, data = mtcars)
#' summary(fit)
#' @export
ridgeReg <- function(formula, data, k = NULL) {
  cl <- match.call()
  formula <- stats::as.formula(formula)
  dsg <- .raise_design(formula, data)
  X <- dsg$X; yc <- dsg$y
  n <- dsg$n; p <- dsg$p
  XtX <- crossprod(X)

  beta_ols <- as.numeric(solve(XtX, crossprod(X, yc)))
  df_ols <- n - p - 1
  sigma2_ols <- sum((yc - X %*% beta_ols)^2) / df_ols

  if (is.null(k)) {
    k <- p * sigma2_ols / sum(beta_ols^2)
  }

  A <- XtX + diag(k, p)
  A_inv <- solve(A)
  beta <- as.numeric(A_inv %*% crossprod(X, yc))
  names(beta) <- colnames(X)

  fitted_c <- as.numeric(X %*% beta)
  residuals <- yc - fitted_c
  df_residual <- n - p - 1
  sigma2 <- sum(residuals^2) / df_residual
  # R^2 for a biased/shrunk estimator is reported as the squared correlation
  # between actual and fitted values (the standard convention for ridge-type
  # estimators), not as 1 - SSE/TSS: the latter is only guaranteed to lie in
  # [0, 1] for OLS, and can behave oddly once shrinkage is introduced.
  r_squared <- stats::cor(yc, fitted_c)^2

  Vhat_slopes <- sigma2 * (A_inv %*% XtX %*% A_inv)
  se_slopes <- sqrt(diag(Vhat_slopes))
  intercept <- dsg$y_mean - sum(dsg$means * beta)
  var_intercept <- sigma2 / n + as.numeric(t(dsg$means) %*% Vhat_slopes %*% dsg$means)
  se_intercept <- sqrt(var_intercept)

  coefficients <- c(`(Intercept)` = intercept, beta)
  se <- c(`(Intercept)` = se_intercept, se_slopes)
  # No t/p-value columns: ridge shrinkage has no exact (or even asymptotically
  # standard) sampling distribution the way OLS or the raise estimator do, so
  # a formal significance test would overstate what these standard errors
  # actually support. Only point estimates and an approximate SE are shown.
  coef_table <- cbind(Estimate = coefficients, `Approx. SE` = se)

  structure(
    list(coefficients = coefficients, coef_table = coef_table,
         fitted.values = fitted_c + dsg$y_mean, residuals = residuals,
         k = k, sigma = sqrt(sigma2), df.residual = df_residual,
         r.squared = r_squared, call = cl, formula = dsg$formula, xlevels = dsg$xlevels, n = n, p = p,
         x = cbind(`(Intercept)` = rep(1, n), dsg$X_raw), y = yc + dsg$y_mean),
    class = "ridgeReg")
}

#' @export
print.ridgeReg <- function(x, ...) {
  cat("Ordinary Ridge Regression (k = ", signif(x$k, 5), ")\n\n", sep = "")
  cat("Call:\n"); print(x$call)
  cat("\nCoefficients:\n")
  print(round(x$coefficients, 5))
  invisible(x)
}

#' @export
summary.ridgeReg <- function(object, ...) structure(object, class = c("summary.ridgeReg", class(object)))

#' @export
print.summary.ridgeReg <- function(x, ...) {
  cat("Ordinary Ridge Regression\n\n")
  cat("Call:\n"); print(x$call); cat("\n")
  cat(sprintf("Biasing parameter k: %.5f\n\n", x$k))
  cat("Residuals:\n"); print(summary(as.numeric(x$residuals)))
  cat("\nCoefficients (point estimates only -- see note below):\n")
  print(round(x$coef_table, 5))
  cat(sprintf("\nResidual standard error: %.4f on %d degrees of freedom\n",
              x$sigma, x$df.residual))
  cat(sprintf("R-squared (squared correlation of actual and fitted values): %.4f\n",
              x$r.squared))
  invisible(x)
}

#' @export
coef.ridgeReg <- function(object, ...) object$coefficients

#' @export
fitted.ridgeReg <- function(object, ...) object$fitted.values

#' @export
residuals.ridgeReg <- function(object, ...) object$residuals

#' @export
predict.ridgeReg <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(object$fitted.values)
  X <- .new_design_matrix(object$formula, newdata, object$xlevels)
  b <- object$coefficients
  as.numeric(b["(Intercept)"] + X[, names(b)[-1], drop = FALSE] %*% b[-1])
}

#' Diagnostic plots for a ridgeReg object
#' @param x a \code{ridgeReg} object.
#' @param which subset of 1:2.
#' @param ... unused.
#' @return Invisibly returns the fitted \code{ridgeReg} object `x`.
#'   Called for its side effect of drawing diagnostic plots.
#' @export
plot.ridgeReg <- function(x, which = 1:2, ...) {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  if (length(which) > 1) graphics::par(mfrow = c(1, length(which)))
  fv <- x$fitted.values; rs <- x$residuals; std_rs <- rs / x$sigma
  if (1 %in% which) {
    graphics::plot(fv, rs, xlab = "Fitted values", ylab = "Residuals",
                   main = "Residuals vs Fitted", pch = 20, col = "steelblue4")
    graphics::abline(h = 0, lty = 2, col = "grey40")
  }
  if (2 %in% which) {
    stats::qqnorm(std_rs, main = "Normal Q-Q", pch = 20, col = "steelblue4")
    stats::qqline(std_rs, col = "firebrick3", lwd = 2)
  }
  invisible(x)
}
