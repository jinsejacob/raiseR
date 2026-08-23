#' Ordinary Liu Regression
#'
#' Fits Liu Regression (Liu, 1993) on a mean-centred design:
#' \deqn{\hat\beta_{Liu} = (X'X + I)^{-1}(X'X + dI)\hat\beta_{OLS}}
#' where \eqn{0 < d < 1} is the biasing parameter.
#'
#' Three ways of choosing \code{d} are available, computed from the
#' canonical (eigenvector) form of the OLS fit:
#' \describe{
#'   \item{\code{"dopt"} (default)}{The exact MSE-minimizing \eqn{d}, derived
#'     from the canonical-form bias-variance tradeoff directly. Numerically
#'     better-behaved than \code{"dmm"} when predictors are close to exactly
#'     collinear (see below).}
#'   \item{\code{"dmm"}}{The original Liu (1993) plug-in formula,
#'     \eqn{d = 1 - \hat\sigma^2 \sum 1/[\lambda_i(\lambda_i+1)] \big/ \sum
#'     \hat\alpha_i^2/(\lambda_i+1)^2}, where \eqn{\lambda_i} are the
#'     eigenvalues of \eqn{X'X} and \eqn{\hat\alpha_i} are the OLS
#'     coefficients rotated into the eigenvector basis. This formula divides
#'     by eigenvalues directly, so on near-exactly singular designs it can
#'     return an extreme, meaningless value; a warning is issued in that
#'     regime and \code{"dopt"} is suggested instead.}
#'   \item{a numeric value}{Used directly as \code{d}, bypassing both
#'     formulas.}
#' }
#' Both \code{"dmm"} and \code{"dopt"} are plug-in estimates built from
#' \eqn{\hat\sigma^2} and \eqn{\hat\beta_{OLS}}; in finite samples,
#' especially under multicollinearity or a low signal-to-noise ratio, they
#' can fall outside the \eqn{(0, 1)} range the underlying theory assumes.
#' \code{clip = TRUE} (the default) snaps an out-of-range estimate to the
#' nearest boundary (0 or 1) rather than using it as computed -- turn this
#' off only if you want to inspect the raw, unclipped formula output, since
#' an unclipped out-of-range d can make predictions considerably worse.
#'
#' @param formula a two-sided formula.
#' @param data a data frame.
#' @param d either \code{"dopt"} (default), \code{"dmm"}, or a fixed
#'   numeric value in principle taken in \eqn{(0, 1)} (values outside this
#'   range are allowed, with a warning when \code{clip = FALSE}).
#' @param clip if \code{TRUE} (default), an automatically-computed \code{d}
#'   ("dopt" or "dmm") outside \eqn{[0, 1]} is snapped to the nearer
#'   boundary. Ignored when \code{d} is supplied numerically. Set to
#'   \code{FALSE} to inspect the raw, unclipped formula output instead.
#' @return an object of class \code{"liuReg"}.
#' @examples
#' fit <- liuReg(mpg ~ disp + hp + wt + drat, data = mtcars)
#' summary(fit)
#' fit_mm <- liuReg(mpg ~ disp + hp + wt + drat, data = mtcars, d = "dmm")
#' fit_manual <- liuReg(mpg ~ disp + hp + wt + drat, data = mtcars, d = 0.5)
#' @export
liuReg <- function(formula, data, d = c("dopt", "dmm"), clip = TRUE) {
  cl <- match.call()
  formula <- stats::as.formula(formula)
  dsg <- .raise_design(formula, data)
  X <- dsg$X; yc <- dsg$y
  n <- dsg$n; p <- dsg$p
  XtX <- crossprod(X)

  beta_ols <- as.numeric(solve(XtX, crossprod(X, yc)))
  # sigma2_ols feeding the d formulas uses n - p (not n - p - 1) to exactly
  # match Liu (1993) / the reference liureg package's own convention,
  # verified numerically against it; this is distinct from the *fitted Liu
  # model's own* reported df.residual/sigma further below, which follows
  # this package's usual n - p - 1 convention (see ?ridgeReg).
  df_ols <- n - p
  sigma2_ols <- sum((yc - X %*% beta_ols)^2) / df_ols

  d_method <- NULL
  if (is.character(d)) {
    d_method <- match.arg(d)
    eig <- eigen(XtX, symmetric = TRUE)
    lambda <- eig$values
    alpha <- as.numeric(t(eig$vectors) %*% beta_ols)
    if (d_method == "dmm" && min(lambda) < 1e-8 * max(lambda)) {
      warning("The design is near-exactly singular (smallest eigenvalue is ",
              "~0 relative to the largest); the 'dmm' formula divides by ",
              "eigenvalues directly and can produce an extreme, meaningless ",
              "d in this regime. Consider d = \"dopt\" (numerically ",
              "better-behaved here) or a manual d.", call. = FALSE)
    }
    if (d_method == "dmm") {
      num <- sum(1 / (lambda * (lambda + 1)))
      den <- sum(alpha^2 / (lambda + 1)^2)
      d_val <- 1 - sigma2_ols * (num / den)
    } else {
      signal <- sum(alpha^2 / (lambda + 1)^2)
      noise  <- sum(sigma2_ols / (lambda + 1)^2)
      penalty <- sum(sigma2_ols / (lambda * (lambda + 1)^2))
      d_val <- (signal - noise) / (signal + penalty)
    }
    if (clip) {
      d_val <- min(max(d_val, 0), 1)
    } else if (d_val < 0 || d_val > 1) {
      warning(sprintf(
        "The automatically-computed d (%s = %.4f) falls outside (0, 1); ",
        d_method, d_val),
        "using it as-is since clip = FALSE. Set clip = TRUE to snap it to ",
        "[0, 1], or supply a fixed numeric d.", call. = FALSE)
    }
  } else {
    d_val <- as.numeric(d)
    if (d_val < 0 || d_val > 1) {
      warning(sprintf("d = %.4f falls outside the theoretical (0, 1) range.", d_val),
              call. = FALSE)
    }
  }

  I_p <- diag(p)
  beta <- as.numeric(solve(XtX + I_p, (XtX + d_val * I_p) %*% beta_ols))
  names(beta) <- colnames(X)

  fitted_c <- as.numeric(X %*% beta)
  residuals <- yc - fitted_c
  df_residual <- n - p - 1
  sigma2 <- sum(residuals^2) / df_residual
  # See ?ridgeReg for the rationale: R^2 as squared correlation of actual
  # and fitted values, the standard convention for a shrunk estimator.
  r_squared <- stats::cor(yc, fitted_c)^2

  A_inv <- solve(XtX + I_p)
  G <- A_inv %*% (XtX + d_val * I_p)
  Vhat_slopes <- sigma2 * (G %*% solve(XtX) %*% t(G))
  se_slopes <- sqrt(diag(Vhat_slopes))
  intercept <- dsg$y_mean - sum(dsg$means * beta)
  var_intercept <- sigma2 / n + as.numeric(t(dsg$means) %*% Vhat_slopes %*% dsg$means)
  se_intercept <- sqrt(var_intercept)

  coefficients <- c(`(Intercept)` = intercept, beta)
  se <- c(`(Intercept)` = se_intercept, se_slopes)
  # No t/p-value columns -- see ?ridgeReg for the rationale; the same
  # applies to Liu regression's biased, non-standard sampling distribution.
  coef_table <- cbind(Estimate = coefficients, `Approx. SE` = se)

  structure(
    list(coefficients = coefficients, coef_table = coef_table,
         fitted.values = fitted_c + dsg$y_mean, residuals = residuals,
         d = d_val, d_method = if (is.null(d_method)) "manual" else d_method,
         sigma = sqrt(sigma2), df.residual = df_residual,
         r.squared = r_squared, call = cl, formula = dsg$formula,
         xlevels = dsg$xlevels, n = n, p = p,
         x = cbind(`(Intercept)` = rep(1, n), dsg$X_raw), y = yc + dsg$y_mean),
    class = "liuReg")
}

#' @export
print.liuReg <- function(x, ...) {
  cat("Ordinary Liu Regression (d = ", signif(x$d, 5), ", ", x$d_method, ")\n\n", sep = "")
  cat("Call:\n"); print(x$call)
  cat("\nCoefficients:\n")
  print(round(x$coefficients, 5))
  invisible(x)
}

#' @export
summary.liuReg <- function(object, ...) structure(object, class = c("summary.liuReg", class(object)))

#' @export
print.summary.liuReg <- function(x, ...) {
  cat("Ordinary Liu Regression\n\n")
  cat("Call:\n"); print(x$call); cat("\n")
  cat(sprintf("Biasing parameter d: %.5f  (%s)\n\n", x$d, x$d_method))
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
coef.liuReg <- function(object, ...) object$coefficients

#' @export
fitted.liuReg <- function(object, ...) object$fitted.values

#' @export
residuals.liuReg <- function(object, ...) object$residuals

#' @export
predict.liuReg <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(object$fitted.values)
  X <- .new_design_matrix(object$formula, newdata, object$xlevels)
  b <- object$coefficients
  as.numeric(b["(Intercept)"] + X[, names(b)[-1], drop = FALSE] %*% b[-1])
}

#' Diagnostic plots for a liuReg object
#' @param x a \code{liuReg} object.
#' @param which subset of 1:2.
#' @param ... unused.
#' @export
plot.liuReg <- function(x, which = 1:2, ...) {
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
