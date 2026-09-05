#' Raise Regression
#'
#' Fits Raise Regression, an inference-preserving alternative to Ridge
#' Regression for combating multicollinearity. Two strategies are available:
#' \describe{
#'   \item{\code{"simultaneous"}}{Simultaneous Raise Regression (SRR) of
#'     Jacob and Varadharajan (2022) <doi:10.1007/s11135-022-01557-9>. A QR
#'     decomposition of the (mean-centred) design matrix is used to compute
#'     the Sequential Variance Inflation Factor (SVIF) for every predictor,
#'     and every predictor at or above \code{threshold} is raised in a
#'     single step.}
#'   \item{\code{"sequential"}}{The original single-variable raise strategy
#'     of Jacob and Varadharajan (2023) <doi:10.13189/ms.2023.110106>: the
#'     predictor with the largest ordinary VIF is raised, VIFs are
#'     recomputed, and the process repeats until every VIF is below
#'     \code{threshold}.}
#' }
#' In both cases, by construction, regressing \eqn{y} on the *raised*
#' design \eqn{\tilde X = Q\tilde R} exactly reproduces the ordinary least
#' squares fitted values \eqn{QQ'y}: raising only reallocates the fitted
#' signal among collinear predictors, it never changes the fitted subspace.
#' The residual standard error, \eqn{R^2}, adjusted \eqn{R^2}, F-statistic
#' and the standard errors underlying the coefficient table are therefore
#' all computed from this raised-design residual and are numerically
#' identical to their ordinary least squares counterparts -- only the
#' individual coefficient estimates change, becoming far more stable, so
#' the usual \eqn{t}- and F-testing machinery remains exactly valid.
#' \code{fitted()}, \code{residuals()} and \code{predict()}, in contrast,
#' apply the raise coefficients to the *original, unraised* predictors, as
#' is required to score new data; these values are close to, but not
#' numerically identical to, the training fit implied by the reported
#' \eqn{R^2} and residual standard error.
#'
#' If a predictor's automatically-computed raise parameter would exceed
#' \code{lambda_max}, it is dropped from the model entirely (with a message)
#' rather than raised to that extreme a value, and everything is
#' recomputed on the remaining predictors; a raise parameter needing to be
#' that large signals a predictor carrying essentially no information the
#' other predictors don't already supply. Dropped predictors get an
#' \code{NA} coefficient, exactly as \code{lm()} does for aliased terms; see
#' \code{$dropped}.
#'
#' @param formula a two-sided formula, as in \code{lm}.
#' @param data a data frame.
#' @param method \code{"sequential"} (default) or \code{"simultaneous"}.
#' @param threshold the VIF/SVIF threshold above which a predictor is
#'   raised (default 10, the conventional cutoff).
#' @param margin safety margin subtracted from \code{threshold} when solving
#'   for the raise parameter for the sequential method, so that the
#'   *post-raise* VIF is strictly, not just marginally, below the threshold
#'   (default 0.04). The simultaneous method solves its raise parameters in
#'   closed form directly against \code{threshold}, so \code{margin} does
#'   not apply to it.
#' @param lambda_max raise parameters above this are treated as "drop the
#'   predictor" rather than "raise it this much" (default 1000, matching
#'   the original implementation's grid bound).
#' @param lambda_step grid resolution in the raise-parameter search, only
#'   used when \code{raise_method = "grid"} (default 0.1); ignored by the
#'   default closed-form solver.
#' @param lambda optional named numeric vector of user-supplied raise
#'   parameters for some or all predictors (the actual multiplier applied,
#'   i.e. what is reported in \code{$k}). Named entries are used directly
#'   instead of being computed automatically, and are never dropped
#'   regardless of \code{lambda_max}.
#' @param raise_method \code{"closed_form"} (default, exact) or
#'   \code{"grid"} (a literal replica of the original grid search, for
#'   reproducibility with earlier published results); only affects
#'   \code{method = "simultaneous"}.
#' @param ... additional arguments, ignored. Present so that generic
#'   callers such as \code{stats::update()} (used internally by
#'   \code{car::ncvTest()}) can pass extra arguments without error.
#' @return an object of class \code{"raiseReg"}. In addition to the usual
#'   \code{coefficients}/\code{fitted.values}/\code{residuals}, it carries
#'   \code{$vif} and \code{$svif} (before/after matrices; SVIF is what
#'   actually drives the raising decision, ordinary VIF depends on *all*
#'   other predictors and can move even for a predictor that was not itself
#'   raised), \code{$cn} (before/after Condition Number), \code{$k} (the
#'   raise parameters actually applied), and \code{$dropped} (predictors
#'   excluded for exceeding \code{lambda_max}).
#' @examples
#' fit <- raiseReg(mpg ~ disp + hp + wt + drat, data = mtcars)
#' summary(fit)
#' fit$vif; fit$svif; fit$cn
#' fit_seq <- raiseReg(mpg ~ disp + hp + wt + drat, data = mtcars, method = "sequential")
#' summary(fit_seq)
#' @export
raiseReg <- function(formula, data, method = c("sequential", "simultaneous"),
                      threshold = 10, margin = 0.04,
                      lambda_max = 1000, lambda_step = 0.1, lambda = NULL,
                      raise_method = c("closed_form", "grid"), ...) {
  method <- match.arg(method)
  raise_method <- match.arg(raise_method)
  cl <- match.call()
  formula <- stats::as.formula(formula)
  dsg <- .raise_design(formula, data)
  X <- dsg$X
  yc <- dsg$y
  vif_before <- .vif_from_matrix(X)
  cn_before <- .condition_number(stats::cor(X))

  if (method == "simultaneous") {
    sv <- .svif_raise(X, threshold = threshold, lambda_max = lambda_max,
                       lambda_step = lambda_step, lambda = lambda, raise_method = raise_method)
    Xraised <- sv$Q %*% sv$Rtilde
    colnames(Xraised) <- setdiff(colnames(X), sv$dropped)
    XtX_tilde <- crossprod(sv$Rtilde)
    k <- sv$k[sv$k > 1 + 1e-8]
    svif_before <- sv$svif_before
    dropped <- sv$dropped
    beta_survivors <- as.numeric(solve(sv$Rtilde, crossprod(sv$Q, yc)))
    names(beta_survivors) <- colnames(Xraised)
  } else {
    sq <- .sequential_raise(X, yc, threshold = threshold, margin = margin,
                             lambda_max = lambda_max, lambda = lambda)
    Xraised <- sq$X
    XtX_tilde <- crossprod(Xraised)
    k <- sq$k
    svif_before <- .svif_raise(X, threshold = threshold)$svif_before
    dropped <- sq$dropped
    beta_survivors <- as.numeric(solve(XtX_tilde, crossprod(Xraised, yc)))
    names(beta_survivors) <- colnames(Xraised)
  }
  survivors <- colnames(Xraised)
  vif_after_surv <- .vif_from_matrix(Xraised)
  svif_after_surv <- .svif_raise(Xraised, threshold = threshold)$svif_before
  cn_after <- .condition_number(stats::cor(Xraised))

  vif_after <- stats::setNames(rep(NA_real_, ncol(X)), colnames(X))
  vif_after[survivors] <- vif_after_surv
  svif_after <- stats::setNames(rep(NA_real_, ncol(X)), colnames(X))
  svif_after[survivors] <- svif_after_surv

  # --- Inference basis -----------------------------------------------------
  # By construction, X_raised %*% beta reproduces the OLS projection of y
  # exactly (raising reallocates fit among collinear predictors without
  # ever changing the fitted subspace); this residual is what a valid
  # sigma^2 for Var(beta) must be based on, and it is numerically identical
  # to plugging beta into an ordinary OLS fit on the same centred design
  # restricted to the surviving predictors.
  fitted_inference <- as.numeric(Xraised %*% beta_survivors)
  resid_inference <- yc - fitted_inference
  n <- dsg$n
  p_used <- length(survivors)
  df_residual <- n - p_used - 1
  SSE <- sum(resid_inference^2)
  TSS <- sum(yc^2)
  sigma2 <- SSE / df_residual
  r_squared <- 1 - SSE / TSS
  adj_r_squared <- 1 - (1 - r_squared) * (n - 1) / df_residual
  f_stat <- ((TSS - SSE) / p_used) / sigma2
  f_pvalue <- stats::pf(f_stat, p_used, df_residual, lower.tail = FALSE)

  XtX_tilde_inv <- solve(XtX_tilde)
  Vhat_slopes <- sigma2 * XtX_tilde_inv
  se_survivors <- sqrt(diag(Vhat_slopes))
  means_surv <- dsg$means[survivors]
  var_intercept <- sigma2 / n + as.numeric(t(means_surv) %*% Vhat_slopes %*% means_surv)
  se_intercept <- sqrt(var_intercept)
  intercept <- dsg$y_mean - sum(means_surv * beta_survivors)

  beta_full <- stats::setNames(rep(NA_real_, ncol(X)), colnames(X))
  beta_full[survivors] <- beta_survivors
  se_full <- stats::setNames(rep(NA_real_, ncol(X)), colnames(X))
  se_full[survivors] <- se_survivors

  coefficients <- c(`(Intercept)` = intercept, beta_full)
  se <- c(`(Intercept)` = se_intercept, se_full)
  t_value <- coefficients / se
  p_value <- 2 * stats::pt(-abs(t_value), df_residual)
  coef_table <- cbind(Estimate = coefficients, `Std. Error` = se,
                       `t value` = t_value, `Pr(>|t|)` = p_value)

  # --- Prediction basis ------------------------------------------------------
  # fitted()/residuals()/predict() apply the raise coefficients to the
  # *original, unraised* predictors, as required to score new data; these
  # are close to, but not numerically identical to, the inference-basis
  # values above.
  fitted_raw <- as.numeric(X[, survivors, drop = FALSE] %*% beta_survivors)
  fitted_values <- fitted_raw + dsg$y_mean
  resid_final <- yc - fitted_raw

  # --- Diagnostics support ---------------------------------------------------
  # Since regressing y on the raised design reproduces the OLS projection of
  # y exactly (see the class documentation), the raised-design fit's hat
  # matrix, residuals, and sigma^2 -- the *inference basis* above -- are
  # identical to those of an ordinary lm() fit on the surviving predictors.
  # cooks.distance()/dfbetas()/covRatio() are therefore computed from this
  # same basis, so they stay consistent with the reported SEs/t-values, and
  # are numerically identical to fitting lm(y ~ survivors) directly.
  # h_full = 1/n + h_centered holds exactly because every raised column
  # remains exactly mean-centered (raising is a linear combination of
  # mean-centered vectors, so it cannot reintroduce a nonzero mean).
  hat_centered <- rowSums((Xraised %*% XtX_tilde_inv) * Xraised)
  hat_full <- 1 / n + hat_centered
  x_raw_full <- cbind(`(Intercept)` = rep(1, n), dsg$X_raw[, survivors, drop = FALSE])
  y_orig <- yc + dsg$y_mean

  structure(
    list(coefficients = coefficients, coef_table = coef_table,
         fitted.values = fitted_values, residuals = resid_final,
         vif = cbind(before = vif_before, after = vif_after),
         svif = cbind(before = svif_before, after = svif_after),
         cn = c(before = cn_before, after = cn_after),
         k = k, dropped = dropped,
         sigma = sqrt(sigma2), df.residual = df_residual,
         r.squared = r_squared, adj.r.squared = adj_r_squared,
         fstatistic = c(value = f_stat, numdf = p_used, dendf = df_residual),
         f.pvalue = f_pvalue, threshold = threshold, method = method,
         call = cl, formula = dsg$formula, xlevels = dsg$xlevels, n = n, p = p_used,
         x = x_raw_full, y = y_orig, hat = hat_full,
         residuals_inference = resid_inference, fitted_inference = fitted_inference),
    class = "raiseReg")
}

#' @export
print.raiseReg <- function(x, ...) {
  cat("Raise Regression (", x$method, ")\n\n", sep = "")
  cat("Call:\n"); print(x$call)
  cat("\nCoefficients:\n")
  print(round(x$coefficients, 5))
  invisible(x)
}

#' @export
summary.raiseReg <- function(object, ...) {
  structure(object, class = c("summary.raiseReg", class(object)))
}

#' @export
print.summary.raiseReg <- function(x, ...) {
  cat("Raise Regression (", x$method, " raise, threshold = ", x$threshold, ")\n\n", sep = "")
  cat("Call:\n"); print(x$call); cat("\n")

  cat("Residuals:\n")
  print(summary(as.numeric(x$residuals)))

  cat("\nCoefficients:\n")
  stats::printCoefmat(x$coef_table, digits = 5, signif.stars = TRUE, na.print = "NA")

  if (length(x$dropped) > 0) {
    cat("\nDropped (raise parameter would have exceeded lambda_max):\n  ",
        paste(x$dropped, collapse = ", "), "\n", sep = "")
  }
  if (length(x$k) > 0) {
    cat("\nRaised variables and raise parameters:\n")
    print(round(x$k, 4))
  } else if (length(x$dropped) == 0) {
    cat("\nNo predictor's SVIF exceeded the threshold; no raising was necessary.\n")
  }

  cat(sprintf("\nResidual standard error: %.4f on %d degrees of freedom\n",
              x$sigma, x$df.residual))
  cat(sprintf("Multiple R-squared: %.4f,  Adjusted R-squared: %.4f\n",
              x$r.squared, x$adj.r.squared))
  cat(sprintf("F-statistic: %.3f on %d and %d DF,  p-value: %.4g\n",
              x$fstatistic["value"], x$fstatistic["numdf"], x$fstatistic["dendf"], x$f.pvalue))
  cat(sprintf("Condition Number: %.4f -> %.4f\n", x$cn["before"], x$cn["after"]))
  invisible(x)
}

#' @export
coef.raiseReg <- function(object, ...) object$coefficients

#' @export
fitted.raiseReg <- function(object, ...) object$fitted.values

#' @export
residuals.raiseReg <- function(object, type, ...) {
  # Plain residuals(fit) calls keep the documented, existing behaviour
  # (prediction-basis residuals, from applying the raise coefficients to
  # the original predictors -- see ?raiseReg). A `type` argument is only
  # ever supplied by generic-dispatch callers like car::ncvTest() and
  # lmtest::bptest(), which need residuals consistent with the reported
  # sigma/hatvalues/df.residual (the inference basis) rather than the
  # prediction basis, so it is used here purely as that signal.
  if (missing(type)) object$residuals else object$residuals_inference
}

#' @export
predict.raiseReg <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(object$fitted.values)
  X <- .new_design_matrix(object$formula, newdata, object$xlevels)
  b <- object$coefficients
  keep <- names(b)[-1][!is.na(b[-1])]
  as.numeric(b["(Intercept)"] + X[, keep, drop = FALSE] %*% b[keep])
}

#' Diagnostic plots for a raiseReg object
#'
#' Produces the same family of goodness-of-fit / diagnostic plots as
#' \code{plot.lm}: residuals vs fitted, a normal Q-Q plot of the residuals,
#' a scale-location plot, and a before/after bar chart of VIF values.
#' @param x a \code{raiseReg} object.
#' @param which subset of 1:4 selecting which plots to draw.
#' @param ... passed on to the underlying plotting functions.
#' @return Invisibly returns the fitted \code{raiseReg} object `x`.
#'   Called for its side effect of drawing diagnostic plots.
#' @export
plot.raiseReg <- function(x, which = 1:4, ...) {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  n_plots <- length(which)
  if (n_plots > 1) graphics::par(mfrow = c(ceiling(n_plots / 2), min(2, n_plots)))

  fv <- x$fitted.values
  rs <- x$residuals
  std_rs <- rs / x$sigma

  if (1 %in% which) {
    graphics::plot(fv, rs, xlab = "Fitted values", ylab = "Residuals",
                   main = "Residuals vs Fitted", pch = 20, col = "steelblue4")
    graphics::abline(h = 0, lty = 2, col = "grey40")
    graphics::lines(stats::lowess(fv, rs), col = "firebrick3", lwd = 2)
  }
  if (2 %in% which) {
    stats::qqnorm(std_rs, main = "Normal Q-Q", pch = 20, col = "steelblue4")
    stats::qqline(std_rs, col = "firebrick3", lwd = 2)
  }
  if (3 %in% which) {
    graphics::plot(fv, sqrt(abs(std_rs)), xlab = "Fitted values",
                   ylab = expression(sqrt("|Standardized residuals|")),
                   main = "Scale-Location", pch = 20, col = "steelblue4")
    graphics::lines(stats::lowess(fv, sqrt(abs(std_rs))), col = "firebrick3", lwd = 2)
  }
  if (4 %in% which) {
    vb <- x$vif[, "before"]; va <- x$vif[, "after"]
    ord <- order(vb)
    mat <- rbind(before = vb[ord], after = va[ord])
    graphics::barplot(mat, beside = TRUE, horiz = TRUE, las = 1,
                       col = c("grey70", "steelblue3"), border = NA,
                       xlab = "VIF", main = "VIF before / after raising",
                       legend.text = TRUE, args.legend = list(x = "bottomright", bty = "n"))
    graphics::abline(v = x$threshold, lty = 2, col = "firebrick3")
  }
  invisible(x)
}
