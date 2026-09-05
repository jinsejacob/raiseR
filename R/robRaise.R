#' Robust Raise Regression (RRM)
#'
#' Fits Robust Raise Regression: observations are downweighted by Tukey's
#' biweight function applied to their Stahel-Donoho projection outlyingness
#' in the joint (response, predictors) space, the weighted design is
#' mean-centred using the corresponding weighted means, and every predictor
#' whose weighted Sequential Variance Inflation Factor (SVIF) is at or above
#' \code{threshold} is raised via a weighted QR decomposition -- combining
#' the Raise Regression's inference-preserving handling of collinearity with
#' resistance to outliers in \eqn{X}, \eqn{Y} and \eqn{XY} space.
#'
#' Exact finite-sample inference is used throughout, following the sandwich
#' covariance \eqn{\Sigma_{RRM} = \sigma^2 \tilde{R}_w^{-1}(Q_w'WQ_w)\tilde{R}_w^{-T}},
#' the effective residual degrees of freedom
#' \eqn{\nu_2 = tr[(I-P_w)W]}, the exactly-unbiased variance estimate
#' \eqn{\hat\sigma^2_w = SSE_w/\nu_2}, and a Satterthwaite-Welch degrees of
#' freedom correction \eqn{\nu^*} for the resulting Wald \eqn{t}-tests, as
#' derived from first principles for this estimator.
#'
#' If a predictor's automatically-computed raise parameter would exceed
#' \code{lambda_max}, it is dropped from the model entirely (with a message)
#' rather than raised to that extreme a value, and everything -- including
#' the weighted QR decomposition and every inference quantity above -- is
#' recomputed on the remaining predictors. This is not just a convenience:
#' down-weighting rows can and does induce exact (not just near) collinearity
#' in low-cardinality predictors once enough rows are zeroed out, which no
#' finite raise parameter can fix. Dropped predictors get an \code{NA}
#' coefficient, exactly as \code{lm()} does for aliased terms; see
#' \code{$dropped}.
#'
#' @inheritParams rvif
#' @param formula a two-sided formula.
#' @param data a data frame.
#' @param method only \code{"simultaneous"} is implemented for the robust
#'   raise regression (there is no robust analogue of the sequential
#'   one-variable-at-a-time strategy in this package); the argument exists
#'   so that passing anything else gives an informative error rather than
#'   R's generic "unused argument" message.
#' @param threshold the SVIF threshold above which a predictor is raised
#'   (default 10).
#' @param margin unused for the robust method (its raise parameters are
#'   solved in closed form directly against \code{threshold}); retained
#'   only for argument-signature symmetry with \code{raiseReg()}.
#' @param lambda_max raise parameters above this are treated as "drop the
#'   predictor" rather than "raise it this much" (default 1000).
#' @param lambda_step grid resolution in the raise-parameter search, only
#'   used when \code{raise_method = "grid"} (default 0.1); ignored by the
#'   default closed-form solver.
#' @param lambda optional named numeric vector of user-supplied raise
#'   parameters for some or all predictors (the actual multiplier applied,
#'   i.e. what is reported in \code{$k}). Named entries are used directly
#'   instead of being computed automatically, and are never dropped
#'   regardless of \code{lambda_max}.
#' @param raise_method \code{"closed_form"} (default, exact) or
#'   \code{"grid"} (a literal replica of the original grid search).
#' @return an object of class \code{"robRaise"}. In addition to the usual
#'   \code{coefficients}/\code{fitted.values}/\code{residuals}, it carries
#'   \code{$vif} and \code{$svif} (before/after matrices, on the weighted
#'   design), \code{$cn} (before/after Condition Number of the weighted
#'   design), \code{$k} (the raise parameters actually applied), and
#'   \code{$dropped} (predictors excluded for exceeding \code{lambda_max}).
#' @examples
#' \donttest{
#' fit <- robRaise(mpg ~ disp + hp + wt + drat, data = mtcars, seed = 1)
#' summary(fit)
#' fit$vif; fit$svif; fit$cn
#' }
#' @export
robRaise <- function(formula, data, method = "simultaneous",
                    cc = 4.685, seed = NULL, threshold = 10, margin = 0.04,
                    lambda_max = 1000, lambda_step = 0.1, lambda = NULL,
                    raise_method = c("closed_form", "grid")) {
  if (!identical(method, "simultaneous")) {
    stop("robRaise() only implements the simultaneous (QR/SVIF-based) raise ",
         "strategy; there is no method = \"", method, "\" option. Use ",
         "raiseReg(..., method = \"sequential\") for a non-robust sequential fit.",
         call. = FALSE)
  }
  raise_method <- match.arg(raise_method)
  cl <- match.call()
  formula <- stats::as.formula(formula)
  raw <- .raw_design(formula, data)
  X_raw <- raw$X
  y <- raw$y
  n <- raw$n
  p <- raw$p

  out <- .projection_outlyingness(cbind(y, X_raw), seed = seed)
  w <- .outlyingness_weight(out, cc = cc)
  if (sum(w > 0) <= p + 1) {
    stop("Fewer than p + 1 observations received nonzero weight; ",
         "the robust fit is not identified. Consider a larger cc.", call. = FALSE)
  }

  wmx <- colSums(X_raw * w) / sum(w)
  wmy <- sum(y * w) / sum(w)
  Xc <- sweep(X_raw, 2, wmx, "-")
  yc <- y - wmy

  sw <- sqrt(w)
  Xw <- Xc * sw           # row-scaling: recycling matches nrow(Xc) == length(sw)
  yw <- yc * sw

  vif_before <- .vif_from_matrix(Xw)
  cn_before <- .condition_number(stats::cor(Xw))

  sv <- .svif_raise(Xw, threshold = threshold, lambda_max = lambda_max,
                     lambda_step = lambda_step, lambda = lambda, raise_method = raise_method)
  Qw <- sv$Q
  Rtilde <- sv$Rtilde
  survivors <- setdiff(colnames(X_raw), sv$dropped)
  dropped <- sv$dropped
  beta_survivors <- as.numeric(solve(Rtilde, crossprod(Qw, yw)))
  names(beta_survivors) <- survivors
  k <- sv$k[sv$k > 1 + 1e-8]

  fitted_w <- as.numeric(Qw %*% crossprod(Qw, yw))
  resid_w <- yw - fitted_w
  SSE_w <- sum(resid_w^2)
  TSS_w <- sum(yw^2)
  r_squared_w <- 1 - SSE_w / TSS_w

  M <- crossprod(Qw, w * Qw)              # Qw' W Qw  (p_used x p_used)
  M2 <- crossprod(Qw, (w^2) * Qw)         # Qw' W^2 Qw (p_used x p_used)
  nu2 <- sum(w) - sum(diag(M))
  if (nu2 <= 0) stop("Effective residual degrees of freedom are non-positive.", call. = FALSE)
  sigma2_w <- SSE_w / nu2
  trace_NW_sq <- sum(w^2) - 2 * sum(diag(M2)) + sum(diag(M %*% M))
  nu_star <- if (trace_NW_sq > 0) nu2^2 / trace_NW_sq else nu2

  Rtilde_inv <- solve(Rtilde)
  Sigma_RRM <- sigma2_w * (Rtilde_inv %*% M %*% t(Rtilde_inv))
  se_survivors <- sqrt(diag(Sigma_RRM))

  wmx_surv <- wmx[survivors]
  intercept <- wmy - sum(wmx_surv * beta_survivors)
  var_intercept <- sigma2_w / sum(w) + as.numeric(t(wmx_surv) %*% Sigma_RRM %*% wmx_surv)
  se_intercept <- sqrt(var_intercept)

  beta_full <- stats::setNames(rep(NA_real_, p), colnames(X_raw))
  beta_full[survivors] <- beta_survivors
  se_full <- stats::setNames(rep(NA_real_, p), colnames(X_raw))
  se_full[survivors] <- se_survivors

  coefficients <- c(`(Intercept)` = intercept, beta_full)
  se <- c(`(Intercept)` = se_intercept, se_full)
  t_value <- coefficients / se
  p_value <- 2 * stats::pt(-abs(t_value), nu_star)
  coef_table <- cbind(Estimate = coefficients, `Std. Error` = se,
                       `t value` = t_value, `Pr(>|t|)` = p_value)

  fitted_values <- as.numeric(intercept + X_raw[, survivors, drop = FALSE] %*% beta_survivors)
  residuals <- y - fitted_values

  # Two distinct diagnostics, reported separately because they can tell
  # very different stories: ordinary (correlation-matrix) weighted VIF sees
  # a predictor's dependence on *all* other predictors, while weighted SVIF
  # -- what actually drives the raising decision -- only sees its
  # dependence on the predictors *before* it in the QR order. A predictor
  # can show a high weighted VIF and a low weighted SVIF (mostly explained
  # by columns ahead of it in the order) and therefore never get raised.
  Xraised <- Qw %*% Rtilde
  colnames(Xraised) <- survivors
  vif_after_surv <- .vif_from_matrix(Xraised)
  svif_before <- sv$svif_before
  svif_after_surv <- .svif_raise(Xraised, threshold = threshold)$svif_before
  cn_after <- .condition_number(stats::cor(Xraised))

  vif_after <- stats::setNames(rep(NA_real_, p), colnames(X_raw))
  vif_after[survivors] <- vif_after_surv
  svif_after <- stats::setNames(rep(NA_real_, p), colnames(X_raw))
  svif_after[survivors] <- svif_after_surv

  structure(
    list(coefficients = coefficients, coef_table = coef_table,
         fitted.values = fitted_values, residuals = residuals,
         weights = w, outlyingness = out, n_downweighted = sum(w == 0),
         vif = cbind(before = vif_before, after = vif_after),
         svif = cbind(before = svif_before, after = svif_after),
         cn = c(before = cn_before, after = cn_after),
         k = k, dropped = dropped,
         sigma = sqrt(sigma2_w), nu2 = nu2, nu_star = nu_star,
         df.residual = nu_star,
         r.squared = r_squared_w, threshold = threshold,
         call = cl, formula = raw$formula, xlevels = raw$xlevels, n = n, p = length(survivors),
         x = cbind(`(Intercept)` = rep(1, n), X_raw[, survivors, drop = FALSE]), y = y),
    class = "robRaise")
}

#' @export
print.robRaise <- function(x, ...) {
  cat("Robust Raise Regression (RRM)\n\n")
  cat("Call:\n"); print(x$call)
  cat("\nCoefficients:\n")
  print(round(x$coefficients, 5))
  invisible(x)
}

#' @export
summary.robRaise <- function(object, ...) {
  structure(object, class = c("summary.robRaise", class(object)))
}

#' @export
print.summary.robRaise <- function(x, ...) {
  cat("Robust Raise Regression (RRM), threshold = ", x$threshold,
      "\n\n", sep = "")
  cat("Call:\n"); print(x$call); cat("\n")

  cat(sprintf("%d of %d observations downweighted to exactly 0 (min weight %.3f, mean weight %.3f)\n\n",
              x$n_downweighted, x$n, min(x$weights), mean(x$weights)))

  cat("Residuals (raw scale):\n")
  print(summary(as.numeric(x$residuals)))

  cat("\nCoefficients (Wald tests, Satterthwaite-Welch df):\n")
  stats::printCoefmat(x$coef_table, digits = 5, signif.stars = TRUE, na.print = "NA")

  if (length(x$dropped) > 0) {
    cat("\nDropped (raise parameter would have exceeded lambda_max):\n  ",
        paste(x$dropped, collapse = ", "), "\n", sep = "")
  }
  if (length(x$k) > 0) {
    cat("\nRaised variables and raise parameters:\n")
    print(round(x$k, 4))
  } else if (length(x$dropped) == 0) {
    cat("\nNo predictor exceeded the weighted-SVIF threshold; no raising was necessary.\n")
  }

  cat(sprintf("\nRobust residual scale: %.4f\n", x$sigma))
  cat(sprintf("Effective residual degrees of freedom (nu2): %.2f\n", x$nu2))
  cat(sprintf("Satterthwaite-Welch degrees of freedom (nu*): %.2f\n", x$nu_star))
  cat(sprintf("Weighted R-squared: %.4f\n", x$r.squared))
  cat(sprintf("Condition Number (weighted design): %.4f -> %.4f\n", x$cn["before"], x$cn["after"]))
  invisible(x)
}

#' @export
coef.robRaise <- function(object, ...) object$coefficients

#' @export
fitted.robRaise <- function(object, ...) object$fitted.values

#' @export
residuals.robRaise <- function(object, ...) object$residuals

#' @export
weights.robRaise <- function(object, ...) object$weights

#' @export
df.residual.robRaise <- function(object, ...) object$df.residual

#' @exportS3Method car::ncvTest
ncvTest.robRaise <- function(model, var.formula, ...) {
  # Same algorithm as ncvTest.raiseReg (see ?raiseReg-diagnostics); robRaise
  # has a single fitted/residual basis (unlike raiseReg's prediction vs.
  # inference distinction), so plain residuals()/fitted() are used directly.
  sumry <- summary(model)
  resid_p <- stats::residuals(model, type = "pearson")
  S.sq <- stats::df.residual(model) * (sumry$sigma)^2 / sum(!is.na(resid_p))
  U <- (resid_p^2) / S.sq
  if (missing(var.formula)) {
    mod <- stats::lm(U ~ model$fitted.values)
    var.formula <- ~fitted.values
    df <- 1
  } else {
    form <- stats::as.formula(paste("U ~ ", as.character(var.formula)[[2]]))
    dat <- as.data.frame(model$x[, -1, drop = FALSE])
    dat$U <- U
    mod <- stats::lm(form, data = dat)
    df <- sum(!is.na(stats::coefficients(mod))) - 1
  }
  SS <- stats::anova(mod)$"Sum Sq"
  RegSS <- sum(SS) - SS[length(SS)]
  Chisq <- RegSS / 2
  result <- list(formula = var.formula, formula.name = "Variance",
                  ChiSquare = Chisq, Df = df,
                  p = stats::pchisq(Chisq, df, lower.tail = FALSE),
                  test = "Non-constant Variance Score Test")
  class(result) <- "chisqTest"
  result
}

#' @export
predict.robRaise <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(object$fitted.values)
  X <- .new_design_matrix(object$formula, newdata, object$xlevels)
  b <- object$coefficients
  keep <- names(b)[-1][!is.na(b[-1])]
  as.numeric(b["(Intercept)"] + X[, keep, drop = FALSE] %*% b[keep])
}

#' Diagnostic plots for an robRaise object
#'
#' Residuals-vs-fitted, a normal Q-Q plot of the residuals with points sized
#' by their downweighting, a weight-vs-outlyingness plot, and a before/after
#' weighted-VIF bar chart.
#' @param x an \code{robRaise} object.
#' @param which subset of 1:4 selecting which plots to draw.
#' @param ... passed on to the underlying plotting functions.
#' @return Invisibly returns the fitted \code{robRaise} object `x`.
#'   Called for its side effect of drawing diagnostic plots.
#' @export
plot.robRaise <- function(x, which = 1:4, ...) {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  n_plots <- length(which)
  if (n_plots > 1) graphics::par(mfrow = c(ceiling(n_plots / 2), min(2, n_plots)))

  fv <- x$fitted.values
  rs <- x$residuals
  std_rs <- rs / x$sigma
  ptcol <- grDevices::adjustcolor(ifelse(x$weights == 0, "firebrick3", "steelblue4"), alpha.f = 0.8)

  if (1 %in% which) {
    graphics::plot(fv, rs, xlab = "Fitted values", ylab = "Residuals",
                   main = "Residuals vs Fitted", pch = 20, col = ptcol)
    graphics::abline(h = 0, lty = 2, col = "grey40")
  }
  if (2 %in% which) {
    stats::qqnorm(std_rs, main = "Normal Q-Q", pch = 20, col = ptcol)
    stats::qqline(std_rs, col = "grey40", lwd = 2)
  }
  if (3 %in% which) {
    graphics::plot(x$outlyingness, x$weights, xlab = "Projection outlyingness",
                   ylab = "Tukey biweight", main = "Downweighting", pch = 20, col = ptcol)
  }
  if (4 %in% which) {
    vb <- x$vif[, "before"]; va <- x$vif[, "after"]
    ord <- order(vb)
    mat <- rbind(before = vb[ord], after = va[ord])
    graphics::barplot(mat, beside = TRUE, horiz = TRUE, las = 1,
                       col = c("grey70", "steelblue3"), border = NA,
                       xlab = "Weighted VIF", main = "Weighted VIF before / after raising",
                       legend.text = TRUE, args.legend = list(x = "bottomright", bty = "n"))
    graphics::abline(v = x$threshold, lty = 2, col = "firebrick3")
  }
  invisible(x)
}
