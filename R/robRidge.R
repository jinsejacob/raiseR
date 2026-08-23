#' Robust Ridge Regression
#'
#' Fits a robust version of Ridge Regression using one of two strategies.
#'
#' \describe{
#'   \item{\code{type = "MM"} (default)}{A Tukey biweight M-estimate on the
#'     mean-centred design provides a seed coefficient vector and residual
#'     scale from which the biasing parameter
#'     \eqn{k = p\hat\sigma_R^2 / \hat\beta_R'\hat\beta_R} is computed (unless
#'     \code{k} is supplied directly), and the final estimate is obtained by
#'     iteratively reweighted ridge regression with observation weights from
#'     Tukey's biweight function applied to the MAD-standardized residuals
#'     at each step (Yohai, 1987; the same MM philosophy as \code{robRaise()}).}
#'   \item{\code{type = "SDO"}}{Observations are weighted once, up front, by
#'     Tukey's biweight function applied to their Stahel-Donoho projection
#'     outlyingness in the joint (response, predictors) space -- exactly as
#'     in \code{robRaise()} -- and ordinary ridge regression is then applied
#'     to the resulting weighted design and response.}
#' }
#'
#' @param formula a two-sided formula.
#' @param data a data frame.
#' @param k optional fixed biasing parameter; if \code{NULL} (default) it
#'   is estimated from the data (from the seed M-estimate for
#'   \code{type = "MM"}, or from the weighted OLS fit for
#'   \code{type = "SDO"}), following the same formula as \code{ridgeReg()}.
#' @param type \code{"MM"} (default) or \code{"SDO"}.
#' @param cc tuning constant for Tukey's biweight (default 4.685). Used by
#'   both types (IRLS reweighting for \code{"MM"}; outlyingness weighting
#'   for \code{"SDO"}).
#' @param maxit maximum number of IRLS iterations for \code{type = "MM"}
#'   (default 200).
#' @param tol convergence tolerance for \code{type = "MM"} (default 1e-4).
#' @param seed optional integer for reproducible direction sampling; only
#'   used when \code{type = "SDO"}.
#' @return an object of class \code{"robRidge"}.
#' @examples
#' fit_mm <- robRidge(mpg ~ disp + hp + wt + drat, data = mtcars)
#' summary(fit_mm)
#' \dontrun{
#' fit_sdo <- robRidge(mpg ~ disp + hp + wt + drat, data = mtcars,
#'                      type = "SDO", seed = 1)
#' summary(fit_sdo)
#' }
#' @export
robRidge <- function(formula, data, k = NULL, type = c("MM", "SDO"),
                      cc = 4.685, maxit = 200, tol = 1e-4,
                      seed = NULL) {
  type <- match.arg(type)
  cl <- match.call()
  formula <- stats::as.formula(formula)

  if (type == "SDO") {
    raw <- .raw_design(formula, data)
    X_raw <- raw$X; y <- raw$y; n <- raw$n; p <- raw$p

    joint <- cbind(X_raw, .response = y)
    out <- .projection_outlyingness(joint, seed = seed)
    w <- .outlyingness_weight(out, cc = cc)
    if (sum(w) < p + 1) {
      warning("Fewer than p+1 observations received nonzero weight; ",
              "the SDO robust ridge fit may be unreliable.", call. = FALSE)
    }
    wsum <- sum(w)
    means_w <- colSums(X_raw * w) / wsum
    y_mean_w <- sum(y * w) / wsum
    Xc <- sweep(X_raw, 2, means_w, "-")
    yc <- y - y_mean_w

    sw <- sqrt(w)
    Xw <- Xc * sw
    yw <- yc * sw
    XtX <- crossprod(Xw)
    # beta_ols here is only a plug-in seed for the default k formula (Hoerl
    # et al.), not the final estimate -- if downweighting has pushed the
    # *unpenalized* weighted design into exact singularity (a real risk once
    # multicollinear predictors lose several informative rows), fall back to
    # a small ridge-stabilized solve purely for this seed step.
    beta_ols <- tryCatch(
      as.numeric(solve(XtX, crossprod(Xw, yw))),
      error = function(e) {
        warning("The unpenalized weighted design is computationally singular ",
                "after SDO downweighting; stabilizing the seed OLS step with ",
                "a small ridge penalty to estimate the default k.", call. = FALSE)
        jitter <- 1e-6 * mean(diag(XtX))
        as.numeric(solve(XtX + diag(jitter, p), crossprod(Xw, yw)))
      })
    df_ols <- wsum - p - 1
    sigma2_ols <- sum((yw - Xw %*% beta_ols)^2) / df_ols
    if (is.null(k)) k <- p * sigma2_ols / sum(beta_ols^2)

    A <- XtX + diag(k, p)
    A_inv <- solve(A)
    beta <- as.numeric(A_inv %*% crossprod(Xw, yw))
    names(beta) <- colnames(X_raw)

    fitted_c <- as.numeric(Xc %*% beta)
    df_residual <- max(wsum - p - 1, 1)
    resid_w <- (yc - fitted_c) * sw
    sigma2 <- sum(resid_w^2) / df_residual
    Vhat_slopes <- sigma2 * (A_inv %*% XtX %*% A_inv)
    se_slopes <- sqrt(diag(Vhat_slopes))

    intercept <- y_mean_w - sum(means_w * beta)
    var_intercept <- sigma2 / wsum + as.numeric(t(means_w) %*% Vhat_slopes %*% means_w)
    se_intercept <- sqrt(var_intercept)

    fitted_values <- fitted_c + y_mean_w
    residuals_raw <- y - fitted_values
    r_squared_w <- stats::cor(y, fitted_values)^2
    iterations <- NA_integer_
    xlevels <- raw$xlevels
    formula_out <- raw$formula
  } else {
    # --- type == "MM" ------------------------------------------------------
    dsg <- .raise_design(formula, data)
    X <- dsg$X; yc <- dsg$y
    n <- dsg$n; p <- dsg$p

    seed_fit <- .rlm_bisquare(X, yc, cc = cc, maxit = maxit)
    if (is.null(k)) k <- p * seed_fit$scale^2 / sum(seed_fit$coefficients^2)

    b_old <- seed_fit$coefficients
    res <- seed_fit$residuals
    w <- seed_fit$weights
    iter <- 0L
    repeat {
      iter <- iter + 1L
      A <- crossprod(X, w * X) + diag(k, p)
      b_new <- as.numeric(solve(A, crossprod(X, w * yc)))
      res <- yc - as.numeric(X %*% b_new)
      s <- stats::mad(res)
      if (s > 1e-12) {
        u <- res / s
        w <- .tukey_weight(u, cc = cc)
      }
      denom <- ifelse(abs(b_old) > 1e-8, b_old, 1)
      delta <- max(abs((b_new - b_old) / denom))
      b_old <- b_new
      if (delta < tol || iter >= maxit) break
    }
    beta <- b_old
    names(beta) <- colnames(X)

    A <- crossprod(X, w * X) + diag(k, p)
    A_inv <- solve(A)
    neff <- sum(w)
    df_residual <- max(neff - p - 1, 1)
    sigma2 <- sum(w * res^2) / df_residual
    XtWX <- crossprod(X, w * X)
    Vhat_slopes <- sigma2 * (A_inv %*% XtWX %*% A_inv)
    se_slopes <- sqrt(diag(Vhat_slopes))

    intercept <- dsg$y_mean - sum(dsg$means * beta)
    var_intercept <- sigma2 / neff + as.numeric(t(dsg$means) %*% Vhat_slopes %*% dsg$means)
    se_intercept <- sqrt(var_intercept)

    fitted_values <- as.numeric(X %*% beta) + dsg$y_mean
    residuals_raw <- (yc + dsg$y_mean) - fitted_values
    r_squared_w <- stats::cor(yc + dsg$y_mean, fitted_values)^2
    iterations <- iter
    xlevels <- dsg$xlevels
    formula_out <- dsg$formula
  }

  coefficients <- c(`(Intercept)` = intercept, beta)
  se <- c(`(Intercept)` = se_intercept, se_slopes)
  coef_table <- cbind(Estimate = coefficients, `Approx. SE` = se)

  structure(
    list(coefficients = coefficients, coef_table = coef_table,
         fitted.values = fitted_values, residuals = residuals_raw,
         weights = w, k = k, type = type, sigma = sqrt(sigma2),
         df.residual = df_residual, r.squared = r_squared_w,
         iterations = iterations, call = cl, formula = formula_out,
         xlevels = xlevels, n = n, p = p,
         x = if (type == "SDO") cbind(`(Intercept)` = rep(1, n), X_raw) else cbind(`(Intercept)` = rep(1, n), dsg$X_raw),
         y = if (type == "SDO") y else (yc + dsg$y_mean)),
    class = "robRidge")
}

#' @export
print.robRidge <- function(x, ...) {
  cat("Robust Ridge Regression (", x$type, ", k = ", signif(x$k, 5), ")\n\n", sep = "")
  cat("Call:\n"); print(x$call)
  cat("\nCoefficients:\n")
  print(round(x$coefficients, 5))
  invisible(x)
}

#' @export
summary.robRidge <- function(object, ...) structure(object, class = c("summary.robRidge", class(object)))

#' @export
print.summary.robRidge <- function(x, ...) {
  if (x$type == "MM") {
    cat("Robust Ridge Regression (MM, ", x$iterations, " IRLS iterations)\n\n", sep = "")
  } else {
    cat("Robust Ridge Regression (SDO)\n\n")
  }
  cat("Call:\n"); print(x$call); cat("\n")
  cat(sprintf("Biasing parameter k: %.5f\n", x$k))
  cat(sprintf("%d of %d observations downweighted to exactly 0\n\n",
              sum(x$weights == 0), length(x$weights)))
  cat("Residuals (raw scale):\n"); print(summary(as.numeric(x$residuals)))
  cat("\nCoefficients (point estimates only -- see note below):\n")
  print(round(x$coef_table, 5))
  cat(sprintf("\nRobust residual scale: %.4f on %.1f effective degrees of freedom\n",
              x$sigma, x$df.residual))
  cat(sprintf("R-squared (squared correlation of actual and fitted values): %.4f\n",
              x$r.squared))
  invisible(x)
}

#' @export
coef.robRidge <- function(object, ...) object$coefficients

#' @export
fitted.robRidge <- function(object, ...) object$fitted.values

#' @export
residuals.robRidge <- function(object, ...) object$residuals

#' @export
predict.robRidge <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(object$fitted.values)
  X <- .new_design_matrix(object$formula, newdata, object$xlevels)
  b <- object$coefficients
  as.numeric(b["(Intercept)"] + X[, names(b)[-1], drop = FALSE] %*% b[-1])
}

#' Diagnostic plots for an robRidge object
#' @param x an \code{robRidge} object.
#' @param which subset of 1:2.
#' @param ... unused.
#' @export
plot.robRidge <- function(x, which = 1:2, ...) {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  if (length(which) > 1) graphics::par(mfrow = c(1, length(which)))
  fv <- x$fitted.values; rs <- x$residuals; std_rs <- rs / x$sigma
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
  invisible(x)
}
