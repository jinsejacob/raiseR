#' Robust Liu Regression
#'
#' Fits a robust version of Liu Regression using one of two strategies,
#' mirroring \code{robRidge()}.
#'
#' \describe{
#'   \item{\code{type = "MM"} (default)}{Following Filzmoser and Kurnaz's MM-Liu
#'     estimator: \code{MASS::rlm(..., method = "MM")} supplies both the
#'     plug-in coefficient vector \eqn{\hat\beta_{MM}} and observation
#'     weights \eqn{W} (from the MM-estimator's own IRLS iterations), and
#'     \deqn{\hat\beta_{MM\text{-}Liu} = (X'WX+I)^{-1}(X'WX+dI)\hat\beta_{MM}.}
#'     The biasing parameter \code{d} minimizes
#'     \eqn{\text{MSE}(d) = \text{bias}(d)'\text{bias}(d) + tr(\text{cov}(d))}
#'     where \eqn{\text{bias}(d) = (d-1)(X'WX+I)^{-1}\hat\beta_{MM}} and
#'     \eqn{\text{cov}(d) = G(d)\,\hat V\,G(d)'} with \eqn{G(d) = (X'WX+I)^{-1}(X'WX+dI)}
#'     and \eqn{\hat V} the (asymptotic, sandwich-type) covariance of
#'     \eqn{\hat\beta_{MM}} reported by \code{MASS::rlm()}. Because
#'     \eqn{G(d)} is affine in \code{d}, this objective is an exact
#'     quadratic in \code{d} and is minimized here in closed form (not by
#'     numerical search).}
#'   \item{\code{type = "SDO"}}{Observations are weighted once, up front,
#'     by Tukey's biweight function applied to their Stahel-Donoho
#'     projection outlyingness in the joint (response, predictors) space --
#'     exactly as in \code{robRidge(type = "SDO")} -- and ordinary Liu
#'     regression (\code{liuReg()}'s \code{dopt}/\code{dmm} formulas) is
#'     then applied to the resulting weighted design and response.}
#' }
#'
#' @param formula a two-sided formula.
#' @param data a data frame.
#' @param type \code{"MM"} (default) or \code{"SDO"}.
#' @param d for \code{type = "SDO"}: either \code{"dopt"} (default, the
#'   more numerically robust of the two -- see \code{?liuReg}),
#'   \code{"dmm"}, or a fixed numeric value. For \code{type = "MM"}: either
#'   \code{NULL} (default, the closed-form MSE-minimizing value described
#'   above -- the MM-specific analogue of "dopt", so the \code{"dmm"}/
#'   \code{"dopt"} strings are not applicable and raise an error) or a
#'   fixed numeric value.
#' @param clip if \code{TRUE} (default), an automatically-computed \code{d}
#'   is snapped to \eqn{[0, 1]} if it falls outside that range.
#' @param cc tuning constant for Tukey's biweight (default 4.685). Used by
#'   both types.
#' @param maxit,acc passed to \code{MASS::rlm()}; only used when
#'   \code{type = "MM"}.
#' @param seed optional integer for reproducible direction sampling; only
#'   used when \code{type = "SDO"}.
#' @return an object of class \code{"robLiu"}.
#' @references
#' Filzmoser, P. and Kurnaz, F.S. A Robust Liu Regression Estimator.
#' Yohai, V.J. (1987). High breakdown-point and high efficiency robust
#' estimates for regression. \emph{Annals of Statistics}.
#' @examples
#' fit_mm <- robLiu(mpg ~ disp + hp + wt + drat, data = mtcars)
#' summary(fit_mm)
#' \donttest{
#' fit_sdo <- robLiu(mpg ~ disp + hp + wt + drat, data = mtcars,
#'                    type = "SDO", seed = 1)
#' summary(fit_sdo)
#' }
#' @export
robLiu <- function(formula, data, type = c("MM", "SDO"), d = NULL, clip = TRUE,
                    cc = 4.685, maxit = 200, acc = 1e-4,
                    seed = NULL) {
  type <- match.arg(type)
  cl <- match.call()
  formula <- stats::as.formula(formula)

  if (type == "MM") {
    dsg <- .raise_design(formula, data)
    X <- dsg$X; yc <- dsg$y
    n <- dsg$n; p <- dsg$p

    rfit <- MASS::rlm(X, yc, method = "MM", psi = MASS::psi.bisquare, c = cc,
                       maxit = maxit, acc = acc)
    beta_mm <- as.numeric(stats::coef(rfit))
    names(beta_mm) <- colnames(X)
    w <- rfit$w
    V <- stats::vcov(rfit)                      # asymptotic covariance of beta_mm

    A <- crossprod(X, w * X)
    I_p <- diag(p)
    G1 <- solve(A + I_p)                        # symmetric
    G0 <- G1 %*% A

    if (is.null(d)) {
      c2 <- as.numeric(t(beta_mm) %*% G1 %*% G1 %*% beta_mm)
      a2 <- sum(diag(G1 %*% V %*% G1))
      a1 <- 2 * sum(diag(G1 %*% V %*% t(G0)))
      d_val <- (2 * c2 - a1) / (2 * (c2 + a2))
      d_method <- "MM (closed-form MSE-minimizing)"
      if (clip) {
        d_val <- min(max(d_val, 0), 1)
      } else if (d_val < 0 || d_val > 1) {
        warning(sprintf("The MM biasing parameter (d = %.4f) falls outside (0, 1); ",
                         d_val), "using it as-is since clip = FALSE.", call. = FALSE)
      }
    } else if (is.character(d)) {
      stop("For type = \"MM\", d must be NULL (the default -- this uses the ",
           "closed-form MSE-minimizing value, which is the MM-specific ",
           "analogue of \"dopt\") or a fixed numeric value. The \"dmm\"/\"dopt\" ",
           "string options are only meaningful for type = \"SDO\", where they ",
           "select between the classical Liu (1993) formulas computed on the ",
           "weighted data.", call. = FALSE)
    } else {
      d_val <- as.numeric(d)
      d_method <- "manual"
    }

    G <- G1 %*% (A + d_val * I_p)
    beta <- as.numeric(G %*% beta_mm)
    names(beta) <- colnames(X)

    Vhat_slopes <- G %*% V %*% t(G)
    se_slopes <- sqrt(diag(Vhat_slopes))
    neff <- sum(w)
    df_residual <- max(neff - p - 1, 1)
    fitted_c <- as.numeric(X %*% beta)
    residuals_raw <- (yc + dsg$y_mean) - (fitted_c + dsg$y_mean)
    sigma2 <- sum(w * residuals_raw^2) / df_residual

    intercept <- dsg$y_mean - sum(dsg$means * beta)
    var_intercept <- sigma2 / neff + as.numeric(t(dsg$means) %*% Vhat_slopes %*% dsg$means)
    se_intercept <- sqrt(var_intercept)

    fitted_values <- fitted_c + dsg$y_mean
    r_squared_w <- stats::cor(yc + dsg$y_mean, fitted_values)^2
    xlevels <- dsg$xlevels
    formula_out <- dsg$formula
  } else {
    # --- type == "SDO": weight the data, then ordinary Liu regression ------
    raw <- .raw_design(formula, data)
    X_raw <- raw$X; y <- raw$y; n <- raw$n; p <- raw$p

    joint <- cbind(X_raw, .response = y)
    out <- .projection_outlyingness(joint, seed = seed)
    w <- .outlyingness_weight(out, cc = cc)
    if (sum(w) < p + 1) {
      warning("Fewer than p+1 observations received nonzero weight; ",
              "the SDO robust Liu fit may be unreliable.", call. = FALSE)
    }
    wsum <- sum(w)
    means_w <- colSums(X_raw * w) / wsum
    y_mean_w <- sum(y * w) / wsum
    Xc <- sweep(X_raw, 2, means_w, "-")
    yc_ <- y - y_mean_w
    sw <- sqrt(w)
    Xw <- Xc * sw
    yw <- yc_ * sw

    XtX <- crossprod(Xw)
    # Both the seed OLS coefficient and the later variance sandwich need
    # XtX^{-1}; if severe multicollinearity plus SDO downweighting has
    # pushed XtX to (near-)exact singularity, stabilize once here and reuse
    # the same stabilized inverse everywhere below, rather than patching
    # each solve() call separately with a different fallback.
    XtX_inv <- tryCatch(
      solve(XtX),
      error = function(e) {
        warning("The unpenalized weighted design is computationally singular ",
                "after SDO downweighting; stabilizing with a small ridge ",
                "penalty for the seed OLS step and variance calculations.", call. = FALSE)
        jitter <- 1e-6 * mean(diag(XtX))
        solve(XtX + diag(jitter, p))
      })
    beta_ols <- as.numeric(XtX_inv %*% crossprod(Xw, yw))
    df_ols <- wsum - p                          # matches liuReg()'s convention
    sigma2_ols <- sum((yw - Xw %*% beta_ols)^2) / df_ols

    d_method <- NULL
    if (is.null(d)) d <- "dopt"
    if (is.character(d)) {
      d_method <- match.arg(d, c("dmm", "dopt"))
      eig <- eigen(XtX, symmetric = TRUE)
      lambda <- eig$values
      alpha <- as.numeric(t(eig$vectors) %*% beta_ols)
      if (d_method == "dmm" && min(lambda) < 1e-8 * max(lambda)) {
        warning("The (possibly SDO-weighted) design is near-exactly singular ",
                "(smallest eigenvalue is ~0 relative to the largest); the ",
                "'dmm' formula divides by eigenvalues directly and can ",
                "produce an extreme, meaningless d in this regime. Consider ",
                "d = \"dopt\" (numerically better-behaved here) or a manual d.",
                call. = FALSE)
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
        warning(sprintf("The automatically-computed d (%s = %.4f) falls outside (0, 1); ",
                         d_method, d_val), "using it as-is since clip = FALSE.", call. = FALSE)
      }
    } else {
      d_val <- as.numeric(d)
      d_method <- "manual"
    }

    I_p <- diag(p)
    beta <- as.numeric(solve(XtX + I_p, (XtX + d_val * I_p) %*% beta_ols))
    names(beta) <- colnames(X_raw)

    fitted_c <- as.numeric(Xc %*% beta)
    df_residual <- max(wsum - p - 1, 1)
    resid_w <- (yc_ - fitted_c) * sw
    sigma2 <- sum(resid_w^2) / df_residual

    A_inv <- solve(XtX + I_p)
    G <- A_inv %*% (XtX + d_val * I_p)
    Vhat_slopes <- sigma2 * (G %*% XtX_inv %*% t(G))
    se_slopes <- sqrt(diag(Vhat_slopes))

    intercept <- y_mean_w - sum(means_w * beta)
    var_intercept <- sigma2 / wsum + as.numeric(t(means_w) %*% Vhat_slopes %*% means_w)
    se_intercept <- sqrt(var_intercept)

    fitted_values <- fitted_c + y_mean_w
    residuals_raw <- y - fitted_values
    r_squared_w <- stats::cor(y, fitted_values)^2
    d_val <- d_val
    xlevels <- raw$xlevels
    formula_out <- raw$formula
  }

  coefficients <- c(`(Intercept)` = intercept, beta)
  se <- c(`(Intercept)` = se_intercept, se_slopes)
  coef_table <- cbind(Estimate = coefficients, `Approx. SE` = se)

  structure(
    list(coefficients = coefficients, coef_table = coef_table,
         fitted.values = fitted_values, residuals = residuals_raw,
         weights = w, d = d_val, d_method = if (is.null(d_method)) "manual" else d_method,
         type = type, sigma = sqrt(sigma2), df.residual = df_residual,
         r.squared = r_squared_w, call = cl, formula = formula_out,
         xlevels = xlevels, n = n, p = p,
         x = cbind(`(Intercept)` = rep(1, n), if (type == "MM") dsg$X_raw else X_raw),
         y = if (type == "MM") (yc + dsg$y_mean) else y),
    class = "robLiu")
}

#' @export
print.robLiu <- function(x, ...) {
  cat("Robust Liu Regression (", x$type, ", d = ", signif(x$d, 5), ")\n\n", sep = "")
  cat("Call:\n"); print(x$call)
  cat("\nCoefficients:\n")
  print(round(x$coefficients, 5))
  invisible(x)
}

#' @export
summary.robLiu <- function(object, ...) structure(object, class = c("summary.robLiu", class(object)))

#' @export
print.summary.robLiu <- function(x, ...) {
  cat("Robust Liu Regression (", x$type, ")\n\n", sep = "")
  cat("Call:\n"); print(x$call); cat("\n")
  cat(sprintf("Biasing parameter d: %.5f  (%s)\n", x$d, x$d_method))
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
coef.robLiu <- function(object, ...) object$coefficients

#' @export
fitted.robLiu <- function(object, ...) object$fitted.values

#' @export
residuals.robLiu <- function(object, ...) object$residuals

#' @export
predict.robLiu <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) return(object$fitted.values)
  X <- .new_design_matrix(object$formula, newdata, object$xlevels)
  b <- object$coefficients
  as.numeric(b["(Intercept)"] + X[, names(b)[-1], drop = FALSE] %*% b[-1])
}

#' Diagnostic plots for a robLiu object
#' @param x a \code{robLiu} object.
#' @param which subset of 1:2.
#' @param ... unused.
#' @return Invisibly returns the fitted \code{robLiu} object `x`.
#'   Called for its side effect of drawing diagnostic plots.
#' @export
plot.robLiu <- function(x, which = 1:2, ...) {
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
