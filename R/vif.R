#' Variance Inflation Factor (Ordinary and Robust)
#'
#' Computes the Variance Inflation Factor (VIF) as the diagonal elements of
#' the inverse of the correlation matrix of the predictors (Marquaridt,
#' 1970). When \code{type = "O"} (the default) the ordinary Pearson
#' correlation matrix is used. When \code{type = "R"} the Robust VIF (RVIF)
#' of Jacob and Varadharajan (2024, Sankhya B, 86(2), 845-871) is returned
#' instead: observations are weighted by Tukey's biweight function applied to
#' their Stahel-Donoho projection outlyingness, and a weighted correlation
#' matrix is formed from those weights.
#'
#' Both types accept the same three forms of \code{object}:
#' \itemize{
#'   \item a \strong{two-sided formula} plus \code{data}: e.g.
#'     \code{vif(Sepal.Length ~ ., data = iris)}.
#'   \item a \strong{fitted model object} from this package (any of
#'     \code{raiseReg}, \code{ridgeReg}, \code{liuReg}, \code{robRaise},
#'     \code{robRidge}, \code{robLiu}) or a plain \code{lm} object.
#'   \item a \strong{numeric matrix} of predictors only (no response column).
#' }
#'
#' @param object a two-sided formula, a fitted model object, or a numeric
#'   predictor matrix. See Details.
#' @param data a data frame; required when \code{object} is a formula.
#' @param type \code{"O"} (ordinary, default) or \code{"R"} (robust).
#' @param cc Tukey biweight tuning constant (default 4.685). Only used when
#'   \code{type = "R"}.
#' @param seed integer for reproducible direction sampling. Only used when
#'   \code{type = "R"}.
#' @param threshold VIF threshold flagged in \code{print} and \code{plot}
#'   (default 10).
#' @return when \code{type = "O"}: a named numeric vector of class
#'   \code{"vif"} with \code{threshold} and \code{cn} (Condition Number)
#'   attributes. When \code{type = "R"}: a list of class \code{"rvif"} with
#'   elements \code{$rvif} (the RVIF values), \code{$cn} (robust Condition
#'   Number), \code{$weights}, and \code{$outlyingness}.
#' @references
#' Marquaridt, D.W. (1970). Generalized inverses, ridge regression, biased
#' linear estimation, and nonlinear estimation. \emph{Technometrics}, 12(3),
#' 591-612.
#'
#' Jacob, J. and Varadharajan, R. (2024). Robust Variance Inflation Factor:
#' A Promising Approach for Collinearity Diagnostics in the Presence of
#' Outliers. \emph{Sankhya B}, 86(2), 845-871.
#' @examples
#' ## Formula interface (works with dot formulas and factors)
#' vif(Sepal.Length ~ ., data = iris)
#'
#' ## Fitted model interface
#' fit <- raiseReg(Sepal.Length ~ ., data = iris)
#' vif(fit)
#'
#' ## lm interface
#' vif(lm(Sepal.Length ~ ., data = iris))
#'
#' \donttest{
#' ## type = "R" (robust VIF) requires the mrfDepth package
#' vif(Sepal.Length ~ ., data = iris, type = "R", seed = 1)
#' vif(fit, type = "R", seed = 1)
#' }
#' @export
vif <- function(object, data = NULL, type = c("O", "R"), cc = 4.685,
                 seed = NULL, threshold = 10) {
  type <- match.arg(type)

  # ---- extract predictor matrix from any input form ----------------------
  mm <- .vif_model_matrix(object, data)

  # ---- ordinary VIF ------------------------------------------------------
  if (type == "O") {
    R <- stats::cor(mm)
    if (any(!is.finite(R))) {
      Rinv <- diag(NA_real_, ncol(mm))
    } else {
      Rinv <- try(solve(R), silent = TRUE)
      if (inherits(Rinv, "try-error")) Rinv <- diag(NA_real_, ncol(mm))
    }
    v <- diag(Rinv)
    names(v) <- colnames(mm)
    return(structure(v, class = "vif", threshold = threshold,
                     cn = .condition_number(R)))
  }

  # ---- robust VIF (RVIF) -------------------------------------------------
  out <- .projection_outlyingness(mm, seed = seed)
  w <- .outlyingness_weight(out, cc = cc)
  if (sum(w) < ncol(mm)) {
    warning("Fewer than p observations received nonzero weight; RVIF may ",
            "be unreliable.", call. = FALSE)
  }
  Rw <- .wcor(mm, w)
  Rinv <- try(solve(Rw), silent = TRUE)
  if (inherits(Rinv, "try-error")) {
    warning("The weighted correlation matrix is computationally singular; ",
            "RVIF values set to NA.", call. = FALSE)
    v <- stats::setNames(rep(NA_real_, ncol(mm)), colnames(mm))
  } else {
    v <- diag(Rinv)
    names(v) <- colnames(mm)
  }
  structure(list(rvif = v, cn = .condition_number(Rw), weights = w,
                 outlyingness = out, call = match.call()),
            class = "rvif", threshold = threshold)
}

# Internal helper: extract a predictor-only design matrix from any input.
.vif_model_matrix <- function(object, data) {
  if (is.matrix(object) || is.data.frame(object)) {
    mm <- as.matrix(object)
    storage.mode(mm) <- "double"
    return(mm)
  }
  if (inherits(object, "formula")) {
    if (is.null(data))
      stop("`data` is required when `object` is a formula.", call. = FALSE)
    mf  <- stats::model.frame(object, data = data)
    mm  <- stats::model.matrix(object, data = data)
    if ("(Intercept)" %in% colnames(mm))
      mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
    if (ncol(mm) < 2)
      stop("VIF requires at least two predictors.", call. = FALSE)
    return(mm)
  }
  # fitted model objects from this package or plain lm
  pkg_classes <- c("raiseReg","ridgeReg","robRaise","robRidge","liuReg","robLiu")
  if (inherits(object, c(pkg_classes, "lm"))) {
    # For package objects $x is a full design matrix (with intercept column);
    # for lm we call model.matrix().
    if (inherits(object, pkg_classes)) {
      x <- object$x
    } else {
      x <- stats::model.matrix(object)
    }
    if ("(Intercept)" %in% colnames(x))
      x <- x[, colnames(x) != "(Intercept)", drop = FALSE]
    # Drop NA-coefficient columns (predictors that were aliased/dropped)
    if (!is.null(object$coefficients)) {
      b <- object$coefficients
      b <- b[names(b) != "(Intercept)"]
      keep <- names(b)[!is.na(b)]
      if (length(keep) > 0 && all(keep %in% colnames(x)))
        x <- x[, keep, drop = FALSE]
    }
    if (ncol(x) < 2)
      stop("VIF requires at least two predictors.", call. = FALSE)
    return(x)
  }
  stop("`object` must be a formula, a fitted model object from raiseReg, ",
       "ridgeReg, liuReg, robRaise, robRidge, robLiu, or lm, or a numeric ",
       "predictor matrix.", call. = FALSE)
}

# rvif() kept as a convenience alias pointing to vif(type = "R")
#' @rdname vif
#' @export
rvif <- function(object, data = NULL, cc = 4.685, seed = NULL,
                  threshold = 10) {
  vif(object = object, data = data, type = "R", cc = cc,
      seed = seed, threshold = threshold)
}

#' @export
print.vif <- function(x, ...) {
  thr <- attr(x, "threshold")
  cat("Ordinary Variance Inflation Factor (correlation-matrix based)\n\n")
  print.default(round(unclass(x), 4))
  flagged <- names(x)[x >= thr]
  if (length(flagged) > 0) {
    cat(sprintf("\nVariables at or above the threshold (%.0f): %s\n",
                thr, paste(flagged, collapse = ", ")))
  } else {
    cat(sprintf("\nNo variable reaches the threshold of %.0f.\n", thr))
  }
  cat(sprintf("Condition Number: %.4f\n", attr(x, "cn")))
  invisible(x)
}

#' @export
plot.vif <- function(x, ...) {
  thr <- attr(x, "threshold")
  ord <- order(unclass(x))
  vals <- unclass(x)[ord]
  cols <- ifelse(vals >= thr, "firebrick3", "steelblue3")
  op <- graphics::par(mar = c(5, max(6, max(nchar(names(vals))) * 0.6), 3, 2))
  on.exit(graphics::par(op))
  graphics::barplot(vals, horiz = TRUE, las = 1, col = cols, border = NA,
                    xlab = "VIF", main = "Ordinary Variance Inflation Factor")
  graphics::abline(v = thr, lty = 2, col = "grey30")
  invisible(x)
}

#' @export
print.rvif <- function(x, ...) {
  thr <- attr(x, "threshold")
  cat("Robust Variance Inflation Factor (Jacob & Varadharajan, 2024)\n")
  cat(sprintf("%d of %d observations downweighted to 0\n\n",
              sum(x$weights == 0), length(x$weights)))
  print.default(round(x$rvif, 4))
  flagged <- names(x$rvif)[x$rvif >= thr]
  if (length(flagged) > 0) {
    cat(sprintf("\nVariables at or above the threshold (%.0f): %s\n",
                thr, paste(flagged, collapse = ", ")))
  } else {
    cat(sprintf("\nNo variable reaches the threshold of %.0f.\n", thr))
  }
  cat(sprintf("Condition Number (weighted): %.4f\n", x$cn))
  invisible(x)
}

#' @export
plot.rvif <- function(x, ...) {
  thr <- attr(x, "threshold")
  ord <- order(x$rvif)
  vals <- x$rvif[ord]
  cols <- ifelse(vals >= thr, "firebrick3", "steelblue3")
  op <- graphics::par(mar = c(5, max(6, max(nchar(names(vals))) * 0.6), 3, 2))
  on.exit(graphics::par(op))
  graphics::barplot(vals, horiz = TRUE, las = 1, col = cols, border = NA,
                    xlab = "RVIF", main = "Robust Variance Inflation Factor")
  graphics::abline(v = thr, lty = 2, col = "grey30")
  invisible(x)
}
