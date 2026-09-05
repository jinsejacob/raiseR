#' Condition Number and Condition Indices
#'
#' Computes the overall Condition Number (Belsley, 1991) together with the
#' full vector of Condition Indices, one per principal axis of the
#' (correlation-matrix) eigenstructure of the predictors -- either in the
#' ordinary sense (\code{type = "O"}) or a robust, outlier-resistant sense
#' (\code{type = "R"}).
#'
#' For \code{type = "R"}, weights are computed from Tukey's biweight
#' function applied to each observation's Stahel-Donoho projection
#' outlyingness -- the same scheme used throughout this package's robust
#' methods (\code{rvif()}, \code{robRaise()}, \code{robRidge()},
#' \code{robLiu()}) -- with one important difference: outlyingness here is
#' computed jointly over the predictor \emph{and} response space (i.e. from
#' \code{cbind(X, y)}, not \code{X} alone), following Jacob and Varadharajan
#' (2024). This is because a point that is only an outlier in the direction
#' of \eqn{y} (a pure vertical outlier) can distort the fitted model just as
#' much as an outlier in the predictor space, and the resulting weights
#' should downweight both kinds -- the weighted correlation matrix used for
#' the Condition Number itself is then formed from the predictors alone,
#' using those joint-space weights.
#'
#' @param formula a two-sided formula.
#' @param data a data frame.
#' @param type \code{"O"} (ordinary, default) or \code{"R"}
#'   (robust/weighted).
#' @param cc tuning constant for Tukey's biweight (default 4.685). Only
#'   used when \code{type = "R"}.
#' @param seed optional integer for reproducible direction sampling. Only
#'   used when \code{type = "R"}.
#' @return an object of class \code{"cn"}, a list with \code{$cn} (the
#'   Condition Number, a single number -- the largest Condition Index) and
#'   \code{$ci} (the full named vector of Condition Indices). Both are
#'   printed by default; access either individually via \code{$cn} /
#'   \code{$ci} if only one is needed.
#' @examples
#' cn(mpg ~ disp + hp + wt + drat, data = mtcars)
#' \donttest{
#' ## type = "R" requires the mrfDepth package
#' cn(mpg ~ disp + hp + wt + drat, data = mtcars, type = "R")
#' }
#' @export
cn <- function(formula, data, type = c("O", "R"), cc = 4.685, seed = NULL) {
  type <- match.arg(type)
  formula <- stats::as.formula(formula)
  mf <- stats::model.frame(formula, data = data)
  y <- stats::model.response(mf)
  mm <- stats::model.matrix(formula, data = data)
  if ("(Intercept)" %in% colnames(mm)) mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  if (ncol(mm) < 2) stop("Condition Number requires at least two predictors.", call. = FALSE)

  w <- NULL
  if (type == "O") {
    R <- stats::cor(mm)
  } else {
    if (!is.numeric(y)) {
      stop("The robust Condition Number requires a numeric response (weights ",
           "are computed jointly over the predictor and response space).", call. = FALSE)
    }
    joint <- cbind(mm, .response = y)
    out <- .projection_outlyingness(joint, seed = seed)
    w <- .outlyingness_weight(out, cc = cc)
    if (sum(w) < ncol(mm)) {
      warning("Fewer than p observations received nonzero weight; the robust ",
              "Condition Number may be unreliable.", call. = FALSE)
    }
    R <- .wcor(mm, w)
  }

  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  ev <- ev[ev > .Machine$double.eps]
  if (length(ev) == 0) stop("The correlation matrix is computationally singular.", call. = FALSE)
  ci <- sort(sqrt(max(ev) / ev))
  names(ci) <- paste0("Dim", seq_along(ci))

  structure(list(cn = unname(max(ci)), ci = ci, type = type, weights = w,
                 formula = formula, call = match.call()),
            class = "cn")
}

#' @export
print.cn <- function(x, ...) {
  cat(if (x$type == "O") "Condition Number (ordinary)\n\n" else
        "Condition Number (robust / weighted)\n\n")
  cat(sprintf("Condition Number: %.4f\n\n", x$cn))
  cat("Condition Indices:\n")
  print(round(x$ci, 4))
  if (x$type == "R") {
    cat(sprintf("\n%d of %d observations downweighted to exactly 0\n",
                sum(x$weights == 0), length(x$weights)))
  }
  invisible(x)
}
