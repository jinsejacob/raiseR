#' Scale a Dataset (Four Conventions)
#'
#' Centers and scales the columns of a numeric matrix or data frame using
#' one of four conventions. Deliberately named \code{scaleDat()} rather
#' than \code{scale()} to avoid masking \code{base::scale()}; it still
#' keeps base R's \code{center}/\code{scale} logical-flag interface,
#' extended with a \code{type} argument for three additional conventions.
#'
#' \describe{
#'   \item{\code{"classical"} (default)}{Subtract the column mean, divide
#'     by the column standard deviation -- i.e. \code{base::scale()}'s own
#'     algorithm.}
#'   \item{\code{"weighted"}}{Subtract a weighted mean, divide by a
#'     weighted standard deviation, with weights from Tukey's biweight
#'     function applied to each row's Stahel-Donoho projection
#'     outlyingness -- the same weighting scheme used throughout this
#'     package's robust methods (\code{rvif()}, \code{robRaise()},
#'     \code{robRidge()}, \code{robLiu()}).}
#'   \item{\code{"median"}}{Subtract the column median, divide by the
#'     column MADN -- the median absolute deviation scaled by 1.4826 so
#'     that it estimates the standard deviation consistently under
#'     normality; this is \code{stats::mad()}'s default.}
#'   \item{\code{"range"}}{Subtract the column minimum and divide by the
#'     column range (max - min), giving values in [0, 1].}
#' }
#'
#' For every type, the base R \code{center}/\code{scale} logical flags are
#' still honoured: \code{center = FALSE} skips the subtraction step (using
#' 0), and \code{scale = FALSE} skips the division step (using 1).
#'
#' @param x a numeric matrix or data frame (or an object coercible to a
#'   numeric matrix).
#' @param center logical; if \code{TRUE} (default), subtract the
#'   type-appropriate center statistic; if \code{FALSE}, don't center.
#' @param scale logical; if \code{TRUE} (default), divide by the
#'   type-appropriate scale statistic; if \code{FALSE}, don't scale.
#' @param type one of \code{"classical"} (default), \code{"weighted"},
#'   \code{"median"}, \code{"range"}.
#' @param cc tuning constant for Tukey's biweight (default 4.685). Only
#'   used when \code{type = "weighted"}.
#' @param seed optional integer for reproducible direction sampling. Only
#'   used when \code{type = "weighted"}.
#' @return the scaled matrix, with attributes \code{"scaled:center"} and
#'   \code{"scaled:scale"} recording what was subtracted/divided by
#'   (mirroring \code{base::scale()}'s own attributes), plus, for
#'   \code{type = "weighted"}, a \code{"weights"} attribute holding the
#'   per-observation weights used.
#' @examples
#' scaleDat(mtcars)                    # same algorithm as base::scale(mtcars)
#' scaleDat(mtcars, scale = FALSE)     # same algorithm as base::scale(mtcars, scale = FALSE)
#' scaleDat(mtcars, type = "median")   # median/MADN
#' scaleDat(mtcars, type = "range")    # min-max to [0, 1]
#' \donttest{
#' scaleDat(mtcars, type = "weighted") # Stahel-Donoho + Tukey biweight
#' }
#' @export
scaleDat <- function(x, center = TRUE, scale = TRUE,
                   type = c("classical", "weighted", "median", "range"),
                   cc = 4.685, seed = NULL) {
  type <- match.arg(type)
  nm <- colnames(x)
  xm <- as.matrix(x)
  storage.mode(xm) <- "double"
  if (is.null(colnames(xm))) colnames(xm) <- nm

  w <- NULL
  if (type == "classical") {
    ctr <- colMeans(xm)
    scl <- apply(xm, 2, stats::sd)
  } else if (type == "median") {
    ctr <- apply(xm, 2, stats::median)
    scl <- apply(xm, 2, stats::mad)                # MADN by default (constant = 1.4826)
    zero_mad <- scl < 1e-12
    if (any(zero_mad)) {
      # MAD is exactly 0 whenever >= 50% of a column's values are tied at
      # the median -- very common for binary/indicator or other
      # low-cardinality columns (e.g. 0/1 dummies). Fall back to the
      # interquartile range (also robust, but tolerant of more ties) for
      # just those columns, and only if that too is ~0 leave the column
      # unscaled, rather than silently dividing by zero into Inf/NaN.
      iqr <- apply(xm[, zero_mad, drop = FALSE], 2, stats::IQR) / 1.349
      still_zero <- iqr < 1e-12
      scl[zero_mad] <- iqr
      scl[zero_mad][still_zero] <- 1
      if (any(!still_zero)) {
        warning("Column(s) ", paste(names(scl)[zero_mad][!still_zero], collapse = ", "),
                " have zero MAD (>= 50% of values tied at the median -- common ",
                "for binary/indicator columns); using the interquartile range ",
                "instead for just those columns.", call. = FALSE)
      }
      if (any(still_zero)) {
        warning("Column(s) ", paste(names(scl)[zero_mad][still_zero], collapse = ", "),
                " have both zero MAD and zero IQR (heavily constant values); ",
                "left unscaled (divisor = 1) for those columns rather than ",
                "producing Inf/NaN.", call. = FALSE)
      }
    }
  } else if (type == "range") {
    mins <- apply(xm, 2, min)
    maxs <- apply(xm, 2, max)
    ctr <- mins
    scl <- maxs - mins
  } else {                                          # "weighted"
    out <- .projection_outlyingness(xm, seed = seed)
    w <- .outlyingness_weight(out, cc = cc)
    wsum <- sum(w)
    ctr <- colSums(xm * w) / wsum
    Xc <- sweep(xm, 2, ctr, "-")
    scl <- sqrt(colSums(w * Xc^2) / (wsum - 1))
  }

  if (!center) ctr <- stats::setNames(rep(0, ncol(xm)), colnames(xm))
  if (!scale)  scl <- stats::setNames(rep(1, ncol(xm)), colnames(xm))

  if (scale && any(scl < 1e-12 | !is.finite(scl))) {
    warning("At least one column has ~zero (or non-finite) scale; ",
            "resulting values may be Inf/NaN for that column.", call. = FALSE)
  }
  out_mat <- sweep(sweep(xm, 2, ctr, "-"), 2, scl, "/")
  attr(out_mat, "scaled:center") <- ctr
  attr(out_mat, "scaled:scale") <- scl
  if (!is.null(w)) attr(out_mat, "weights") <- w
  out_mat
}
