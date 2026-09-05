#' Condition Number (Belsley 1991): the square root of the ratio of the
#' largest to smallest eigenvalue of the sum-of-squares-and-cross-products
#' matrix for the predictors. Takes the correlation (or weighted
#' correlation) matrix directly, since that ratio is invariant to the
#' column scaling used to build it -- so the same correlation matrix VIF is
#' already computed from can be reused here at no extra cost.
#' @keywords internal
#' @noRd
.condition_number <- function(R) {
  if (any(!is.finite(R))) {
    bad <- colnames(R)[apply(R, 2, function(col) any(!is.finite(col)))]
    warning("Condition Number is undefined: the correlation matrix has ",
            "non-finite entries, most likely because ",
            if (length(bad) > 0) {
              paste0(if (length(bad) > 1) "columns " else "column ",
                     paste(bad, collapse = ", "), " have")
            } else "one or more columns have",
            " zero variance -- constant among the (possibly downweighted) ",
            "observations used. Returning NA.", call. = FALSE)
    return(NA_real_)
  }
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  ev <- ev[ev > .Machine$double.eps]
  if (length(ev) == 0) return(NA_real_)
  sqrt(max(ev) / min(ev))
}

#' @keywords internal
#' @noRd
.tukey_weight <- function(u, cc = 4.685) {
  w <- (1 - (u / cc)^2)^2
  w[abs(u) > cc] <- 0
  w
}

#' Convert raw projection outlyingness into Tukey biweight weights.
#'
#' When \code{engine = "internal"}, the maximum, over many random projection
#' directions, of a robust per-direction z-score is systematically larger
#' Convert raw projection outlyingness into Tukey biweight weights.
#' @keywords internal
#' @noRd
.outlyingness_weight <- function(out, cc = 4.685) {
  .tukey_weight(out, cc = cc)
}

#' Stahel-Donoho / projection outlyingness
#'
#' Computes projection outlyingness via \code{mrfDepth::outlyingness()} --
#' the implementation used throughout Jacob & Varadharajan's published
#' methodology.
#'
#' @param X a numeric matrix (n x p), n >= p.
#' @param seed optional integer for reproducible direction sampling. The
#'   global RNG state is restored on exit.
#' @return a numeric vector of length n with the outlyingness of each row.
#' @keywords internal
#' @noRd
.projection_outlyingness <- function(X, seed = NULL) {
  if (!requireNamespace("mrfDepth", quietly = TRUE)) {
    stop("This package's robust methods require the 'mrfDepth' package. ",
         "Install it with install.packages(\"mrfDepth\").", call. = FALSE)
  }
  X <- as.matrix(X)
  if (!is.null(seed)) {
    # Run the (randomised) outlyingness computation in a local, temporary
    # RNG scope so that a user-supplied seed gives reproducible results
    # without modifying the user's global random-number state. This uses
    # withr::with_seed(), which saves and restores the RNG state internally
    # and never writes to the global environment (CRAN policy compliant).
    out <- withr::with_seed(
      seed,
      mrfDepth::outlyingness(X, X, options = list(type = "Rotation"))
    )
  } else {
    out <- mrfDepth::outlyingness(X, X, options = list(type = "Rotation"))
  }
  as.numeric(out$outlyingnessX)
}

#' @keywords internal
#' @noRd
.wcor <- function(X, w) {
  X <- as.matrix(X)
  if (ncol(X) < 2L) {
    R <- matrix(1, 1L, 1L, dimnames = list(colnames(X), colnames(X)))
    return(R)
  }
  wsum <- sum(w)
  xbar <- colSums(X * w) / wsum
  Xc <- sweep(X, 2, xbar, "-")
  S <- crossprod(Xc, Xc * w)
  d <- sqrt(diag(S))
  R <- S / outer(d, d)
  dimnames(R) <- list(colnames(X), colnames(X))
  R
}

#' Capture per-factor levels from a data frame of predictors, for
#' consistent dummy-encoding of new data at predict() time (mirrors
#' \code{lm}'s own \code{xlevels} mechanism).
#' @keywords internal
#' @noRd
.get_xlevels <- function(Xf) {
  if (ncol(Xf) == 0) return(list())
  lapply(as.data.frame(lapply(Xf, factor)), levels)
}

#' Dummy-encode a factor data frame against a fixed, stored set of levels
#' (standard treatment contrasts, intercept column dropped), so that a
#' small or single-row \code{newdata} always produces the same columns, in
#' the same order, as at fitting time -- even if some levels are absent
#' from that particular batch of new rows.
#' @keywords internal
#' @noRd
.encode_factors <- function(Xf, xlevels) {
  if (ncol(Xf) == 0 || length(xlevels) == 0) return(NULL)
  for (nm in names(Xf)) {
    Xf[[nm]] <- factor(Xf[[nm]], levels = xlevels[[nm]])
  }
  Xf_dum <- stats::model.matrix(~ ., data = Xf)
  Xf_dum[, colnames(Xf_dum) != "(Intercept)", drop = FALSE]
}

#' Build a design matrix for new data from a fitted model's formula and
#' stored factor levels; shared by every \code{predict.*} method in the
#' package.
#' @keywords internal
#' @noRd
.new_design_matrix <- function(formula, newdata, xlevels) {
  rhs <- stats::delete.response(stats::terms(formula))
  mf <- stats::model.frame(rhs, data = newdata, na.action = stats::na.pass)
  is_num <- vapply(mf, is.numeric, logical(1))
  Xn <- if (any(is_num)) as.matrix(mf[, is_num, drop = FALSE]) else NULL
  Xf <- mf[, !is_num, drop = FALSE]
  Xf_dum <- .encode_factors(Xf, xlevels)
  X <- cbind(Xn, Xf_dum)
  colnames(X) <- c(colnames(Xn), colnames(Xf_dum))
  X
}

#' Build a fully mean-centred design matrix (numeric predictors and
#' dummy-encoded factor columns alike) and a mean-centred response.
#'
#' Every predictor column -- numeric or 0/1 dummy -- is centred by its own
#' mean before the raise transform is applied. This is what makes the
#' intercept recovered afterwards (\code{y_mean - drop(means \%*\% beta)})
#' mathematically identical to the intercept of an ordinary \code{lm()} fit
#' on the raw data: centring only the numeric block (and leaving dummy
#' columns raw, as in some early drafts of this method) silently changes the
#' fitted model whenever a factor predictor has an unequal split across
#' levels, so every column is centred here for correctness.
#' @keywords internal
#' @noRd
.raise_design <- function(formula, data) {
  mf <- stats::model.frame(formula, data = data)
  y <- stats::model.response(mf)
  if (!is.numeric(y)) stop("The response variable must be numeric.", call. = FALSE)
  preds <- mf[, -1, drop = FALSE]
  if (ncol(preds) < 1) stop("At least one predictor is required.", call. = FALSE)

  is_num <- vapply(preds, is.numeric, logical(1))
  Xn <- preds[, is_num, drop = FALSE]
  Xf <- preds[, !is_num, drop = FALSE]

  Xn_mat <- if (ncol(Xn) > 0) as.matrix(Xn) else NULL
  xlevels <- .get_xlevels(Xf)
  if (ncol(Xf) > 0) {
    Xf_typed <- as.data.frame(lapply(Xf, factor))
    Xf_dum <- .encode_factors(Xf_typed, xlevels)
  } else {
    Xf_dum <- NULL
  }

  X_raw <- cbind(Xn_mat, Xf_dum)
  if (is.null(X_raw)) stop("Could not build a design matrix from the supplied formula.", call. = FALSE)
  colnames(X_raw) <- c(colnames(Xn_mat), colnames(Xf_dum))

  means <- colMeans(X_raw)
  X <- sweep(X_raw, 2, means, "-")

  y_mean <- mean(y)
  y_c <- y - y_mean

  list(X = X, y = y_c, means = means, y_mean = y_mean,
       numeric_names = colnames(Xn_mat), factor_names = colnames(Xf_dum),
       xlevels = xlevels, n = nrow(X), p = ncol(X),
       formula = stats::formula(mf), X_raw = X_raw)
}

#' VIF (correlation-matrix based) directly from a numeric predictor matrix.
#' @keywords internal
#' @noRd
.vif_from_matrix <- function(X) {
  if (ncol(X) < 2) {
    v <- rep(1, ncol(X))
    names(v) <- colnames(X)
    return(v)
  }
  R <- stats::cor(X)
  if (any(!is.finite(R))) {
    # A constant column (zero variance) makes correlation undefined for it.
    # Report VIF = NA for the affected column(s) rather than letting solve()
    # fail unpredictably on a matrix containing NaN.
    v <- rep(NA_real_, ncol(X))
    names(v) <- colnames(X)
    return(v)
  }
  Rinv <- try(solve(R), silent = TRUE)
  if (inherits(Rinv, "try-error")) {
    v <- rep(NA_real_, ncol(X))
  } else {
    v <- diag(Rinv)
  }
  names(v) <- colnames(X)
  v
}

#' Small self-contained IRLS M-estimator with Tukey's biweight, used as a
#' drop-in replacement for \code{MASS::rlm(..., psi = psi.bisquare)} so that
#' the package has no hard dependency on \pkg{MASS}. Operates on an
#' intercept-free design (the caller is expected to have already centred
#' \code{X} and \code{y}).
#' @keywords internal
#' @noRd
.rlm_bisquare <- function(X, y, cc = 4.685, maxit = 200, tol = 1e-6) {
  n <- nrow(X); p <- ncol(X)
  b <- as.numeric(solve(crossprod(X), crossprod(X, y)))
  for (it in seq_len(maxit)) {
    res <- y - as.numeric(X %*% b)
    s <- stats::mad(res)
    if (s < 1e-12) break
    u <- res / s
    w <- .tukey_weight(u, cc = cc)
    if (sum(w) < p + 1) break
    b_new <- as.numeric(solve(crossprod(X, w * X), crossprod(X, w * y)))
    denom <- ifelse(abs(b) > 1e-8, b, 1)
    delta <- max(abs((b_new - b) / denom))
    b <- b_new
    if (delta < tol) break
  }
  res <- y - as.numeric(X %*% b)
  s <- stats::mad(res)
  list(coefficients = b, residuals = res, scale = s, weights = w, iterations = it)
}

#' Sequential (one-variable-at-a-time) raise, following the raise-parameter
#' selection strategy of Jacob and Varadharajan (2023)
#' <doi:10.13189/ms.2023.110106>, generalised to an arbitrary VIF threshold
#' and operating on an already mean-centred design so that dummy-coded
#' factor predictors are handled consistently with the simultaneous method.
#'
#' @param lambda_max if the automatically-computed raise parameter for a
#'   column would exceed this, the column is dropped from the design
#'   instead (see \code{.svif_raise} for the rationale, which applies
#'   equally here). Also applied to an exact (RSS ~ 0) dependency, which no
#'   finite raise parameter can fix.
#' @param lambda optional named numeric vector of user-supplied raise
#'   parameters for some or all columns; named entries are applied directly
#'   (skipping both the automatic VIF-driven selection and the
#'   \code{lambda_max} check) before the automatic loop runs on whatever
#'   remains.
#' @keywords internal
#' @noRd
.sequential_raise <- function(X, y, threshold = 10, margin = 0.04,
                               lambda_max = 1000, lambda = NULL) {
  Cval <- threshold - margin
  Xcur <- X
  k <- numeric(0)
  raised <- character(0)
  dropped <- character(0)
  manual_names <- if (!is.null(lambda) && !is.null(names(lambda))) names(lambda) else character(0)

  # Exact (machine-precision) linear dependencies must be removed up front.
  # No finite raise parameter can fix them, and they also make the ordinary
  # VIF undefined (a singular correlation matrix), which would otherwise
  # silently stop the loop below before anything is dropped -- leaving a
  # rank-deficient design for the caller to choke on. Rank-revealing QR
  # identifies exactly the aliased columns, the same way lm() does.
  qr_full <- qr(Xcur)
  if (qr_full$rank < ncol(Xcur)) {
    keep_idx <- sort(qr_full$pivot[seq_len(qr_full$rank)])
    aliased <- setdiff(colnames(Xcur), colnames(Xcur)[keep_idx])
    aliased <- setdiff(aliased, manual_names)   # never auto-drop a user-specified column
    if (length(aliased) > 0) {
      dropped <- c(dropped, aliased)
      Xcur <- Xcur[, setdiff(colnames(Xcur), aliased), drop = FALSE]
    }
  }

  for (nm in manual_names) {
    if (!(nm %in% colnames(Xcur))) next
    idx <- match(nm, colnames(Xcur))
    other <- setdiff(seq_len(ncol(Xcur)), idx)
    aux <- stats::lm.fit(Xcur[, other, drop = FALSE], Xcur[, idx])
    kk <- lambda[[nm]]
    Xcur[, idx] <- aux$fitted.values + sqrt(kk) * aux$residuals
    k[nm] <- kk
    raised <- c(raised, nm)
  }

  iter <- 0L
  max_iter <- 20L * ncol(Xcur)
  repeat {
    if (iter >= max_iter || ncol(Xcur) <= 1) break
    v <- .vif_from_matrix(Xcur)
    v[colnames(Xcur) %in% c(manual_names, dropped)] <- -Inf
    if (all(!is.finite(v)) || max(v) < threshold) break
    idx <- which.max(v)
    nm <- colnames(Xcur)[idx]
    other <- setdiff(seq_len(ncol(Xcur)), idx)
    aux <- stats::lm.fit(Xcur[, other, drop = FALSE], Xcur[, idx])
    fit_aux <- aux$fitted.values
    res_aux <- aux$residuals
    ESS <- sum(fit_aux^2)          # mean(fit_aux) == 0 because all columns are centred
    RSS <- sum(res_aux^2)
    if (RSS < 1e-12) {             # exact dependency: no finite raise parameter fixes this
      dropped <- c(dropped, nm)
      Xcur <- Xcur[, -idx, drop = FALSE]
      iter <- iter + 1L
      next
    }
    kk <- (ESS + RSS) / (Cval * RSS)
    if (kk > lambda_max) {
      dropped <- c(dropped, nm)
      Xcur <- Xcur[, -idx, drop = FALSE]
      iter <- iter + 1L
      next
    }
    Xcur[, idx] <- fit_aux + sqrt(kk) * res_aux
    k[nm] <- kk
    raised <- c(raised, nm)
    iter <- iter + 1L
  }
  if (iter >= max_iter) {
    warning("Sequential raise stopped after ", max_iter,
            " iterations without bringing all VIFs below the threshold.", call. = FALSE)
  }
  list(X = Xcur, k = k, raised = raised, dropped = dropped, iterations = iter)
}

#' Raw (uncentred) numeric + dummy-encoded design matrix and response,
#' used as the starting point for weighted (robust) raise regression, where
#' centring must be done with weighted means computed after the outlyingness
#' weights are known.
#' @keywords internal
#' @noRd
.raw_design <- function(formula, data) {
  mf <- stats::model.frame(formula, data = data)
  y <- stats::model.response(mf)
  if (!is.numeric(y)) stop("The response variable must be numeric.", call. = FALSE)
  preds <- mf[, -1, drop = FALSE]
  if (ncol(preds) < 1) stop("At least one predictor is required.", call. = FALSE)

  is_num <- vapply(preds, is.numeric, logical(1))
  Xn <- preds[, is_num, drop = FALSE]
  Xf <- preds[, !is_num, drop = FALSE]
  Xn_mat <- if (ncol(Xn) > 0) as.matrix(Xn) else NULL
  xlevels <- .get_xlevels(Xf)
  if (ncol(Xf) > 0) {
    Xf_typed <- as.data.frame(lapply(Xf, factor))
    Xf_dum <- .encode_factors(Xf_typed, xlevels)
  } else {
    Xf_dum <- NULL
  }
  X <- cbind(Xn_mat, Xf_dum)
  if (is.null(X)) stop("Could not build a design matrix from the supplied formula.", call. = FALSE)
  colnames(X) <- c(colnames(Xn_mat), colnames(Xf_dum))
  list(X = X, y = as.numeric(y), xlevels = xlevels, n = nrow(X), p = ncol(X),
       formula = stats::formula(mf))
}

#' Closed-form raise parameter: exact solution of VIFraise(lambda) = threshold.
#' Returns Inf when tss is so close to 1 that no finite lambda achieves it
#' (the caller compares this against lambda_max to decide whether to drop
#' the column instead of raising it).
#' @keywords internal
#' @noRd
.raise_lambda_closed_form <- function(tss, threshold) {
  if (!is.finite(tss) || tss <= 0) return(1)
  vifraise0 <- (1 + tss) + tss / (1 - tss)
  if (vifraise0 < threshold) return(1)
  denom <- (1 - tss) * (threshold - 1 - tss)
  lam2 <- if (denom > 0) tss / denom - 1 else Inf
  lam2 <- max(lam2, 0)
  sqrt(lam2) + 1
}

#' Literal replica of the original grid search for the raise parameter
#' (Jacob & Varadharajan's published R code), offered as an explicit,
#' opt-in alternative to the closed-form solution for anyone who needs
#' byte-for-byte reproducibility with earlier published results. Numerically
#' equivalent to the closed form up to the grid's own resolution.
#' @keywords internal
#' @noRd
.raise_lambda_grid <- function(tss, threshold, lambda_max, lambda_step) {
  if (!is.finite(tss) || tss <= 0) return(1)
  lambda <- seq(0, lambda_max, by = lambda_step)
  U <- (1 + lambda^2) * (1 - tss^2) + tss
  L <- (1 + lambda^2) * (1 - tss)
  vifraise <- U / L
  sindex <- sum(vifraise >= threshold)
  slambda <- lambda[sindex + 1]
  if (is.na(slambda)) return(Inf)  # never dropped below threshold within [0, lambda_max]
  slambda + 1
}

#' Sequential Variance Inflation Factor + QR-based raising of every flagged
#' column in one pass (Simultaneous Raise Regression, Jacob & Varadharajan
#' 2022 <doi:10.1007/s11135-022-01557-9>).
#'
#' If a column's required raise parameter would exceed \code{lambda_max},
#' the column is not raised to that extreme value; instead it is *dropped*
#' from the design and the whole procedure is recomputed on the remaining
#' columns. A raise parameter needing to be that large signals a predictor
#' carrying essentially no information beyond what the other predictors
#' already supply -- distinguishing it as still "raisable in principle" and
#' forcing an enormous coefficient scale-down is less useful than simply
#' recognising it as redundant and excluding it. This can cascade (dropping
#' one column changes the QR context for the remaining ones), so the
#' recomputation repeats until no remaining column needs more than
#' \code{lambda_max}.
#' @param X numeric design matrix.
#' @param threshold SVIF threshold.
#' @param lambda_max raise parameters above this are treated as "drop the
#'   column" rather than "raise it this much" (default 1000).
#' @param lambda_step grid resolution, only used when \code{raise_method =
#'   "grid"} (default 0.1, matching the original implementation).
#' @param lambda optional named numeric vector of user-supplied raise
#'   parameters (the actual multiplier k applied to the diagonal of R, i.e.
#'   what is reported in \code{$k}/\code{$svif}) for some or all columns;
#'   named entries override the automatic threshold/lambda_max logic
#'   entirely for that column (it is never dropped, regardless of
#'   \code{lambda_max}).
#' @param raise_method "closed_form" (default, exact) or "grid" (literal
#'   replica of the original grid search, for reproducibility).
#' @param dropped internal accumulator of dropped column names across
#'   recursive calls; not set by the user.
#' @keywords internal
#' @noRd
.svif_raise <- function(X, threshold = 10, lambda_max = 1000, lambda_step = 0.1,
                         lambda = NULL, raise_method = c("closed_form", "grid"),
                         dropped = character(0), orig_svif = NULL) {
  raise_method <- match.arg(raise_method)
  qrd <- qr(X)
  Q <- qr.Q(qrd)
  R <- qr.R(qrd)
  d <- ncol(X)
  cn <- colnames(X)

  tss <- numeric(d)
  for (i in seq_len(d)) {
    tss[i] <- sum(R[seq_len(i), i]^2)
    tss[i] <- 1 - (R[i, i]^2) / tss[i]
  }
  svif_now <- 1 / (1 - tss)
  names(svif_now) <- cn
  # svif_before, as returned to the caller, always reflects the *original*
  # (pre-drop) design -- captured only on the outermost call -- so that a
  # dropped variable's SVIF is still visible to explain why it was dropped;
  # svif_now (below) is this recursion level's own SVIF, used only to decide
  # what happens at this level.
  if (is.null(orig_svif)) orig_svif <- svif_now

  k <- numeric(d)
  need_drop <- character(0)
  Rtilde <- R
  for (i in seq_len(d)) {
    nm <- cn[i]
    if (!is.null(lambda) && !is.null(names(lambda)) && nm %in% names(lambda)) {
      k[i] <- lambda[[nm]]
      Rtilde[i, i] <- k[i] * R[i, i]
      next
    }
    slambda <- if (raise_method == "closed_form") {
      .raise_lambda_closed_form(tss[i], threshold)
    } else {
      .raise_lambda_grid(tss[i], threshold, lambda_max, lambda_step)
    }
    if (slambda > lambda_max) {
      need_drop <- c(need_drop, nm)
      next
    }
    k[i] <- slambda
    Rtilde[i, i] <- slambda * R[i, i]
  }

  if (length(need_drop) > 0) {
    if (length(need_drop) >= d) {
      stop("Every remaining predictor (", paste(sQuote(need_drop), collapse = ", "),
           ") would need a raise parameter above lambda_max = ", lambda_max,
           " -- the design cannot be stabilised by dropping columns alone. ",
           "Try a larger lambda_max, supply manual `lambda` values for these ",
           "columns, or remove some predictors yourself.", call. = FALSE)
    }
    keep <- setdiff(cn, need_drop)
    return(.svif_raise(X[, keep, drop = FALSE], threshold = threshold,
                        lambda_max = lambda_max, lambda_step = lambda_step,
                        lambda = lambda, raise_method = raise_method,
                        dropped = c(dropped, need_drop), orig_svif = orig_svif))
  }

  names(k) <- cn
  list(Q = Q, R = R, Rtilde = Rtilde, svif_before = orig_svif, k = k, dropped = dropped)
}
