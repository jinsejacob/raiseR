#' Regression Diagnostics for a raiseReg Object
#'
#' \code{raiseReg} objects support the standard influence-diagnostic suite
#' (\code{hatvalues()}, \code{cooks.distance()}, \code{dfbetas()},
#' \code{covRatio()}) as well as \code{lmtest::bptest()} and
#' \code{car::ncvTest()}. This is possible without approximation because
#' regressing \eqn{y} on the raised design reproduces the ordinary least
#' squares projection of \eqn{y} \emph{exactly} (see \code{?raiseReg}): the
#' hat matrix, residuals and \eqn{\hat\sigma^2} underlying these diagnostics
#' are therefore identical to those of an ordinary \code{lm()} fit on the
#' surviving predictors, and every value returned here is numerically
#' identical to fitting that \code{lm()} directly and calling the same
#' generic on it.
#'
#' These methods are only provided for \code{raiseReg} (the exact, unbiased
#' case). \code{ridgeReg}, \code{liuReg}, and their robust counterparts are
#' deliberately biased estimators without an idempotent hat matrix, so
#' classical influence measures don't have an agreed-on, unambiguous
#' definition for them, and are not provided here; the robust methods
#' (\code{robRaise()}, \code{robRidge()}, \code{robLiu()}) already downweight
#' influential/outlying observations directly via their fitting procedure.
#'
#' @param model,object a \code{raiseReg} object.
#' @param var.formula optional one-sided formula specifying the variables
#'   against which non-constant variance is tested (\code{ncvTest} only);
#'   defaults to the fitted values.
#' @param ... unused (present for generic consistency).
#' @return
#'   \code{hatvalues.raiseReg} returns a named numeric vector of hat
#'   (leverage) values, one per observation.
#'
#'   \code{cooks.distance.raiseReg} returns a named numeric vector of
#'   Cook's distances, one per observation.
#'
#'   \code{dfbetas.raiseReg} returns a numeric matrix of DFBETAS, with one
#'   row per observation and one column per model coefficient.
#'
#'   \code{covRatio} and \code{covRatio.raiseReg} return a named numeric
#'   vector of covariance ratios, one per observation. \code{covRatio} is
#'   the S3 generic; \code{covRatio.default} dispatches to
#'   \code{stats::covratio()}.
#'
#'   \code{ncvTest.raiseReg} returns an object of class \code{"chisqTest"}
#'   (from the \pkg{car} package): a list containing the score-test
#'   statistic \code{ChiSquare}, its degrees of freedom \code{Df}, and the
#'   p-value \code{p} for the test of non-constant error variance.
#' @name raiseReg-diagnostics
NULL

#' @rdname raiseReg-diagnostics
#' @export
hatvalues.raiseReg <- function(model, ...) {
  h <- model$hat
  names(h) <- rownames(model$x)
  h
}

#' @rdname raiseReg-diagnostics
#' @export
cooks.distance.raiseReg <- function(model, ...) {
  h <- model$hat
  e <- model$residuals_inference
  s2 <- model$sigma^2
  p_total <- model$p + 1   # + intercept
  d <- (e^2 * h) / (p_total * s2 * (1 - h)^2)
  names(d) <- rownames(model$x)
  d
}

#' @rdname raiseReg-diagnostics
#' @export
dfbetas.raiseReg <- function(model, ...) {
  X <- model$x
  M <- solve(crossprod(X))
  e <- model$residuals_inference
  h <- model$hat
  s <- model$sigma
  df <- model$df.residual
  # Leave-one-out residual scale, from the standard identity relating the
  # deleted SSE to the retained one: SSE_(i) = SSE - e_i^2/(1-h_i).
  s_i <- sqrt(pmax(((df * s^2) - e^2 / (1 - h)) / (df - 1), 0))

  IF <- M %*% t(X)                          # (p+1) x n
  DFBETA <- sweep(IF, 2, e / (1 - h), "*")  # (p+1) x n
  sqrtMjj <- sqrt(diag(M))
  out <- t(DFBETA) / outer(s_i, sqrtMjj)    # n x (p+1)
  dimnames(out) <- list(rownames(X), colnames(X))
  out
}

#' @rdname raiseReg-diagnostics
#' @export
covRatio <- function(model, ...) UseMethod("covRatio")

#' @rdname raiseReg-diagnostics
#' @export
covRatio.default <- function(model, ...) stats::covratio(model, ...)

#' @rdname raiseReg-diagnostics
#' @export
covRatio.raiseReg <- function(model, ...) {
  X <- model$x
  e <- model$residuals_inference
  h <- model$hat
  s <- model$sigma
  df <- model$df.residual
  p_total <- model$p + 1
  s_i <- sqrt(pmax(((df * s^2) - e^2 / (1 - h)) / (df - 1), 0))
  cr <- (s_i^2 / s^2)^p_total / (1 - h)
  names(cr) <- rownames(X)
  cr
}

#' @rdname raiseReg-diagnostics
#' @export
model.matrix.raiseReg <- function(object, ...) object$x

#' @rdname raiseReg-diagnostics
#' @export
weights.raiseReg <- function(object, ...) NULL

#' @rdname raiseReg-diagnostics
#' @export
df.residual.raiseReg <- function(object, ...) object$df.residual

#' @rdname raiseReg-diagnostics
#' @exportS3Method car::ncvTest
ncvTest.raiseReg <- function(model, var.formula, ...) {
  # A self-contained re-implementation of car::ncvTest()'s algorithm
  # (Breusch-Pagan-style score test for non-constant variance), rather than
  # calling car:::ncvTest.lm() directly, so this package has no hard
  # dependency on car's internals. Registered as an S3 method for car's
  # `ncvTest` generic; requires the car package to be loaded for the
  # generic itself to exist and for `model` to actually be passed to it.
  sumry <- summary(model)
  resid_p <- stats::residuals(model, type = "pearson")
  S.sq <- stats::df.residual(model) * (sumry$sigma)^2 / sum(!is.na(resid_p))
  U <- (resid_p^2) / S.sq
  if (missing(var.formula)) {
    mod <- stats::lm(U ~ model$fitted_inference)
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

#' @rdname raiseReg-diagnostics
#' @export
vcov.raiseReg <- function(object, ...) {
  se <- object$coef_table[, "Std. Error"]
  v <- diag(se^2)
  dimnames(v) <- list(names(se), names(se))
  v
}
