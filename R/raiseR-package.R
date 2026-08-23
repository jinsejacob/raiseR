#' raiseR: Ordinary and Robust Raise, Ridge, and Liu Regression
#'
#' Tools for diagnosing and combating multicollinearity in linear
#' regression models, with ordinary and outlier-resistant (robust) variants
#' throughout: the Raise Method (\code{raiseReg}, \code{robRaise}), Ridge
#' Regression (\code{ridgeReg}, \code{robRidge}), Liu Regression
#' (\code{liuReg}, \code{robLiu}), collinearity diagnostics (\code{vif},
#' \code{rvif}, \code{cn}), and a flexible \code{scale} function.
#'
#' @author Jinse Jacob \email{jinsejacob@@hotmail.com}
#'   (ORCID: 0000-0002-0873-4864), Assistant Professor, Department of
#'   Statistics, Govt. Victoria College, Palakkad.
#'
#' @references
#' Jacob, J. and Varadharajan, R. (2023). Simultaneous raise regression: a
#' novel approach to combating collinearity in linear regression models.
#' \emph{Quality & Quantity}, 57(5), 4365-4386.
#'
#' Jacob, J. and Varadharajan, R. (2023). Raise estimation: An alternative
#' approach in the presence of problematic multicollinearity.
#' \emph{Mathematics and Statistics}, 51-64.
#'
#' Jacob, J. and Varadharajan, R. (2024). Robust Variance Inflation Factor:
#' A Promising Approach for Collinearity Diagnostics in the Presence of
#' Outliers. \emph{Sankhya B}, 86(2), 845-871.
#'
#' Hoerl, A. E. and Kennard, R. W. (1970). Ridge regression: Biased
#' estimation for nonorthogonal problems. \emph{Technometrics}, 12(1), 55-67.
#'
#' Liu, K. (1993). A new class of biased estimate in linear regression.
#' \emph{Communications in Statistics - Theory and Methods}, 22(2), 393-402.
#'
#' Yohai, V. J. (1987). High breakdown-point and high efficiency robust
#' estimates for regression. \emph{The Annals of Statistics}, 15(2), 642-656.
#'
#' Kan, B., Alpu, O. and Yazici, B. (2013). Robust ridge and robust Liu
#' estimator for regression based on the LTS estimator. \emph{Journal of
#' Applied Statistics}, 40(3), 644-655.
#'
#' @importFrom stats cooks.distance dfbetas hatvalues df.residual vcov
#' @importFrom stats model.matrix weights
#' @keywords internal
"_PACKAGE"
