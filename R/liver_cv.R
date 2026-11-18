# -----------------------------------------------------------------------------
# File:    R/liver_cv.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Figure 1.
# Purpose: K-fold cross-validated MSE of ridge regression (glmnet) over a
#          grid of penalties for several responses, and its plot.
# Inputs:  None.
# Outputs: None.
# Usage:   source("R/setup.R"); cv_curves(X, Y, 10^seq(4, -2.5, by = -0.1), 10)
# Source:  fig1/CV Ridge.R.
# -----------------------------------------------------------------------------

#' Draw and discard fold assignments.
#'
#' Advances the random number stream by what `count` calls of
#' `glmnet::cv.glmnet()` draw (one fold assignment each).
#'
#' @param count Number of fold assignments.
#' @param n Number of observations.
#' @param n_folds Number of folds.
#' @return NULL, invisibly.
skip_fold_assignments <- function(count, n, n_folds) {
  for (i in seq_len(count)) sample(rep(seq_len(n_folds), length = n))
  invisible(NULL)
}

#' Cross-validated MSE of ridge regression for every response.
#'
#' The fits draw their folds from the current random number stream, one
#' response after the other.
#'
#' @param X Predictors.
#' @param Y Responses, one column each.
#' @param lambdas Decreasing penalties.
#' @param n_folds Number of folds.
#' @return Long data frame with columns `variable` ("Variable<j>"),
#'   `log_lambda`, `cvm` (cross-validated MSE), `at_minimum` (TRUE where
#'   `cvm` is smallest for the variable) and `log_lambda_min`.
cv_curves <- function(X, Y, lambdas, n_folds) {
  do.call(rbind, lapply(seq_len(ncol(Y)), function(j) {
    fit <- glmnet::cv.glmnet(X, Y[, j], alpha = 0, nfolds = n_folds, lambda = lambdas)
    data.frame(variable = paste0("Variable", j), log_lambda = log(lambdas),
               cvm = fit$cvm, at_minimum = fit$cvm == min(fit$cvm),
               log_lambda_min = log(fit$lambda.min))
  }))
}

#' Plot of the cross-validated MSE against log(lambda), one line per variable
#' and a point at each minimum.
#'
#' @param curves Output of `cv_curves()`.
#' @param y_max Upper limit of the MSE axis.
#' @return A ggplot object.
cv_plot <- function(curves, y_max) {
  curves$variable <- factor(curves$variable, levels = unique(curves$variable))
  ggplot(curves, aes(x = log_lambda, y = cvm, color = variable)) +
    geom_line(linewidth = 1) + theme_bw() +
    geom_point(data = curves[curves$at_minimum, ], color = "blue") +
    xlab("log(lambda)") + ylab("Cross-validated MSE") +
    theme(legend.position = "none") +
    ylim(NA, y_max)
}

#' The minimum of every curve.
#'
#' @param curves Output of `cv_curves()`.
#' @param panel Label of the panel.
#' @return Data frame with columns `panel`, `variable`, `log_lambda_min`,
#'   `min_cvm`.
cv_minima <- function(curves, panel) {
  lowest <- curves[curves$at_minimum & curves$log_lambda == curves$log_lambda_min, ]
  data.frame(panel = panel, variable = lowest$variable,
             log_lambda_min = lowest$log_lambda_min, min_cvm = lowest$cvm)
}
