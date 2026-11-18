# -----------------------------------------------------------------------------
# File:    R/mantel.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Methods (all figures and tables).
# Purpose: Building blocks of the ridge-penalized Mantel test: truncated
#          singular values, the grid of ridge penalties, the test statistic for
#          every penalty, and its permutation and gamma p-values.
# Inputs:  None.
# Outputs: None.
# Usage:   source("R/setup.R")
# Source:  R_fcns.R, R_adaMantCV.R (Feb 2022); RPC_fcns.R (Dec 2021).
# -----------------------------------------------------------------------------

#' Fraction of values strictly below each value (as `dplyr::percent_rank()`).
#'
#' @param x Numeric vector.
#' @return Numeric vector of the same length, in [0, 1].
percent_rank <- function(x) {
  (match(x, sort(x)) - 1) / (length(x) - 1)
}

#' SVD of a matrix with its small singular values set to zero.
#'
#' @param X Numeric matrix.
#' @param q Quantile of the singular values below which they are set to zero.
#' @return List with `u` (left singular vectors), `d` (singular values) and
#'   `d_tilde` (the truncated singular values).
truncated_svd <- function(X, q) {
  X_svd <- svd(X)
  d <- X_svd$d
  list(u = X_svd$u, d = d, d_tilde = d * (d >= quantile(d, probs = q)[[1]]))
}

#' Grid of ridge penalties.
#'
#' `n_negative - 1` values from the lower bound up to -sqrt(eps), then zero,
#' then `n_positive` values from sqrt(eps) to `upper`.
#'
#' @param d_tilde Truncated singular values of the predictors.
#' @param n,p Numbers of rows and columns of the predictors.
#' @param n_negative,n_positive Sizes of the negative and positive parts.
#' @param upper Largest penalty.
#' @return Increasing numeric vector of penalties.
lambda_grid <- function(d_tilde, n, p, n_negative, n_positive, upper = 1000) {
  eps <- sqrt(.Machine$double.eps)
  lower <- -(sqrt(p / n) - 1) * (sqrt(p / n) - 1) * min(d_tilde[d_tilde != 0])
  if (lower >= 0) {
    return(sort(c(0, seq(lower, upper, length.out = n_positive))))
  }
  negative <- seq(min(lower, -eps), max(lower, -eps), length.out = n_negative - 1)
  positive <- seq(min(upper, eps), max(upper, eps), length.out = n_positive)
  sort(c(negative, 0, positive))
}

#' Ridge weights max(0, d^2 / (lambda + d^2)).
#'
#' At lambda = 0 every weight is one, also for the singular values that were
#' set to zero.
#'
#' @param d_tilde Truncated singular values.
#' @param lambdas Ridge penalties.
#' @return Matrix with one row per singular value and one column per penalty.
ridge_weights <- function(d_tilde, lambdas) {
  d2 <- d_tilde^2
  weights <- pmax(d2 / outer(d2, lambdas, "+"), 0)
  weights[, lambdas == 0] <- 1
  weights
}

#' Mantel statistics T[b, l] = sum(weights[, l] * Z[, b]^2).
#'
#' @param Z Matrix U'Y with one column per data set, observed data first.
#' @param weights Ridge weights, one column per penalty.
#' @param lambdas Ridge penalties.
#' @param exact FALSE uses one matrix product. TRUE computes every statistic
#'   as its own inner product, in the operation order of the original code;
#'   the column of lambda = 0 is always computed that way.
#' @return Matrix with one row per data set and one column per penalty.
mantel_statistics <- function(Z, weights, lambdas, exact = FALSE) {
  Z2 <- Z^2
  statistics <- crossprod(Z2, weights)
  data_sets <- seq_len(ncol(Z2))
  original_order <- if (exact) seq_along(lambdas) else which(lambdas == 0)
  for (l in original_order) {
    w <- weights[, l]
    statistics[, l] <- vapply(data_sets, function(b) drop(crossprod(w, Z2[, b])), numeric(1))
  }
  statistics
}

#' Permutation p-value of every data set within each penalty.
#'
#' @param statistics Matrix of statistics, one row per data set.
#' @return Matrix of the same shape: one minus the fraction of statistics of
#'   the column strictly below each entry.
permutation_pvalues <- function(statistics) {
  1 - apply(statistics, 2, percent_rank)
}

#' Permutation p-values of the observed data only.
#'
#' @param statistics Matrix of statistics, observed data in row 1.
#' @return `permutation_pvalues(statistics)[1, ]`, one value per penalty.
observed_pvalues <- function(statistics) {
  observed <- rep(statistics[1, ], each = nrow(statistics))
  1 - colSums(statistics < observed) / (nrow(statistics) - 1)
}

#' Gamma approximation of the permutation distribution within each penalty.
#'
#' Shape and rate come from the mean and variance of the permuted statistics
#' (rows 2 onwards), computed with `na.rm = TRUE` as in the original code.
#'
#' @param statistics Matrix of statistics, observed data in row 1.
#' @param rows Rows for which p-values are returned.
#' @return Matrix of upper tail probabilities, one row per element of `rows`.
gamma_pvalues <- function(statistics, rows = seq_len(nrow(statistics))) {
  permuted <- statistics[-1, , drop = FALSE]
  center <- apply(permuted, 2, mean, na.rm = TRUE)
  spread <- apply(permuted, 2, var, na.rm = TRUE)
  pgamma(statistics[rows, , drop = FALSE],
         shape = rep(center^2 / spread, each = length(rows)),
         rate = rep(center / spread, each = length(rows)), lower.tail = FALSE)
}
