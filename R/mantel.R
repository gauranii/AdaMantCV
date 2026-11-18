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

# Fraction of values strictly below each value (dplyr::percent_rank).
percent_rank <- function(x) {
  (match(x, sort(x)) - 1) / (length(x) - 1)
}

# SVD of X with the singular values below their q-th quantile set to zero.
truncated_svd <- function(X, q) {
  X_svd <- svd(X)
  d <- X_svd$d
  list(u = X_svd$u, d = d, d_tilde = d * (d >= quantile(d, probs = q)[[1]]))
}

# Grid of ridge penalties: n_negative - 1 values from the lower bound up to
# -sqrt(eps), then zero, then n_positive values from sqrt(eps) to `upper`.
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

# Weights max(0, d^2 / (lambda + d^2)), one column per penalty. At lambda = 0
# every weight is one, also for the singular values that were set to zero.
ridge_weights <- function(d_tilde, lambdas) {
  d2 <- d_tilde^2
  weights <- pmax(d2 / outer(d2, lambdas, "+"), 0)
  weights[, lambdas == 0] <- 1
  weights
}

# Statistic T[b, l] = sum(weights[, l] * Z[, b]^2) for every data set b (the
# columns of Z, observed data first) and every penalty l.
# exact = FALSE uses one matrix product. exact = TRUE computes every statistic
# as its own inner product, in the operation order of the original code; the
# column of lambda = 0 is always computed that way.
mantel_statistics <- function(Z, weights, lambdas, exact = FALSE) {
  statistics <- crossprod(Z^2, weights)
  original_order <- if (exact) seq_along(lambdas) else which(lambdas == 0)
  for (l in original_order) {
    statistics[, l] <- apply(Z, 2, function(z) crossprod(weights[, l], z^2))
  }
  statistics
}

# Permutation p-value of every data set within each penalty: one minus the
# fraction of statistics strictly below its own.
permutation_pvalues <- function(statistics) {
  1 - apply(statistics, 2, percent_rank)
}

# Same as permutation_pvalues(statistics)[1, ], for the observed data only.
observed_pvalues <- function(statistics) {
  observed <- rep(statistics[1, ], each = nrow(statistics))
  1 - colSums(statistics < observed) / (nrow(statistics) - 1)
}

# Gamma approximation of the permutation distribution of each column. Shape
# and rate come from the mean and variance of the permuted statistics (rows
# 2 onwards); p-values are returned for the rows in `rows`.
# mean() and var() are called with na.rm = TRUE, as in the original code.
gamma_pvalues <- function(statistics, rows = seq_len(nrow(statistics))) {
  center <- apply(statistics[-1, , drop = FALSE], 2, mean, na.rm = TRUE)
  spread <- apply(statistics[-1, , drop = FALSE], 2, var, na.rm = TRUE)
  pgamma(statistics[rows, , drop = FALSE],
         shape = rep(center^2 / spread, each = length(rows)),
         rate = rep(center / spread, each = length(rows)), lower.tail = FALSE)
}
