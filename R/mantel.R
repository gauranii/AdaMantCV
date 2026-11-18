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
  (rank(x, ties.method = "min") - 1) / (length(x) - 1)
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
# At lambda = 0 the statistic equals ||Y||^2 for every permutation in exact
# arithmetic when U is square, so its ranks come from rounding. That column is
# therefore computed one permutation at a time, exactly as the original code
# did. `exact = TRUE` computes every column that way.
mantel_statistics <- function(Z, weights, lambdas, exact = FALSE) {
  statistics <- crossprod(Z^2, weights)
  original_form <- if (exact) seq_along(lambdas) else which(lambdas == 0)
  for (l in original_form) {
    statistics[, l] <- apply(Z, 2, function(z) crossprod(weights[, l], z^2))
  }
  statistics
}

# Relative tolerance below which two statistics count as equal (a tie) in the
# "paper" method. Rounding differences are near 1e-15 relative; real ones are
# far larger than 1e-8.
TIE_TOLERANCE <- sqrt(.Machine$double.eps)

# Permutation p-value of every data set within each penalty.
# method = "original": 1 - percent rank, as the original code (ties at
#   lambda = 0 are broken by rounding).
# method = "paper": the formula of the paper, P(b) = #{c : T(c) >= T(b)} / (B + 1)
#   over c = 0, ..., B, with values within TIE_TOLERANCE counted as ties.
permutation_pvalues <- function(statistics, method = c("original", "paper")) {
  if (match.arg(method) == "original") return(1 - apply(statistics, 2, percent_rank))
  apply(statistics, 2, function(t) {
    sorted <- sort(t)
    below <- findInterval(t - TIE_TOLERANCE * abs(t), sorted, left.open = TRUE)
    (length(t) - below) / length(t)
  })
}

# Same as permutation_pvalues(statistics, method)[1, ], for the observed data only.
observed_pvalues <- function(statistics, method = c("original", "paper")) {
  observed <- rep(statistics[1, ], each = nrow(statistics))
  if (match.arg(method) == "paper") {
    return(colSums(statistics >= observed - TIE_TOLERANCE * abs(observed)) / nrow(statistics))
  }
  1 - colSums(statistics < observed) / (nrow(statistics) - 1)
}

# Gamma approximation of the permutation distribution of each column. Shape
# and rate come from the mean and variance of the permuted statistics (rows
# 2 onwards); p-values are returned for the rows in `rows`.
# na.rm = TRUE is kept from the original although nothing is missing: var()
# takes a different code path with it and rounds differently.
# method = "paper": where the permuted statistics are all equal (relative
# spread below TIE_TOLERANCE, as at lambda = 0 when n < p) the null
# distribution is a point mass and P(T >= t) = 1.
gamma_pvalues <- function(statistics, rows = seq_len(nrow(statistics)),
                          method = c("original", "paper")) {
  center <- apply(statistics[-1, , drop = FALSE], 2, mean, na.rm = TRUE)
  spread <- apply(statistics[-1, , drop = FALSE], 2, var, na.rm = TRUE)
  p <- pgamma(statistics[rows, , drop = FALSE],
              shape = rep(center^2 / spread, each = length(rows)),
              rate = rep(center / spread, each = length(rows)), lower.tail = FALSE)
  if (match.arg(method) == "paper") {
    p[, sqrt(spread) <= TIE_TOLERANCE * abs(center)] <- 1
  }
  p
}
