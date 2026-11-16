# -----------------------------------------------------------------------------
# File:    R/adamant_cv.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Tables 1 to 4, Figure 3.
# Purpose: The adaptive Mantel test with cross-validation (AdaMantCV): the
#          ridge penalty is chosen by GCV or LOOCV for the observed data and
#          for every permutation, and the p-value at the chosen penalty is
#          ranked among the permutations.
# Inputs:  None.
# Outputs: None.
# Usage:   source("R/setup.R"); adamant_cv(X, Y, q = 0.5, 50, 200, 1000, 8891)
# Source:  R_adaMantCV.R (Feb 2022).
# -----------------------------------------------------------------------------

# Percent rank as computed by the C++ "hybrid" evaluation of dplyr 0.8.5, which
# the original runs used for the final adaptive p-value: each distinct value,
# in ascending order, gets the running sum of (number of ties) / (m - 1). The
# sum is built by repeated addition, so a value such as 0.05 can come out one
# rounding step above 50 / 1000; rejections (p <= alpha) depend on this.
percent_rank_cumulative <- function(x) {
  distinct <- sort(unique(x))
  position <- match(x, distinct)
  ties <- tabulate(position, length(distinct))
  cumsum(c(0, ties / (length(x) - 1)))[position]
}

# X: n x p predictors; Y: one response of length n; q: quantile below which
# singular values are set to zero; n_negative, n_positive: sizes of the
# penalty grid; n_perms: permutations; seed: seed of the permutations.
# Row 1 of every per-permutation result is the observed data.
# exact = TRUE computes every statistic in the original operation order (bit
# for bit); FALSE is the faster matrix form (differences near 1e-14).
adamant_cv <- function(X, Y, q, n_negative, n_positive, n_perms, seed,
                       upper = 1000, exact = TRUE) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n <- nrow(X)
  stopifnot(nrow(Y) == n, ncol(Y) == 1)

  X_svd <- truncated_svd(X, q)
  lambdas <- lambda_grid(X_svd$d_tilde, n, ncol(X), n_negative, n_positive, upper)

  set.seed(seed)
  Y_perm <- cbind(Y, vapply(seq_len(n_perms), function(b) Y[sample.int(n), ],
                            numeric(n)))
  Z <- apply(Y_perm, 2, function(y) crossprod(X_svd$u, y))

  statistics <- mantel_statistics(Z, ridge_weights(X_svd$d_tilde, lambdas),
                                  lambdas, exact)
  p_perm <- permutation_pvalues(statistics)
  p_gamma <- gamma_pvalues(statistics)
  cv <- ridge_cv(X, Y_perm, lambdas, variant = "feb2022")

  rows <- seq_len(n_perms + 1)
  best_gcv <- apply(cv$gcv, 1, which.min)
  best_loo <- apply(cv$loo, 1, which.min)
  adaptive <- function(pvalues, best) {
    percent_rank_cumulative(pvalues[cbind(rows, best)])[1]
  }
  list(p_values = c(gcv = adaptive(p_perm, best_gcv),
                    loo = adaptive(p_perm, best_loo),
                    ga_gcv = adaptive(p_gamma, best_gcv),
                    ga_loo = adaptive(p_gamma, best_loo)),
       lambdas = lambdas,
       lambda_gcv = lambdas[best_gcv], lambda_loo = lambdas[best_loo],
       gcv = cv$gcv, loo = cv$loo,
       statistics = statistics, p_perm = p_perm, p_gamma = p_gamma)
}
