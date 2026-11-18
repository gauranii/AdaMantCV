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

# Percent rank as the compiled code of dplyr 0.8.5 (used by the original runs)
# computes it: each distinct value, in ascending order, gets the running sum
# of (number of values equal to it) / (m - 1).
percent_rank_cumulative <- function(x) {
  distinct <- sort(unique(x))
  position <- match(x, distinct)
  counts <- tabulate(position, length(distinct))
  cumsum(c(0, counts / (length(x) - 1)))[position]
}

# Adaptive p-value: the percent rank of the observed data (first element)
# among the p-values of all data sets, each taken at its own chosen penalty.
adaptive_pvalue <- function(chosen_pvalues) percent_rank_cumulative(chosen_pvalues)[1]

# X: n x p predictors; Y: one response of length n; q: quantile below which
# singular values are set to zero; n_negative, n_positive: sizes of the
# penalty grid; n_perms: permutations; seed: seed of the permutations.
# Row 1 of every per-permutation result is the observed data.
# exact = TRUE computes every statistic in the operation order of the
# original code; FALSE is the faster matrix form.
# observed_lambda0: optional c(perms_below = , gamma = ), values of the
#   observed data at lambda = 0 saved by the original run: the number of
#   permuted statistics below the observed one and the gamma p-value. When
#   given, the two p-values of the observed data at lambda = 0 are taken from
#   them.
adamant_cv <- function(X, Y, q, n_negative, n_positive, n_perms, seed,
                       upper = 1000, exact = TRUE, observed_lambda0 = NULL) {
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
  if (!is.null(observed_lambda0)) {
    p_perm[1, lambdas == 0] <- 1 - observed_lambda0[["perms_below"]] / n_perms
    p_gamma[1, lambdas == 0] <- observed_lambda0[["gamma"]]
  }
  cv <- ridge_cv(X, Y_perm, lambdas, variant = "feb2022")

  rows <- seq_len(n_perms + 1)
  best_gcv <- apply(cv$gcv, 1, which.min)
  best_loo <- apply(cv$loo, 1, which.min)
  adaptive <- function(pvalues, best) adaptive_pvalue(pvalues[cbind(rows, best)])
  list(p_values = c(gcv = adaptive(p_perm, best_gcv),
                    loo = adaptive(p_perm, best_loo),
                    ga_gcv = adaptive(p_gamma, best_gcv),
                    ga_loo = adaptive(p_gamma, best_loo)),
       lambdas = lambdas,
       lambda_gcv = lambdas[best_gcv], lambda_loo = lambdas[best_loo],
       gcv = cv$gcv, loo = cv$loo,
       statistics = statistics, p_perm = p_perm, p_gamma = p_gamma)
}
