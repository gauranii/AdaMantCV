# -----------------------------------------------------------------------------
# File:    R/power_curve.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Figures 2 and 5 to 8.
# Purpose: Power of the ridge-penalized Mantel test and the GCV and LOOCV
#          scores as functions of the ridge penalty, for one response.
# Inputs:  None.
# Outputs: results/power_<tag>_<settings>.csv, results/gcv_<tag>_<grid>.csv
#          (written by curve_table()).
# Usage:   source("R/setup.R"); power_curve(X, y, 500, 2000, 0.5, 1000, 1000, 0.05, 8891)
# Source:  RPC_fcns.R (lambda_values, power_ridge, power_lambda,
#          lambda_optimal; Dec 2021).
# -----------------------------------------------------------------------------

POWER_QUANTILE <- 0.5  # singular values below their median are set to zero

# One replication: noisy response, its statistic and those of n_perms
# permutations. Returns the permutation (row 1) and gamma (row 2) p-values of
# the observed data for every penalty.
# Quirk kept: a permutation permutes the ROWS OF U, not the response.
power_run <- function(run, U, weights, lambdas, y, sigma, n_perms, seed, exact) {
  n <- length(y)
  set.seed(seed + 10 * run)
  Y <- y + rnorm(n, mean = 0, sd = sigma)
  Z <- cbind(crossprod(U, Y),
             vapply(seq_len(n_perms), function(b) drop(crossprod(U[sample.int(n), ], Y)),
                    numeric(ncol(U))))
  statistics <- mantel_statistics(Z, weights, lambdas, exact)
  rbind(observed_pvalues(statistics), gamma_pvalues(statistics, rows = 1))
}

# Power and mean -log10 p-value for every penalty over n_runs replications,
# for the permutation test and its gamma approximation (MOM). Columns as in
# the archived CV_ch<k>.csv files, without the two cross-validation columns.
power_curve <- function(X, y, n_negative, n_positive, sigma, n_perms, n_runs,
                        alpha, seed, upper = 1000, cores = 1, exact = FALSE) {
  X_svd <- truncated_svd(X, POWER_QUANTILE)
  lambdas <- lambda_grid(X_svd$d_tilde, nrow(X), ncol(X), n_negative, n_positive, upper)
  weights <- ridge_weights(X_svd$d_tilde, lambdas)

  runs <- parallel::mclapply(seq_len(n_runs), power_run, U = X_svd$u,
                             weights = weights, lambdas = lambdas, y = y,
                             sigma = sigma, n_perms = n_perms, seed = seed,
                             exact = exact, mc.cores = cores)
  failed <- vapply(runs, inherits, logical(1), what = "try-error")
  if (any(failed)) stop(sum(failed), " replications failed: ", runs[failed][[1]])
  pvals <- t(vapply(runs, function(run) run[1, ], numeric(length(lambdas))))
  pvals_mom <- t(vapply(runs, function(run) run[2, ], numeric(length(lambdas))))

  minus_log10 <- function(p) -log10(p + sqrt(.Machine$double.eps))
  power <- function(p) colSums(p <= alpha, na.rm = TRUE) / n_runs
  data.frame(lambdas = lambdas,
             ave_logpval = apply(minus_log10(pvals), 2, mean, na.rm = TRUE),
             ave_logpval_MOM = apply(minus_log10(pvals_mom), 2, mean, na.rm = TRUE),
             sd_logpval = apply(minus_log10(pvals), 2, sd, na.rm = TRUE),
             sd_logpval_MOM = apply(minus_log10(pvals_mom), 2, sd, na.rm = TRUE),
             power = power(pvals),
             power_MOM = power(pvals_mom),
             power_check = colSums(1 - (pvals > alpha), na.rm = TRUE) / n_runs,
             sign = 2 + sign(lambdas))
}

# GCV and LOOCV scores of ridge regression of one noisy copy of the response
# for every penalty of the same grid.
gcv_curve <- function(X, y, n_negative, n_positive, sigma, seed, upper = 1000) {
  X_svd <- truncated_svd(X, POWER_QUANTILE)
  lambdas <- lambda_grid(X_svd$d_tilde, nrow(X), ncol(X), n_negative, n_positive, upper)
  set.seed(seed)
  Y <- y + rnorm(length(y), mean = 0, sd = sigma)
  cv <- ridge_cv(X, Y, lambdas, variant = "dec2021")
  data.frame(lambdas = lambdas, gcv = drop(cv$gcv), loo = drop(cv$loo))
}

# Power and cross-validation curves of one response as one table (the columns
# of the archived CV_ch<k>.csv). Both parts are stored in results/ under
# `tag` and the settings in `opts` (runs, perms, negative, positive, cores),
# so figures that share a response share the work.
curve_table <- function(tag, X, y, opts) {
  power <- cached_csv(
    results_path(sprintf("power_%s_%s.csv", tag, settings_tag(opts, opts$runs))),
    function() power_curve(X, y, opts$negative, opts$positive, SIGMA, opts$perms,
                           opts$runs, ALPHA, SEED, UPPER, opts$cores))
  cv <- cached_csv(
    results_path(sprintf("gcv_%s_grid%dx%d.csv", tag, opts$negative, opts$positive)),
    function() gcv_curve(X, y, opts$negative, opts$positive, SIGMA, SEED, UPPER))
  cbind(power, cv[c("gcv", "loo")])
}
