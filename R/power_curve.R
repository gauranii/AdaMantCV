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

# Columns of a curve table: the penalties, the mean -log10 p-value and the
# power of the permutation test and of its gamma approximation (MOM), and the
# GCV and LOOCV scores.
CURVE_COLUMNS <- c("lambdas", "ave_logpval", "ave_logpval_MOM", "power", "power_MOM",
                   "gcv", "loo")

#' One replication of the power study.
#'
#' A noisy response, its statistic and those of `n_perms` permutations. A
#' permutation permutes the rows of U, as in the 2021 code, which is the same
#' as permuting the response: crossprod(U[s, ], Y) = U' Y[order(s)], and
#' order(s) is again a uniform random permutation.
#'
#' @param run Replication number; the random stream starts at seed
#'   `seed + 10 * run`.
#' @param U Left singular vectors of the predictors.
#' @param weights,lambdas,exact As for `mantel_statistics()`.
#' @param y Response without noise.
#' @param sigma Sd of the noise.
#' @param n_perms Number of permutations.
#' @param seed Base seed.
#' @return Matrix with two rows, the permutation and the gamma p-values of
#'   the observed data, and one column per penalty.
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

#' Power and mean -log10 p-value of the Mantel test for every penalty.
#'
#' Each replication seeds itself, so the result does not depend on `cores`.
#'
#' @param X Predictors, n x p.
#' @param y Response without noise.
#' @param n_negative,n_positive Sizes of the penalty grid.
#' @param sigma Sd of the noise added in every replication.
#' @param n_perms Permutations per replication.
#' @param n_runs Number of replications.
#' @param alpha Significance level.
#' @param seed Base seed.
#' @param upper Largest penalty.
#' @param cores Number of parallel processes over the replications.
#' @param exact As for `mantel_statistics()`.
#' @param q Quantile below which singular values are set to zero.
#' @return Data frame with one row per penalty and the columns of the
#'   CV_ch<k>.csv files of the original runs, without the two
#'   cross-validation columns: `lambdas`, mean and sd of -log10 p for the
#'   permutation test and its gamma approximation (MOM), `power`,
#'   `power_MOM`, `power_check` and `sign`.
power_curve <- function(X, y, n_negative, n_positive, sigma, n_perms, n_runs,
                        alpha, seed, upper = 1000, cores = 1, exact = FALSE, q = 0.5) {
  X_svd <- truncated_svd(X, q)
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

#' GCV and LOOCV scores of one noisy copy of the response for every penalty.
#'
#' @param X,y,n_negative,n_positive,sigma,seed,upper,q As for `power_curve()`;
#'   the noise is drawn with seed `seed`.
#' @return Data frame with columns `lambdas`, `gcv`, `loo`.
gcv_curve <- function(X, y, n_negative, n_positive, sigma, seed, upper = 1000, q = 0.5) {
  X_svd <- truncated_svd(X, q)
  lambdas <- lambda_grid(X_svd$d_tilde, nrow(X), ncol(X), n_negative, n_positive, upper)
  set.seed(seed)
  Y <- y + rnorm(length(y), mean = 0, sd = sigma)
  cv <- ridge_cv(X, Y, lambdas, variant = "dec2021")
  data.frame(lambdas = lambdas, gcv = drop(cv$gcv), loo = drop(cv$loo))
}

#' Power and cross-validation curves of one response as one table.
#'
#' Both parts are stored in `results_dir` under `tag` and the settings in
#' `opts`, so figures that share a response share the work.
#'
#' @param tag Label of the response in the file names (for example "eeg184").
#' @param X,y Predictors and response.
#' @param opts Options: `runs`, `perms`, `negative`, `positive`, `cores`.
#' @param settings Settings of the paper (`sigma`, `alpha`, `seed`, `upper`).
#' @param results_dir Folder of the stored parts.
#' @return Data frame with the columns `CURVE_COLUMNS`, one row per penalty.
curve_table <- function(tag, X, y, opts, settings, results_dir) {
  power <- cached_csv(
    output_file(results_dir, sprintf("power_%s_%s.csv", tag, settings_tag(opts, opts$runs))),
    function() power_curve(X, y, opts$negative, opts$positive, settings$sigma, opts$perms,
                           opts$runs, settings$alpha, settings$seed, settings$upper,
                           opts$cores))
  cv <- cached_csv(
    output_file(results_dir, sprintf("gcv_%s_grid%dx%d.csv", tag, opts$negative, opts$positive)),
    function() gcv_curve(X, y, opts$negative, opts$positive, settings$sigma, settings$seed,
                         settings$upper))
  cbind(power, cv[c("gcv", "loo")])[CURVE_COLUMNS]
}
