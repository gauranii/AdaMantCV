# -----------------------------------------------------------------------------
# File:    tests/testthat/test-power_curve.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Unit tests of R/power_curve.R: one replication rebuilt from its seed, the power of
#          a response with and without signal, the cross-validation curve
#          and the stored curve table.
# Inputs:  None (files are written to temporary folders).
# Outputs: None.
# Usage:   Rscript tests/testthat.R      # all test files
# -----------------------------------------------------------------------------

source("../../R/setup.R")

# Last GCV score of gcv_curve() on the small problem with seed 100, as first
# computed with R 4.6.1; python/tests/test_power_curve.py has the same.
STORED_GCV_LAST <- 8.049197476488950898

#' The small test problem: 12 observations, 3 predictors, a response in the
#' column space of the predictors.
small_problem <- function() {
  X <- seeded_matrix(12, 3, seed = 31)
  list(X = X, y = drop(X %*% c(2, -1, 1)))
}

test_that("power_run is the replication of seed + 10 * run", {
  problem <- small_problem()
  X_svd <- truncated_svd(problem$X, 0)
  lambdas <- lambda_grid(X_svd$d_tilde, 12, 3, 3, 3)
  weights <- ridge_weights(X_svd$d_tilde, lambdas)
  p <- power_run(run = 2, X_svd$u, weights, lambdas, problem$y, sigma = 0.5, n_perms = 9,
                 seed = 100, exact = FALSE)
  expect_equal(dim(p), c(2, 6))
  expect_true(all(p >= 0 & p <= 1))

  set.seed(100 + 10 * 2)
  Y <- problem$y + rnorm(12, mean = 0, sd = 0.5)
  rows <- lapply(1:9, function(b) sample.int(12))
  Z <- cbind(crossprod(X_svd$u, Y), vapply(rows, function(s) crossprod(X_svd$u[s, ], Y),
                                           numeric(3)))
  statistics <- mantel_statistics(Z, weights, lambdas)
  expect_equal(p[1, ], permutation_pvalues(statistics)[1, ])
  expect_equal(p[2, ], gamma_pvalues(statistics)[1, ])
  # permuting the rows of U is permuting the response
  expect_equal(Z[, 2], drop(crossprod(X_svd$u, Y[order(rows[[1]])])))
  expect_equal(power_run(2, X_svd$u, weights, lambdas, problem$y, 0.5, 9, 100, exact = TRUE), p)
})

test_that("power_curve has one row per penalty and the columns of the original files", {
  problem <- small_problem()
  curve <- power_curve(problem$X, problem$y, n_negative = 3, n_positive = 3, sigma = 0.01,
                       n_perms = 19, n_runs = 4, alpha = 0.05, seed = 100, q = 0)
  expect_named(curve, c("lambdas", "ave_logpval", "ave_logpval_MOM", "sd_logpval",
                        "sd_logpval_MOM", "power", "power_MOM", "power_check", "sign"))
  expect_equal(curve$lambdas, lambda_grid(truncated_svd(problem$X, 0)$d_tilde, 12, 3, 3, 3))
  expect_equal(curve$sign, c(1, 1, 2, 3, 3, 3))
  expect_equal(curve$power, curve$power_check)
  # a response in the column space with almost no noise: at lambda = 0 the
  # observed statistic is above every permuted one in every replication
  expect_equal(curve$power[curve$lambdas == 0], 1)
  expect_equal(curve$ave_logpval[curve$lambdas == 0], -log10(sqrt(.Machine$double.eps)))
  expect_equal(curve$sd_logpval[curve$lambdas == 0], 0)
})

test_that("power_curve has no power for pure noise and does not depend on cores", {
  problem <- small_problem()
  curve <- power_curve(problem$X, rep(0, 12), 3, 3, sigma = 1, n_perms = 19, n_runs = 6,
                       alpha = 0.05, seed = 100)
  expect_true(all(curve$power <= 0.5))
  expect_true(all(curve$power >= 0 & curve$power_MOM >= 0 & curve$power_MOM <= 1))
  expect_identical(power_curve(problem$X, rep(0, 12), 3, 3, 1, 19, 6, 0.05, 100, cores = 2),
                   curve)
  expect_equal(power_curve(problem$X, rep(0, 12), 3, 3, 1, 19, 6, 0.05, 100, exact = TRUE), curve)
})

test_that("gcv_curve scores the response plus the noise of its seed", {
  problem <- small_problem()
  curve <- gcv_curve(problem$X, problem$y, 3, 3, sigma = 0.5, seed = 100)
  expect_named(curve, c("lambdas", "gcv", "loo"))
  lambdas <- lambda_grid(truncated_svd(problem$X, 0.5)$d_tilde, 12, 3, 3, 3)
  expect_equal(curve$lambdas, lambdas)
  set.seed(100)
  Y <- problem$y + rnorm(12, mean = 0, sd = 0.5)
  cv <- ridge_cv(problem$X, Y, lambdas, variant = "dec2021")
  expect_equal(curve$gcv, drop(cv$gcv))
  expect_equal(curve$loo, drop(cv$loo))
  expect_equal(curve$gcv[6], STORED_GCV_LAST, tolerance = 1e-10)
})

test_that("curve_table joins both curves and stores their parts", {
  problem <- small_problem()
  results <- new_folder()
  table <- curve_table("small", problem$X, problem$y, SMALL_OPTS, SMALL_SETTINGS, results)
  expect_named(table, CURVE_COLUMNS)
  expect_equal(nrow(table), 6)
  expect_equal(sort(list.files(results)),
               c("gcv_small_grid3x3.csv", "power_small_runs2_perms9_grid3x3.csv"))
  power <- power_curve(problem$X, problem$y, 3, 3, 0.5, 9, 2, 0.05, 8891, 1000)
  expect_equal(table$power, power$power)
  expect_equal(table$gcv, gcv_curve(problem$X, problem$y, 3, 3, 0.5, 8891, 1000)$gcv)
  # the second call reads the stored parts
  expect_output(again <- curve_table("small", problem$X, problem$y, SMALL_OPTS, SMALL_SETTINGS,
                                     results), "reusing")
  expect_equal(again, table)
})
