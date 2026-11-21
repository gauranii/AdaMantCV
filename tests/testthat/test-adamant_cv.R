# -----------------------------------------------------------------------------
# File:    tests/testthat/test-adamant_cv.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Unit tests of R/adamant_cv.R: the percent rank of the original runs on hand-computed
#          cases, and the test on a small seeded design (structure, the
#          statistic of the observed data in closed form, the permutations of
#          the seed, stored values).
# Inputs:  None.
# Outputs: None.
# Usage:   Rscript tests/testthat.R      # all test files
# -----------------------------------------------------------------------------

source("../../R/mantel.R")
source("../../R/ridge_cv.R")
source("../../R/adamant_cv.R")

# Values of adamant_cv(X, y, 0.5, 4, 5, 19, seed = 7) on the small problem, as
# first computed with R 4.6.1; python/tests/test_adamant_cv.py has the same.
STORED <- list(lambda_first = -1.0926875988102316, statistic_last = 0.32052880115284371,
               gcv_last = 1.1434087615079229, p_gamma_zero = 4.8335366933517551e-05)

#' The small test problem: 20 observations, 5 predictors, signal in the first.
small_problem <- function() {
  X <- seeded_matrix(20, 5, seed = 11)
  list(X = X, y = X[, 1] + 0.5 * drop(seeded_matrix(20, 1, seed = 12)))
}

test_that("percent_rank_cumulative adds the counts of the distinct values", {
  # distinct values 0.1, 0.2, 0.4 with counts 1, 2, 1 of m - 1 = 3
  expect_equal(percent_rank_cumulative(c(0.2, 0.1, 0.2, 0.4)), c(1 / 3, 0, 1 / 3, 1))
  # without ties it is the usual percent rank
  expect_equal(percent_rank_cumulative(c(3, 1, 2)), percent_rank(c(3, 1, 2)))
})

test_that("adaptive_pvalue is the rank of the first element", {
  expect_equal(adaptive_pvalue(c(0.2, 0.1, 0.2, 0.4)), 1 / 3)
  expect_equal(adaptive_pvalue(c(0.01, 0.5, 0.3)), 0)
  expect_equal(adaptive_pvalue(c(0.9, 0.5, 0.3)), 1)
})

test_that("adamant_cv returns every part with the right shape", {
  problem <- small_problem()
  ada <- adamant_cv(problem$X, problem$y, q = 0, n_negative = 4, n_positive = 5,
                    n_perms = 19, seed = 7)
  expect_named(ada, c("p_values", "lambdas", "lambda_gcv", "lambda_loo", "gcv", "loo",
                      "statistics", "p_perm", "p_gamma"))
  expect_named(ada$p_values, TESTS)
  expect_true(all(ada$p_values >= 0 & ada$p_values <= 1))
  expect_length(ada$lambdas, 9)
  expect_equal(sum(ada$lambdas == 0), 1)
  # one row for the observed data and each of the 19 permutations
  for (part in c("gcv", "loo", "statistics", "p_perm", "p_gamma")) {
    expect_equal(dim(ada[[part]]), c(1 + 19, 9))
  }
  expect_length(ada$lambda_gcv, 1 + 19)
  expect_true(all(ada$lambda_gcv %in% ada$lambdas) && all(ada$lambda_loo %in% ada$lambdas))
  expect_true(all(ada$p_perm >= 0 & ada$p_perm <= 1))
  expect_true(all(ada$p_gamma >= 0 & ada$p_gamma <= 1))
})

test_that("row 1 is the observed data and row 2 the first permutation of the seed", {
  problem <- small_problem()
  ada <- adamant_cv(problem$X, problem$y, 0, 4, 5, 19, seed = 7)
  ridge_statistic <- function(y, lambda) {
    drop(t(y) %*% problem$X %*% solve(crossprod(problem$X) + diag(lambda, 5)) %*%
           t(problem$X) %*% y)
  }
  positive <- which(ada$lambdas > 0)
  expect_equal(ada$statistics[1, positive],
               vapply(ada$lambdas[positive], ridge_statistic, numeric(1), y = problem$y))
  set.seed(7)
  permuted <- problem$y[sample.int(20)]
  expect_equal(ada$statistics[2, positive],
               vapply(ada$lambdas[positive], ridge_statistic, numeric(1), y = permuted))
})

test_that("adamant_cv depends on the seed only, and finds the signal", {
  problem <- small_problem()
  ada <- adamant_cv(problem$X, problem$y, 0, 4, 5, 19, seed = 7)
  runif(3)  # the state of the random stream before the call does not matter
  expect_identical(adamant_cv(problem$X, problem$y, 0, 4, 5, 19, seed = 7), ada)
  expect_equal(ada$statistics, adamant_cv(problem$X, problem$y, 0, 4, 5, 19, seed = 7,
                                          exact = FALSE)$statistics)
  other <- adamant_cv(problem$X, problem$y, 0, 4, 5, 19, seed = 8)
  expect_equal(other$statistics[1, ], ada$statistics[1, ])
  expect_false(isTRUE(all.equal(other$statistics[2, ], ada$statistics[2, ])))
  # the observed statistic is above all 19 permuted ones at lambda = 0
  expect_equal(ada$p_perm[1, ada$lambdas == 0], 0)
  expect_equal(ada$p_values[["gcv"]], 0)
})

test_that("observed_lambda0 replaces the p-values of the observed data at lambda = 0", {
  problem <- small_problem()
  ada <- adamant_cv(problem$X, problem$y, 0, 4, 5, 19, seed = 7,
                    observed_lambda0 = c(perms_below = 15, gamma = 0.123))
  zero <- ada$lambdas == 0
  expect_equal(ada$p_perm[1, zero], 1 - 15 / 19)
  expect_equal(ada$p_gamma[1, zero], 0.123)
})

test_that("adamant_cv gives the stored values of the small problem", {
  problem <- small_problem()
  ada <- adamant_cv(problem$X, problem$y, 0.5, 4, 5, 19, seed = 7)
  expect_equal(ada$lambdas[1], STORED$lambda_first, tolerance = 1e-10)
  expect_equal(ada$statistics[1, 9], STORED$statistic_last, tolerance = 1e-10)
  expect_equal(ada$gcv[1, 9], STORED$gcv_last, tolerance = 1e-10)
  expect_equal(ada$p_gamma[1, 5], STORED$p_gamma_zero, tolerance = 1e-8)
})
