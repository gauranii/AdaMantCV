# -----------------------------------------------------------------------------
# File:    tests/testthat/test-mantel.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Unit tests of R/mantel.R: every function on inputs with a known answer
#          (hand-computed values, and the closed form of the ridge statistic).
# Inputs:  None.
# Outputs: None.
# Usage:   Rscript tests/testthat.R      # all test files
# -----------------------------------------------------------------------------

source("../../R/mantel.R")

test_that("percent_rank is the fraction of values strictly below", {
  expect_equal(percent_rank(c(3, 1, 2)), c(1, 0, 0.5))
  expect_equal(percent_rank(c(1, 1, 2)), c(0, 0, 1))
  expect_equal(percent_rank(c(5, 1, 3, 2)), c(3, 0, 2, 1) / 3)
})

test_that("truncated_svd zeroes the singular values below the quantile", {
  X <- diag(c(4, 3, 2, 1))
  expect_equal(truncated_svd(X, 0.5)$d, c(4, 3, 2, 1))
  expect_equal(truncated_svd(X, 0.5)$d_tilde, c(4, 3, 0, 0))
  expect_equal(truncated_svd(X, 0)$d_tilde, c(4, 3, 2, 1))
  expect_equal(truncated_svd(X, 0.25)$d_tilde, c(4, 3, 2, 0))

  X <- seeded_matrix(6, 4, seed = 1)
  X_svd <- truncated_svd(X, 0)
  expect_equal(dim(X_svd$u), c(6, 4))
  expect_equal(crossprod(X_svd$u), diag(4))
  expect_equal(X_svd$d^2, eigen(crossprod(X), symmetric = TRUE)$values)
})

test_that("lambda_grid has the negative part, zero and the positive part", {
  eps <- sqrt(.Machine$double.eps)
  # sqrt(p / n) = 2, smallest kept singular value 3: lower bound -(2 - 1)^2 * 3
  grid <- lambda_grid(c(4, 3, 0, 0), n = 4, p = 16, n_negative = 4, n_positive = 3)
  expect_equal(grid, c(-3, (-3 - eps) / 2, -eps, 0, eps, (eps + 1000) / 2, 1000))
  expect_equal(max(lambda_grid(c(4, 3), 4, 16, 4, 3, upper = 10)), 10)
  # n = p: the lower bound is zero and the grid has no negative part
  expect_equal(lambda_grid(c(2, 1), n = 5, p = 5, n_negative = 4, n_positive = 3),
               c(0, 0, 500, 1000))
})

test_that("ridge_weights are max(0, d^2 / (lambda + d^2)), one at lambda = 0", {
  weights <- ridge_weights(c(2, 1, 0), c(-0.5, 0, 1))
  expect_equal(weights, cbind(c(4 / 3.5, 2, 0), c(1, 1, 1), c(0.8, 0.5, 0)))
  # a penalty below -d^2 gives a negative ratio, which is set to zero
  expect_equal(ridge_weights(c(2, 1), -2)[, 1], c(2, 0))
})

test_that("mantel_statistics are the weighted sums of squares, in both forms", {
  Z <- cbind(c(1, 3), c(2, 4))                 # two data sets
  weights <- cbind(c(1, 1), c(0.5, 0.25))      # two penalties
  expected <- rbind(c(10, 2.75), c(20, 6))     # 1 + 9, 0.5 + 2.25; 4 + 16, 2 + 4
  expect_equal(mantel_statistics(Z, weights, c(0, 1)), expected)
  expect_equal(mantel_statistics(Z, weights, c(0, 1), exact = TRUE), expected)
})

test_that("the statistic is the closed form y'X (X'X + lambda I)^-1 X'y", {
  X <- seeded_matrix(8, 5, seed = 2)
  y <- seeded_matrix(8, 1, seed = 3)
  lambdas <- c(0, 0.5, 3, 1000)
  X_svd <- truncated_svd(X, 0)
  statistics <- mantel_statistics(crossprod(X_svd$u, y), ridge_weights(X_svd$d_tilde, lambdas),
                                  lambdas, exact = TRUE)
  closed_form <- vapply(lambdas, function(lambda) {
    drop(t(y) %*% X %*% solve(crossprod(X) + diag(lambda, 5)) %*% t(X) %*% y)
  }, numeric(1))
  expect_equal(drop(statistics), closed_form)
})

test_that("permutation_pvalues are exact for a small case", {
  statistics <- cbind(c(5, 1, 3, 2), c(1, 2, 3, 4), c(2, 2, 2, 2))
  expected <- cbind(c(0, 1, 1 / 3, 2 / 3), c(1, 2 / 3, 1 / 3, 0), c(1, 1, 1, 1))
  expect_equal(permutation_pvalues(statistics), expected)
  expect_equal(observed_pvalues(statistics), expected[1, ])

  random <- seeded_matrix(20, 6, seed = 4)^2
  p <- permutation_pvalues(random)
  expect_true(all(p >= 0 & p <= 1))
  expect_equal(observed_pvalues(random), p[1, ])
  # without ties every multiple of 1 / 19 appears once per penalty
  expect_equal(sort(p[, 1]), (0:19) / 19)
})

test_that("gamma_pvalues match the gamma with the moments of the permutations", {
  # permuted statistics 1, 2, 3: mean 2, variance 1, so shape 4 and rate 2;
  # P(Gamma(4, 2) > 3) = P(Poisson(6) <= 3) = 61 exp(-6)
  statistics <- cbind(c(3, 1, 2, 3), c(6, 2, 4, 6))
  p <- gamma_pvalues(statistics)
  expect_equal(dim(p), c(4, 2))
  expect_equal(p[1, 1], 61 * exp(-6))
  expect_equal(p[1, 2], 61 * exp(-6))          # the same after doubling every statistic
  expect_equal(p[2, 1], 19 / 3 * exp(-2))      # P(Poisson(2) <= 3)
  expect_equal(gamma_pvalues(statistics, rows = 1), p[1, , drop = FALSE])
  expect_true(all(p >= 0 & p <= 1))
})
