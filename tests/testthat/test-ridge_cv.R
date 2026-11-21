# -----------------------------------------------------------------------------
# File:    tests/testthat/test-ridge_cv.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Unit tests of R/ridge_cv.R: the ridge operators and the cross-validation scores
#          against direct solves and against leave-one-out refits.
# Inputs:  None.
# Outputs: None.
# Usage:   Rscript tests/testthat.R      # all test files
# -----------------------------------------------------------------------------

source("../../R/ridge_cv.R")

#' Leave-one-out mean squared error of ridge regression by refitting.
#'
#' @param X,y Predictors and response.
#' @param penalty Number added to the diagonal of X'X.
#' @return Mean over i of the squared error of observation i predicted from
#'   the fit without it.
refit_loo <- function(X, y, penalty) {
  errors <- vapply(seq_len(nrow(X)), function(i) {
    beta <- solve(crossprod(X[-i, ]) + diag(penalty, ncol(X)), crossprod(X[-i, ], y[-i]))
    y[i] - sum(X[i, ] * beta)
  }, numeric(1))
  mean(errors^2)
}

test_that("matrix_trace adds the diagonal", {
  expect_equal(matrix_trace(diag(c(1, 2, 3))), 6)
  expect_equal(matrix_trace(matrix(1:4, 2, 2)), 5)
})

test_that("ridge_operator is (X'X + n lambda I)^-1 X' when n > p", {
  X <- seeded_matrix(8, 3, seed = 5)
  op <- ridge_operator(X, t(X), crossprod(X), lambda = 0.5, zero_rule = TRUE)
  inverse <- solve(crossprod(X) + diag(8 * 0.5, 3))
  expect_equal(op$inverse, inverse)
  expect_equal(op$W, inverse %*% t(X))
  expect_equal(op$L, X %*% inverse %*% t(X))
})

test_that("ridge_operator uses XX' when n <= p, and I - inverse at lambda = 0", {
  X <- seeded_matrix(4, 6, seed = 6)
  gram <- tcrossprod(X)
  op <- ridge_operator(X, t(X), gram, lambda = 0.5, zero_rule = TRUE)
  expect_equal(op$inverse, solve(gram + diag(4 * 0.5, 4)))
  expect_equal(op$W, t(X) %*% op$inverse)
  expect_equal(op$L, gram %*% op$inverse)

  zero <- ridge_operator(X, t(X), gram, lambda = 0, zero_rule = TRUE)
  expect_equal(zero$L, diag(4) - solve(gram))
  plain <- ridge_operator(X, t(X), gram, lambda = 0, zero_rule = FALSE)
  expect_equal(plain$L, diag(4))               # X X' (X X')^-1
})

test_that("penalty_scores are the GCV formula and the leave-one-out refits", {
  X <- seeded_matrix(10, 3, seed = 7)
  y <- seeded_matrix(10, 1, seed = 8)
  lambda <- 0.3
  scores <- penalty_scores(X, t(X), crossprod(X), lambda, y, feb2022 = TRUE)
  L <- X %*% solve(crossprod(X) + diag(10 * lambda, 3)) %*% t(X)
  residual <- y - L %*% y
  expect_equal(scores$gcv, mean((residual / (1 - sum(diag(L)) / 10))^2))
  expect_equal(scores$loo, refit_loo(X, y, 10 * lambda))
})

test_that("the residual of a later column is the observed response minus its fit", {
  X <- seeded_matrix(10, 3, seed = 7)
  Y <- seeded_matrix(10, 2, seed = 9)
  scores <- penalty_scores(X, t(X), crossprod(X), 0.3, Y, feb2022 = TRUE)
  L <- X %*% solve(crossprod(X) + diag(3, 3)) %*% t(X)
  residual <- Y[, 1] - L %*% Y[, 2]
  expect_equal(scores$loo[2], mean((residual / (1 - diag(L)))^2))
  expect_length(scores$gcv, 2)
})

test_that("ridge_cv has one row per response and one column per penalty", {
  X <- seeded_matrix(10, 3, seed = 7)
  Y <- seeded_matrix(10, 2, seed = 9)
  lambdas <- c(0.1, 0.3, 2)
  cv <- ridge_cv(X, Y, lambdas)
  expect_equal(dim(cv$gcv), c(2, 3))
  expect_equal(dim(cv$loo), c(2, 3))
  expect_equal(cv$loo[1, ], vapply(lambdas, function(lambda) refit_loo(X, Y[, 1], 10 * lambda),
                                   numeric(1)))
  one <- penalty_scores(X, t(X), crossprod(X), 0.3, Y, feb2022 = TRUE)
  expect_equal(cv$gcv[, 2], one$gcv)
  # the two versions of the original code agree away from lambda = 0
  expect_equal(ridge_cv(X, Y, lambdas, variant = "dec2021"), cv)
  expect_error(ridge_cv(X, Y, lambdas, variant = "other"))
})

test_that("at lambda = 0 with n <= p the fit is exact and the Feb 2022 scores are finite", {
  X <- seeded_matrix(4, 6, seed = 6)
  y <- seeded_matrix(4, 1, seed = 10)
  # the smoother I - (XX')^-1 keeps the denominators away from zero, and the
  # residual of an exact fit is zero up to rounding
  feb <- ridge_cv(X, y, c(0, 1), variant = "feb2022")
  expect_true(all(is.finite(feb$gcv)) && all(is.finite(feb$loo)))
  expect_lt(feb$gcv[1, 1], 1e-20)
  expect_lt(feb$loo[1, 1], 1e-20)
  expect_gt(feb$gcv[1, 2], 1e-3)
  dec <- ridge_cv(X, y, c(0, 1), variant = "dec2021")
  expect_equal(dec$gcv[1, 2], feb$gcv[1, 2])
})
