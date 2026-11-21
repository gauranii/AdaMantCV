# -----------------------------------------------------------------------------
# File:    tests/testthat/test-liver_cv.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Unit tests of R/liver_cv.R: the fold assignments drawn and skipped, the
#          cross-validated curves against glmnet with the same folds and
#          against a ridge fit in closed form, their plot and their minima.
# Inputs:  None.
# Outputs: None.
# Usage:   Rscript tests/testthat.R      # all test files
# -----------------------------------------------------------------------------

source("../../R/setup.R")

LAMBDAS <- c(10, 1, 0.1, 0.01)

# cvm of cv_curves() on the small problem with seed 3, as first computed with
# R 4.6.1 and glmnet 5.0; python/tests/test_liver_cv.py has the same.
STORED_CVM <- c(1.9565450234509685, 0.74109508987441752, 0.32795538515868683,
                0.32300771374188475, 1.1202080079946386, 0.61124977338760755,
                0.41355989194038351, 0.41167216021214931)

#' The small test problem: 30 observations, 4 predictors, 2 responses.
small_problem <- function() {
  X <- seeded_matrix(30, 4, seed = 21)
  noise <- seeded_matrix(30, 2, seed = 22)
  list(X = X, Y = cbind(X[, 1] - X[, 2], X[, 3]) + 0.5 * noise)
}

test_that("skip_fold_assignments advances the stream by one fold draw per count", {
  set.seed(5)
  expect_null(skip_fold_assignments(2, n = 10, n_folds = 5))
  after_skipping <- runif(1)
  set.seed(5)
  first <- sample(rep(1:5, length = 10))
  second <- sample(rep(1:5, length = 10))
  expect_equal(after_skipping, runif(1))
  # a fold assignment holds every fold equally often
  expect_equal(as.vector(table(first)), rep(2, 5))
  expect_false(identical(first, second))
  set.seed(5)
  skip_fold_assignments(0, 10, 5)
  expect_equal(sample(rep(1:5, length = 10)), first)
})

test_that("cv_curves is a long table with one curve per response", {
  problem <- small_problem()
  set.seed(3)
  curves <- cv_curves(problem$X, problem$Y, LAMBDAS, n_folds = 5)
  expect_named(curves, c("variable", "log_lambda", "cvm", "at_minimum", "log_lambda_min"))
  expect_equal(curves$variable, rep(c("Variable1", "Variable2"), each = 4))
  expect_equal(curves$log_lambda, rep(log(LAMBDAS), 2))
  expect_true(all(curves$cvm > 0))
  for (variable in unique(curves$variable)) {
    curve <- curves[curves$variable == variable, ]
    expect_equal(sum(curve$at_minimum), 1)
    expect_equal(curve$cvm[curve$at_minimum], min(curve$cvm))
    expect_equal(unique(curve$log_lambda_min), curve$log_lambda[curve$at_minimum])
  }
})

test_that("cv_curves draws the folds of each response from the stream, in turn", {
  problem <- small_problem()
  set.seed(3)
  curves <- cv_curves(problem$X, problem$Y, LAMBDAS, n_folds = 5)
  set.seed(3)
  folds <- lapply(1:2, function(j) sample(rep(1:5, length = 30)))
  for (j in 1:2) {
    fit <- glmnet::cv.glmnet(problem$X, problem$Y[, j], alpha = 0, lambda = LAMBDAS,
                             foldid = folds[[j]])
    expect_equal(curves$cvm[curves$variable == paste0("Variable", j)], fit$cvm)
  }
  set.seed(4)
  expect_false(isTRUE(all.equal(cv_curves(problem$X, problem$Y, LAMBDAS, 5)$cvm, curves$cvm)))
})

test_that("the cross-validated error is that of ridge fits in closed form", {
  # glmnet minimizes sum((y - a - X b)^2) / (2 n) + lambda * |b|^2 / 2 on
  # predictors and response scaled to variance one (divisor n): on centred
  # data b = (X'X + n lambda sd(y) I)^-1 X'y in the units of the response
  problem <- small_problem()
  y <- problem$Y[, 1]
  set.seed(3)
  folds <- sample(rep(1:5, length = 30))
  population_sd <- function(v) sqrt(mean((v - mean(v))^2))
  errors <- vapply(1:5, function(fold) {
    train <- folds != fold
    centre <- colMeans(problem$X[train, ])
    spread <- apply(problem$X[train, ], 2, population_sd)
    Z <- scale(problem$X[train, ], centre, spread)
    b <- solve(crossprod(Z) + diag(sum(train) * LAMBDAS[2] / population_sd(y[train]), 4),
               crossprod(Z, y[train] - mean(y[train])))
    predicted <- mean(y[train]) + scale(problem$X[!train, ], centre, spread) %*% b
    mean((y[!train] - predicted)^2)
  }, numeric(1))
  set.seed(3)
  curves <- cv_curves(problem$X, problem$Y, LAMBDAS, n_folds = 5)
  expect_equal(curves$cvm[2], mean(errors), tolerance = 1e-4)
})

test_that("cv_curves gives the stored values of the small problem", {
  problem <- small_problem()
  set.seed(3)
  curves <- cv_curves(problem$X, problem$Y, LAMBDAS, n_folds = 5)
  expect_equal(curves$cvm, STORED_CVM, tolerance = 1e-10)
})

test_that("cv_plot draws one line per variable and a point at each minimum", {
  problem <- small_problem()
  set.seed(3)
  curves <- cv_curves(problem$X, problem$Y, LAMBDAS, n_folds = 5)
  plot <- cv_plot(curves, y_max = 10)
  expect_s3_class(plot, "ggplot")
  layers <- ggplot_build(plot)$data
  expect_equal(nrow(layers[[1]]), 8)
  expect_equal(length(unique(layers[[1]]$group)), 2)
  expect_equal(sort(layers[[2]]$y), sort(curves$cvm[curves$at_minimum]))
})

test_that("cv_minima has the lowest point of every curve", {
  curves <- data.frame(variable = rep(c("Variable1", "Variable2"), each = 3),
                       log_lambda = rep(c(2, 1, 0), 2), cvm = c(5, 3, 4, 1, 2, 6),
                       at_minimum = c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE),
                       log_lambda_min = rep(c(1, 2), each = 3))
  expect_equal(cv_minima(curves, "A"),
               data.frame(panel = "A", variable = c("Variable1", "Variable2"),
                          log_lambda_min = c(1, 2), min_cvm = c(3, 1)))
})
