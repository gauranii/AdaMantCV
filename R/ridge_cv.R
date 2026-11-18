# -----------------------------------------------------------------------------
# File:    R/ridge_cv.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Methods (all figures and tables).
# Purpose: Generalized and leave-one-out cross-validation scores of ridge
#          regression over a grid of penalties, in the two versions of
#          the original code (Dec 2021 and Feb 2022).
# Inputs:  None.
# Outputs: None.
# Usage:   source("R/setup.R")
# Source:  R_fcns.R (inv_Tikhonov, lambda_optimal; Feb 2022);
#          RPC_fcns.R (lambda_optimal; Dec 2021).
# -----------------------------------------------------------------------------

#' Trace of a matrix, accumulated in a loop in double precision as in the
#' original code.
#'
#' @param A Square numeric matrix.
#' @return The trace.
matrix_trace <- function(A) {
  total <- 0
  for (value in diag(A)) total <- total + value
  total
}

#' Ridge operators for one penalty.
#'
#' @param X Predictors, n x p.
#' @param Xt Transpose of `X`.
#' @param gram X'X when n > p and XX' otherwise.
#' @param lambda Ridge penalty.
#' @param zero_rule TRUE (the Feb 2022 code): for n <= p and lambda = 0 the
#'   smoother is I - inverse; otherwise it is X W.
#' @return List with `W` (maps a response to coefficients), `L` (the smoother
#'   matrix) and `inverse` (the inverse of the penalized Gram matrix).
ridge_operator <- function(X, Xt, gram, lambda, zero_rule) {
  n <- nrow(X)
  inverse <- MASS::ginv(gram + diag(n * lambda, nrow(gram)),
                        tol = sqrt(.Machine$double.eps))
  W <- if (n > ncol(X)) inverse %*% Xt else Xt %*% inverse
  L <- if (zero_rule && n <= ncol(X) && lambda == 0) diag(1, n) - inverse else X %*% W
  list(W = W, L = L, inverse = inverse)
}

#' GCV and LOOCV scores of every response for one penalty.
#'
#' The residual of a column of `Y` is the observed response (the first
#' column) minus the fit of that column.
#'
#' @param X,Xt,gram,lambda As for `ridge_operator()`.
#' @param Y Responses, one column each, observed response first.
#' @param feb2022 TRUE for the Feb 2022 code: `zero_rule` of
#'   `ridge_operator()`, and a GCV score that is NA while its LOOCV partner
#'   is not is computed, with that partner, from the inverse matrix.
#' @return List with `gcv` and `loo`, one score per column of `Y`.
penalty_scores <- function(X, Xt, gram, lambda, Y, feb2022) {
  n <- nrow(X)
  observed <- Y[, 1]
  op <- ridge_operator(X, Xt, gram, lambda, zero_rule = feb2022)
  residual <- observed - X %*% (op$W %*% Y)
  gcv <- colSums((residual / (1 - matrix_trace(op$L) / n))^2) / n
  loo <- colSums((residual / (1 - diag(op$L)))^2) / n

  from_inverse <- feb2022 & is.na(gcv) & !is.na(loo)
  if (any(from_inverse)) {
    numerator <- op$inverse %*% observed
    gcv[from_inverse] <- sum((numerator / (matrix_trace(op$inverse) / n))^2) / n
    loo[from_inverse] <- sum((numerator / diag(op$inverse))^2) / n
  }
  list(gcv = gcv, loo = loo)
}

#' GCV and LOOCV scores of ridge regression over a grid of penalties.
#'
#' @param X Predictors, n x p.
#' @param Y Responses, one column each; the first column is the observed
#'   response.
#' @param lambdas Ridge penalties.
#' @param variant "feb2022" (Tables 1 to 4, Figure 3) or "dec2021" (one
#'   response, Figures 2 and 5 to 8); see `penalty_scores()`.
#' @param cores Number of parallel processes over the penalties.
#' @return List with matrices `gcv` and `loo`: one row per response, one
#'   column per penalty.
ridge_cv <- function(X, Y, lambdas, variant = c("feb2022", "dec2021"), cores = 1) {
  feb2022 <- match.arg(variant) == "feb2022"
  Y <- as.matrix(Y)
  Xt <- t(X)
  gram <- if (nrow(X) > ncol(X)) Xt %*% X else X %*% Xt

  scores <- parallel::mclapply(lambdas, function(lambda) {
    penalty_scores(X, Xt, gram, lambda, Y, feb2022)
  }, mc.cores = cores)
  failed <- vapply(scores, inherits, logical(1), what = "try-error")
  if (any(failed)) stop(sum(failed), " penalties failed: ", scores[failed][[1]])
  column <- function(name) {
    matrix(vapply(scores, `[[`, numeric(ncol(Y)), name), nrow = ncol(Y))
  }
  list(gcv = column("gcv"), loo = column("loo"))
}
