# -----------------------------------------------------------------------------
# File:    R/ridge_cv.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Methods (all figures and tables).
# Purpose: Generalized and leave-one-out cross-validation scores of ridge
#          regression over a grid of penalties, as computed in the two
#          generations of the original code.
# Inputs:  None.
# Outputs: None.
# Usage:   source("R/setup.R")
# Source:  R_fcns.R (inv_Tikhonov, lambda_optimal; Feb 2022);
#          RPC_fcns.R (lambda_optimal; Dec 2021).
# -----------------------------------------------------------------------------

# Trace by plain double precision accumulation, as in the original loop
# (sum() accumulates in extended precision on some platforms).
matrix_trace <- function(A) {
  total <- 0
  for (value in diag(A)) total <- total + value
  total
}

# Ridge operators for one penalty: W maps a response to coefficients and L is
# the smoother matrix. `gram` is X'X when n > p and XX' otherwise.
# Quirk kept from the Feb 2022 code (`zero_rule = TRUE`): for n <= p and
# lambda = 0 the smoother is taken as I - inverse instead of X W.
ridge_operator <- function(X, gram, lambda, zero_rule) {
  n <- nrow(X)
  inverse <- MASS::ginv(gram + diag(n * lambda, nrow(gram)),
                        tol = sqrt(.Machine$double.eps))
  if (n > ncol(X)) {
    W <- inverse %*% t(X)
  } else {
    W <- t(X) %*% inverse
  }
  L <- if (zero_rule && n <= ncol(X) && lambda == 0) diag(1, n) - inverse else X %*% W
  list(W = W, L = L, inverse = inverse)
}

# GCV and LOOCV scores for every response in the columns of Y (rows of the
# result) and every penalty (columns of the result).
# variant = "feb2022" keeps three quirks of the code behind Tables 1 to 4:
#   - the residual of every column uses the FIRST column of Y (the observed
#     response), not the column that was fitted;
#   - the lambda = 0 rule of ridge_operator();
#   - a score that is NA while its LOOCV partner is not is replaced by a score
#     built from the inverse matrix (the original test
#     `is.na(gcv || is.na(loo))` is true in exactly that case).
# variant = "dec2021" (one response, Figs 2 and 5 to 8) has none of the last two.
# variant = "paper" follows the paper: each column is scored with its own
#   residual, the smoother is always X W, and where 1 - tr(L) / n is zero
#   (lambda = 0 with rank(X) = n) the scores are their limits as lambda -> 0
#   (Hastie et al. 2019), built from the inverse matrix.
ridge_cv <- function(X, Y, lambdas, variant = c("feb2022", "dec2021", "paper")) {
  variant <- match.arg(variant)
  feb2022 <- variant == "feb2022"
  Y <- as.matrix(Y)
  n <- nrow(X)
  observed <- Y[, 1]
  gram <- if (n > ncol(X)) t(X) %*% X else X %*% t(X)
  gcv <- loo <- matrix(NA_real_, nrow = ncol(Y), ncol = length(lambdas))

  for (l in seq_along(lambdas)) {
    op <- ridge_operator(X, gram, lambdas[l], zero_rule = feb2022)
    if (variant == "paper") {
      df_left <- 1 - matrix_trace(op$L) / n
      if (abs(df_left) < TIE_TOLERANCE) {
        numerator <- op$inverse %*% Y
        gcv[, l] <- colSums((numerator / (matrix_trace(op$inverse) / n))^2) / n
        loo[, l] <- colSums((numerator / diag(op$inverse))^2) / n
      } else {
        residual <- Y - X %*% (op$W %*% Y)
        gcv[, l] <- colSums((residual / df_left)^2) / n
        loo[, l] <- colSums((residual / (1 - diag(op$L)))^2) / n
      }
      next
    }
    residual <- observed - X %*% (op$W %*% Y)
    gcv[, l] <- colSums((residual / (1 - matrix_trace(op$L) / n))^2) / n
    loo[, l] <- colSums((residual / (1 - diag(op$L)))^2) / n

    replace <- feb2022 & is.na(gcv[, l]) & !is.na(loo[, l])
    if (any(replace)) {
      numerator <- op$inverse %*% observed
      gcv[replace, l] <- sum((numerator / (matrix_trace(op$inverse) / n))^2) / n
      loo[replace, l] <- sum((numerator / diag(op$inverse))^2) / n
    }
  }
  list(gcv = gcv, loo = loo)
}
