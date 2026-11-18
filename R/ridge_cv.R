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

# Trace accumulated in a loop in double precision, as in the original code.
matrix_trace <- function(A) {
  total <- 0
  for (value in diag(A)) total <- total + value
  total
}

# Ridge operators for one penalty: W maps a response to coefficients and L is
# the smoother matrix. `gram` is X'X when n > p and XX' otherwise.
# zero_rule = TRUE (the Feb 2022 code): for n <= p and lambda = 0 the smoother
# is I - inverse; otherwise it is X W.
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
# result) and every penalty (columns of the result). The first column of Y is
# the observed response; the residual of a column is the observed response
# minus the fit of that column.
# variant = "feb2022" (Tables 1 to 4, Figure 3): ridge_operator() with
#   zero_rule = TRUE, and a GCV score that is NA while its LOOCV partner is
#   not is computed, with that partner, from the inverse matrix.
# variant = "dec2021" (one response, Figures 2 and 5 to 8): neither of the two.
ridge_cv <- function(X, Y, lambdas, variant = c("feb2022", "dec2021")) {
  feb2022 <- match.arg(variant) == "feb2022"
  Y <- as.matrix(Y)
  n <- nrow(X)
  observed <- Y[, 1]
  gram <- if (n > ncol(X)) t(X) %*% X else X %*% t(X)
  gcv <- loo <- matrix(NA_real_, nrow = ncol(Y), ncol = length(lambdas))

  for (l in seq_along(lambdas)) {
    op <- ridge_operator(X, gram, lambdas[l], zero_rule = feb2022)
    residual <- observed - X %*% (op$W %*% Y)
    gcv[, l] <- colSums((residual / (1 - matrix_trace(op$L) / n))^2) / n
    loo[, l] <- colSums((residual / (1 - diag(op$L)))^2) / n

    from_inverse <- feb2022 & is.na(gcv[, l]) & !is.na(loo[, l])
    if (any(from_inverse)) {
      numerator <- op$inverse %*% observed
      gcv[from_inverse, l] <- sum((numerator / (matrix_trace(op$inverse) / n))^2) / n
      loo[from_inverse, l] <- sum((numerator / diag(op$inverse))^2) / n
    }
  }
  list(gcv = gcv, loo = loo)
}
