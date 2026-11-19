# -----------------------------------------------------------------------------
# File:    python/adamant/ridge_cv.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Methods (all figures and tables).
# Purpose: Generalized and leave-one-out cross-validation scores of ridge
#          regression over a grid of penalties, in the two versions of
#          the original code (Dec 2021 and Feb 2022).
# Inputs:  None.
# Outputs: None.
# Usage:   from adamant.ridge_cv import ridge_cv
# Source:  R/ridge_cv.R.
# -----------------------------------------------------------------------------
import numpy as np

from . import rbase
from .setup import parallel_map


def matrix_trace(A):
    """Trace of a matrix, accumulated in a loop in double precision as in the
    original code.

    Parameters: A (square matrix). Returns: the trace.
    """
    return rbase.r_sum(np.diag(A))


def ridge_operator(X, Xt, gram, lam, zero_rule):
    """Ridge operators for one penalty.

    Parameters: X (predictors, n x p), Xt (transpose of X), gram (X'X when
    n > p and XX' otherwise), lam (ridge penalty), zero_rule (True, the Feb
    2022 code: for n <= p and lambda = 0 the smoother is I - inverse;
    otherwise it is X W).
    Returns: dict with `W` (maps a response to coefficients), `L` (the
    smoother matrix) and `inverse` (the inverse of the penalized Gram
    matrix).
    """
    n, p = X.shape
    inverse = rbase.ginv(gram + np.diag(np.full(gram.shape[0], n * lam)), tol=rbase.SQRT_EPS)
    W = rbase.matprod(inverse, Xt) if n > p else rbase.matprod(Xt, inverse)
    L = np.eye(n) - inverse if zero_rule and n <= p and lam == 0 else rbase.matprod(X, W)
    return {"W": W, "L": L, "inverse": inverse}


def penalty_scores(X, Xt, gram, lam, Y, feb2022):
    """GCV and LOOCV scores of every response for one penalty.

    The residual of a column of `Y` is the observed response (the first
    column) minus the fit of that column.

    Parameters: X, Xt, gram, lam (as for `ridge_operator()`), Y (responses,
    one column each, observed response first), feb2022 (True for the Feb
    2022 code: `zero_rule` of `ridge_operator()`, and a GCV score that is
    missing while its LOOCV partner is not is computed, with that partner,
    from the inverse matrix).
    Returns: dict with `gcv` and `loo`, one score per column of `Y`.
    """
    n = X.shape[0]
    observed = Y[:, 0]
    op = ridge_operator(X, Xt, gram, lam, zero_rule=feb2022)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        residual = observed[:, None] - rbase.matprod(X, rbase.matprod(op["W"], Y))
        gcv = rbase.col_sums(np.square(residual / (1 - matrix_trace(op["L"]) / n))) / n
        loo = rbase.col_sums(np.square(residual / (1 - np.diag(op["L"]))[:, None])) / n

        from_inverse = feb2022 & np.isnan(gcv) & ~np.isnan(loo)
        if from_inverse.any():
            numerator = rbase.matprod(op["inverse"], observed)[:, 0]
            gcv[from_inverse] = rbase.r_sum(
                np.square(numerator / (matrix_trace(op["inverse"]) / n))) / n
            loo[from_inverse] = rbase.r_sum(np.square(numerator / np.diag(op["inverse"]))) / n
    return {"gcv": gcv, "loo": loo}


def ridge_cv(X, Y, lambdas, variant="feb2022", cores=1):
    """GCV and LOOCV scores of ridge regression over a grid of penalties.

    Parameters: X (predictors, n x p), Y (responses, one column each; the
    first column is the observed response), lambdas (ridge penalties),
    variant ("feb2022" for Tables 1 to 4 and Figure 3, or "dec2021" for one
    response, Figures 2 and 5 to 8; see `penalty_scores()`), cores (number
    of parallel processes over the penalties).
    Returns: dict with matrices `gcv` and `loo`: one row per response, one
    column per penalty.
    """
    if variant not in ("feb2022", "dec2021"):
        raise ValueError(f"unknown variant {variant}")
    feb2022 = variant == "feb2022"
    X = rbase.as_matrix(X)
    Y = rbase.as_matrix(Y)
    Xt = np.asfortranarray(X.T)
    gram = rbase.matprod(Xt, X) if X.shape[0] > X.shape[1] else rbase.matprod(X, Xt)

    scores = parallel_map(lambda lam: penalty_scores(X, Xt, gram, lam, Y, feb2022),
                          lambdas, cores)
    return {"gcv": np.column_stack([score["gcv"] for score in scores]),
            "loo": np.column_stack([score["loo"] for score in scores])}
