# -----------------------------------------------------------------------------
# File:    python/adamant/mantel.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Methods (all figures and tables).
# Purpose: Building blocks of the ridge-penalized Mantel test: truncated
#          singular values, the grid of ridge penalties, the test statistic for
#          every penalty, and its permutation and gamma p-values.
# Inputs:  None.
# Outputs: None.
# Usage:   from adamant.mantel import truncated_svd, lambda_grid
# Source:  R/mantel.R.
# -----------------------------------------------------------------------------
import numpy as np
from numba import njit

from . import rbase
from .rng import FUSED, muladd
from .rbase import SQRT_EPS


def percent_rank(x):
    """Fraction of values strictly below each value (as `dplyr::percent_rank()`).

    Parameters: x (numeric vector).
    Returns: vector of the same length, in [0, 1].
    """
    x = np.asarray(x, dtype=np.float64)
    below = np.searchsorted(np.sort(x[~np.isnan(x)]), x, side="left").astype(np.float64)
    below[np.isnan(x)] = np.nan
    return below / (len(x) - 1)


def truncated_svd(X, q):
    """SVD of a matrix with its small singular values set to zero.

    Parameters: X (matrix), q (quantile of the singular values below which
    they are set to zero).
    Returns: dict with `u` (left singular vectors), `d` (singular values)
    and `d_tilde` (the truncated singular values).
    """
    X_svd = rbase.svd(X)
    d = X_svd["d"]
    return {"u": X_svd["u"], "d": d, "d_tilde": d * (d >= rbase.quantile(d, q))}


def lambda_grid(d_tilde, n, p, n_negative, n_positive, upper=1000):
    """Grid of ridge penalties.

    `n_negative - 1` values from the lower bound up to -sqrt(eps), then zero,
    then `n_positive` values from sqrt(eps) to `upper`.

    Parameters: d_tilde (truncated singular values of the predictors), n, p
    (numbers of rows and columns of the predictors), n_negative, n_positive
    (sizes of the negative and positive parts), upper (largest penalty).
    Returns: increasing vector of penalties.
    """
    eps = SQRT_EPS
    lower = -(np.sqrt(p / n) - 1) * (np.sqrt(p / n) - 1) * np.min(d_tilde[d_tilde != 0])
    if lower >= 0:
        return np.sort(np.concatenate(([0.0], rbase.seq_length(lower, upper, n_positive))))
    negative = rbase.seq_length(min(lower, -eps), max(lower, -eps), n_negative - 1)
    positive = rbase.seq_length(min(upper, eps), max(upper, eps), n_positive)
    return np.sort(np.concatenate((negative, [0.0], positive)))


def ridge_weights(d_tilde, lambdas):
    """Ridge weights max(0, d^2 / (lambda + d^2)).

    At lambda = 0 every weight is one, also for the singular values that were
    set to zero.

    Parameters: d_tilde (truncated singular values), lambdas (penalties).
    Returns: matrix with one row per singular value and one column per
    penalty.
    """
    d2 = d_tilde * d_tilde
    with np.errstate(divide="ignore", invalid="ignore"):
        weights = np.maximum(d2[:, None] / (d2[:, None] + lambdas[None, :]), 0.0)
    weights[:, lambdas == 0] = 1.0
    return np.asfortranarray(weights)


def mantel_statistics(Z, weights, lambdas, exact=False):
    """Mantel statistics T[b, l] = sum(weights[, l] * Z[, b]^2).

    Parameters: Z (matrix U'Y with one column per data set, observed data
    first), weights (ridge weights, one column per penalty), lambdas
    (penalties), exact (False uses one matrix product; True computes every
    statistic as its own inner product, in the operation order of the
    original code; the column of lambda = 0 is always computed that way).
    Returns: matrix with one row per data set and one column per penalty.
    """
    Z2 = np.asfortranarray(Z * Z)
    statistics = rbase.crossprod(Z2, weights)
    original_order = range(len(lambdas)) if exact else np.flatnonzero(lambdas == 0)
    for l in original_order:
        statistics[:, l] = rbase.crossprod_columns(weights[:, l], Z2)[0, :]
    return statistics


def permutation_pvalues(statistics):
    """Permutation p-value of every data set within each penalty.

    Parameters: statistics (matrix, one row per data set).
    Returns: matrix of the same shape: one minus the fraction of statistics
    of the column strictly below each entry.
    """
    ranks = np.column_stack([percent_rank(statistics[:, l])
                             for l in range(statistics.shape[1])])
    return 1 - ranks


def observed_pvalues(statistics):
    """Permutation p-values of the observed data only.

    Parameters: statistics (matrix, observed data in row 1).
    Returns: `permutation_pvalues(statistics)[0, :]`, one value per penalty.
    """
    below = (statistics < statistics[0:1, :]).sum(axis=0).astype(np.float64)
    below[np.isnan(statistics).any(axis=0)] = np.nan
    return 1 - below / (statistics.shape[0] - 1)


@njit(cache=True)
def _column_moments(A, fused):
    """Mean and variance of every column, as R's mean() and var() with
    na.rm = TRUE; `fused` adds each squared deviation with one rounding."""
    centre = np.empty(A.shape[1])
    spread = np.empty(A.shape[1])
    for j in range(A.shape[1]):
        n = 0
        s = 0.0
        for i in range(A.shape[0]):
            if not np.isnan(A[i, j]):
                s += A[i, j]
                n += 1
        if n == 0:
            centre[j] = np.nan
            spread[j] = np.nan
            continue
        s /= n
        if np.isfinite(s):
            t = 0.0
            for i in range(A.shape[0]):
                if not np.isnan(A[i, j]):
                    t += A[i, j] - s
            s += t / n
        total = 0.0
        for i in range(A.shape[0]):
            if not np.isnan(A[i, j]):
                deviation = A[i, j] - s
                if fused:
                    total = muladd(deviation, deviation, total)
                else:
                    total += deviation * deviation
        centre[j] = s
        spread[j] = total / (n - 1) if n > 1 else np.nan
    return centre, spread


def gamma_pvalues(statistics, rows=None):
    """Gamma approximation of the permutation distribution within each penalty.

    Shape and rate come from the mean and variance of the permuted statistics
    (rows 2 onwards), computed without missing values as in the original
    code.

    Parameters: statistics (matrix, observed data in row 1), rows (rows for
    which p-values are returned, zero-based; all by default).
    Returns: matrix of upper tail probabilities, one row per element of
    `rows`.
    """
    rows = np.arange(statistics.shape[0]) if rows is None else np.asarray(rows)
    center, spread = _column_moments(np.ascontiguousarray(statistics[1:, :]), FUSED)
    with np.errstate(divide="ignore", invalid="ignore"):
        return rbase.pgamma_upper(statistics[rows, :], (center * center / spread)[None, :],
                                  (center / spread)[None, :])
