# -----------------------------------------------------------------------------
# File:    python/adamant/adamant_cv.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Tables 1 to 4, Figure 3.
# Purpose: The adaptive Mantel test with cross-validation (AdaMantCV): the
#          ridge penalty is chosen by GCV or LOOCV for the observed data and
#          for every permutation, and the p-value at the chosen penalty is
#          ranked among the permutations.
# Inputs:  None.
# Outputs: None.
# Usage:   from adamant.adamant_cv import adamant_cv
#          adamant_cv(X, Y, 0.5, 50, 200, 1000, 8891)
# Source:  R/adamant_cv.R.
# -----------------------------------------------------------------------------
import numpy as np

from . import rbase, rng
from .mantel import (gamma_pvalues, lambda_grid, mantel_statistics, permutation_pvalues,
                     ridge_weights, truncated_svd)
from .ridge_cv import ridge_cv

# Names of the four tests: penalty chosen by GCV or LOOCV, permutation
# p-values or their gamma approximation (ga).
TESTS = ["gcv", "loo", "ga_gcv", "ga_loo"]


def percent_rank_cumulative(x):
    """Percent rank as the compiled code of dplyr 0.8.5 (used by the original
    runs) computes it.

    Each distinct value, in ascending order, gets the running sum of (number
    of values equal to it) / (m - 1).

    Parameters: x (numeric vector of length m).
    Returns: vector of the same length.
    """
    x = np.asarray(x, dtype=np.float64)
    distinct, position, counts = np.unique(x, return_inverse=True, return_counts=True)
    return np.cumsum(np.concatenate(([0.0], counts / (len(x) - 1))))[position]


def adaptive_pvalue(chosen_pvalues):
    """Adaptive p-value of the observed data.

    Parameters: chosen_pvalues (p-values of all data sets, each at its own
    chosen penalty; the observed data first).
    Returns: the percent rank of the first element among all elements.
    """
    return float(percent_rank_cumulative(chosen_pvalues)[0])


def adamant_cv(X, Y, q, n_negative, n_positive, n_perms, seed, upper=1000, exact=True,
               observed_lambda0=None, cores=1):
    """Adaptive Mantel test with cross-validation.

    Parameters: X (predictors, n x p), Y (one response of length n), q
    (quantile below which singular values are set to zero), n_negative,
    n_positive (sizes of the penalty grid), n_perms (number of
    permutations), seed (seed of the permutations), upper (largest penalty),
    exact (True computes every statistic in the operation order of the
    original code; False is the faster matrix form), observed_lambda0
    (optional dict with perms_below and gamma, values of the observed data
    at lambda = 0 saved by the original run: the number of permuted
    statistics below the observed one and the gamma p-value; when given, the
    two p-values of the observed data at lambda = 0 are taken from them),
    cores (number of parallel processes over the penalties in the
    cross-validation).
    Returns: dict with `p_values` (the four adaptive p-values, named as
    `TESTS`), `lambdas` (the grid), `lambda_gcv` and `lambda_loo` (chosen
    penalties), `gcv` and `loo` (scores), `statistics`, `p_perm` and
    `p_gamma`. Row or element 1 of every per-data-set result is the observed
    data, followed by the permutations.
    """
    X = rbase.as_matrix(X)
    Y = rbase.as_matrix(Y)
    n = X.shape[0]
    assert Y.shape == (n, 1)

    X_svd = truncated_svd(X, q)
    lambdas = lambda_grid(X_svd["d_tilde"], n, X.shape[1], n_negative, n_positive, upper)

    rng.set_seed(seed)
    permutations = rng.permutation_rows(n, n_perms)
    Y_perm = np.asfortranarray(np.column_stack([Y[:, 0], Y[:, 0][permutations].T]))
    Z = rbase.crossprod_columns(X_svd["u"], Y_perm)

    statistics = mantel_statistics(Z, ridge_weights(X_svd["d_tilde"], lambdas), lambdas, exact)
    p_perm = permutation_pvalues(statistics)
    p_gamma = gamma_pvalues(statistics)
    if observed_lambda0 is not None:
        p_perm[0, lambdas == 0] = 1 - observed_lambda0["perms_below"] / n_perms
        p_gamma[0, lambdas == 0] = observed_lambda0["gamma"]
    cv = ridge_cv(X, Y_perm, lambdas, variant="feb2022", cores=cores)

    rows = np.arange(n_perms + 1)
    best_gcv = rbase.which_min_rows(cv["gcv"])
    best_loo = rbase.which_min_rows(cv["loo"])
    return {"p_values": {"gcv": adaptive_pvalue(p_perm[rows, best_gcv]),
                         "loo": adaptive_pvalue(p_perm[rows, best_loo]),
                         "ga_gcv": adaptive_pvalue(p_gamma[rows, best_gcv]),
                         "ga_loo": adaptive_pvalue(p_gamma[rows, best_loo])},
            "lambdas": lambdas,
            "lambda_gcv": lambdas[best_gcv], "lambda_loo": lambdas[best_loo],
            "gcv": cv["gcv"], "loo": cv["loo"],
            "statistics": statistics, "p_perm": p_perm, "p_gamma": p_gamma}
