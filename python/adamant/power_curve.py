# -----------------------------------------------------------------------------
# File:    python/adamant/power_curve.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Figures 2 and 5 to 8.
# Purpose: Power of the ridge-penalized Mantel test and the GCV and LOOCV
#          scores as functions of the ridge penalty, for one response.
# Inputs:  None.
# Outputs: results/power_<tag>_<settings>.csv, results/gcv_<tag>_<grid>.csv
#          (written by curve_table()).
# Usage:   from adamant.power_curve import power_curve
#          power_curve(X, y, 500, 2000, 0.5, 1000, 1000, 0.05, 8891)
# Source:  R/power_curve.R.
# -----------------------------------------------------------------------------
import math

import numpy as np
import pandas as pd
from numba import njit

from . import rbase, rng
from .mantel import (gamma_pvalues, lambda_grid, mantel_statistics, observed_pvalues,
                     ridge_weights, truncated_svd)
from .ridge_cv import ridge_cv
from .rng import FUSED, muladd
from .setup import cached_csv, output_file, parallel_map, settings_tag

# Columns of a curve table: the penalties, the mean -log10 p-value and the
# power of the permutation test and of its gamma approximation (MOM), and the
# GCV and LOOCV scores.
CURVE_COLUMNS = ["lambdas", "ave_logpval", "ave_logpval_MOM", "power", "power_MOM",
                 "gcv", "loo"]


def power_run(run, U, weights, lambdas, y, sigma, n_perms, seed, exact):
    """One replication of the power study.

    A noisy response, its statistic and those of `n_perms` permutations. A
    permutation permutes the rows of U, as in the 2021 code, which is the
    same as permuting the response.

    Parameters: run (replication number; the random stream starts at seed
    `seed + 10 * run`), U (left singular vectors of the predictors),
    weights, lambdas, exact (as for `mantel_statistics()`), y (response
    without noise), sigma (sd of the noise), n_perms (number of
    permutations), seed (base seed).
    Returns: matrix with two rows, the permutation and the gamma p-values of
    the observed data, and one column per penalty.
    """
    n = len(y)
    rng.set_seed(seed + 10 * run)
    Y = y + rng.rnorm(n, 0, sigma)
    Z = np.empty((U.shape[1], n_perms + 1), order="F")
    Z[:, 0] = rbase.crossprod(U, Y)[:, 0]
    for b in range(n_perms):
        Z[:, b + 1] = rbase.crossprod(np.asfortranarray(U[rng.sample_int(n) - 1, :]), Y)[:, 0]
    statistics = mantel_statistics(Z, weights, lambdas, exact)
    return np.vstack([observed_pvalues(statistics), gamma_pvalues(statistics, rows=[0])])


@njit(cache=True)
def _minus_log10_moments(p, shift, fused):
    """Mean and sd over the rows of -log10(p + shift), per column, without
    missing values, where R's long double is a double; `fused` adds each
    squared deviation with one rounding."""
    centre = np.empty(p.shape[1])
    spread = np.empty(p.shape[1])
    values = np.empty(p.shape[0])
    for j in range(p.shape[1]):
        n = 0
        for i in range(p.shape[0]):
            if not np.isnan(p[i, j]):
                values[n] = -math.log10(p[i, j] + shift)
                n += 1
        if n == 0:
            centre[j] = np.nan
            spread[j] = np.nan
            continue
        s = 0.0
        for i in range(n):
            s += values[i]
        s /= n
        if np.isfinite(s):
            t = 0.0
            for i in range(n):
                t += values[i] - s
            s += t / n
        total = 0.0
        for i in range(n):
            deviation = values[i] - s
            if fused:
                total = muladd(deviation, deviation, total)
            else:
                total += deviation * deviation
        centre[j] = s
        spread[j] = math.sqrt(total / (n - 1)) if n > 1 else np.nan
    return centre, spread


@njit(cache=True)
def _minus_log10(p, shift):
    """-log10(p + shift) of every entry; NaN stays NaN."""
    out = np.empty(p.shape)
    for j in range(p.shape[1]):
        for i in range(p.shape[0]):
            out[i, j] = np.nan if np.isnan(p[i, j]) else -math.log10(p[i, j] + shift)
    return out


def minus_log10_moments(p, shift):
    """Mean and sd over the rows of -log10(p + shift), per column, without
    missing values, as R's mean() and sd() compute them.

    Parameters: p (matrix of p-values), shift (added before the logarithm).
    Returns: means, sds (float arrays).
    """
    p = np.ascontiguousarray(p, dtype=np.float64)
    if not rbase.WIDE:
        return _minus_log10_moments(p, shift, FUSED)
    centre, spread = rbase.column_moments(_minus_log10(p, shift))
    return centre, np.sqrt(spread)


def power_curve(X, y, n_negative, n_positive, sigma, n_perms, n_runs, alpha, seed,
                upper=1000, cores=1, exact=False, q=0.5):
    """Power and mean -log10 p-value of the Mantel test for every penalty.

    Each replication seeds itself, so the result does not depend on `cores`.

    Parameters: X (predictors, n x p), y (response without noise),
    n_negative, n_positive (sizes of the penalty grid), sigma (sd of the
    noise added in every replication), n_perms (permutations per
    replication), n_runs (number of replications), alpha (significance
    level), seed (base seed), upper (largest penalty), cores (number of
    parallel processes over the replications), exact (as for
    `mantel_statistics()`), q (quantile below which singular values are set
    to zero).
    Returns: DataFrame with one row per penalty and the columns of the
    CV_ch<k>.csv files of the original runs, without the two
    cross-validation columns: lambdas, mean and sd of -log10 p for the
    permutation test and its gamma approximation (MOM), power, power_MOM,
    power_check and sign.
    """
    X = rbase.as_matrix(X)
    y = np.asarray(y, dtype=np.float64)
    X_svd = truncated_svd(X, q)
    lambdas = lambda_grid(X_svd["d_tilde"], X.shape[0], X.shape[1], n_negative, n_positive,
                          upper)
    weights = ridge_weights(X_svd["d_tilde"], lambdas)

    runs = parallel_map(lambda run: power_run(run, X_svd["u"], weights, lambdas, y, sigma,
                                              n_perms, seed, exact),
                        range(1, n_runs + 1), cores)
    pvals = np.vstack([run[0] for run in runs])
    pvals_mom = np.vstack([run[1] for run in runs])

    with np.errstate(invalid="ignore"):
        def power(p):
            return (p <= alpha).sum(axis=0) / n_runs

        ave, sd = minus_log10_moments(pvals, rbase.SQRT_EPS)
        ave_mom, sd_mom = minus_log10_moments(pvals_mom, rbase.SQRT_EPS)
        return pd.DataFrame({
            "lambdas": lambdas,
            "ave_logpval": ave,
            "ave_logpval_MOM": ave_mom,
            "sd_logpval": sd,
            "sd_logpval_MOM": sd_mom,
            "power": power(pvals),
            "power_MOM": power(pvals_mom),
            "power_check": np.nansum(1 - (pvals > alpha), axis=0) / n_runs,
            "sign": 2 + np.sign(lambdas)})


def gcv_curve(X, y, n_negative, n_positive, sigma, seed, upper=1000, q=0.5):
    """GCV and LOOCV scores of one noisy copy of the response for every
    penalty.

    Parameters: X, y, n_negative, n_positive, sigma, seed, upper, q (as for
    `power_curve()`; the noise is drawn with seed `seed`).
    Returns: DataFrame with columns lambdas, gcv, loo.
    """
    X = rbase.as_matrix(X)
    y = np.asarray(y, dtype=np.float64)
    X_svd = truncated_svd(X, q)
    lambdas = lambda_grid(X_svd["d_tilde"], X.shape[0], X.shape[1], n_negative, n_positive,
                          upper)
    rng.set_seed(seed)
    Y = y + rng.rnorm(len(y), 0, sigma)
    cv = ridge_cv(X, Y, lambdas, variant="dec2021")
    return pd.DataFrame({"lambdas": lambdas, "gcv": cv["gcv"][0, :], "loo": cv["loo"][0, :]})


def curve_table(tag, X, y, opts, settings, results_dir):
    """Power and cross-validation curves of one response as one table.

    Both parts are stored in `results_dir` under `tag` and the settings in
    `opts`, so figures that share a response share the work.

    Parameters: tag (label of the response in the file names, for example
    "eeg184"), X, y (predictors and response), opts (options: runs, perms,
    negative, positive, cores), settings (settings of the paper: sigma,
    alpha, seed, upper), results_dir (folder of the stored parts).
    Returns: DataFrame with the columns `CURVE_COLUMNS`, one row per penalty.
    """
    power = cached_csv(
        output_file(results_dir, f"power_{tag}_{settings_tag(opts, opts['runs'])}.csv"),
        lambda: power_curve(X, y, opts["negative"], opts["positive"], settings["sigma"],
                            opts["perms"], opts["runs"], settings["alpha"], settings["seed"],
                            settings["upper"], opts["cores"]))
    cv = cached_csv(
        output_file(results_dir, f"gcv_{tag}_grid{opts['negative']}x{opts['positive']}.csv"),
        lambda: gcv_curve(X, y, opts["negative"], opts["positive"], settings["sigma"],
                          settings["seed"], settings["upper"]))
    return pd.concat([power, cv[["gcv", "loo"]]], axis=1)[CURVE_COLUMNS]
