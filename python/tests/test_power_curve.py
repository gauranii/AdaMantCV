# -----------------------------------------------------------------------------
# File:    python/tests/test_power_curve.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Unit tests of python/adamant/power_curve.py: one replication rebuilt from its seed,
#          the mean and sd of -log10(p), the power of a response with and
#          without signal, the cross-validation curve and the stored curve
#          table; the cases of
#          tests/testthat/test-power_curve.R.
# Inputs:  None (files are written to temporary folders).
# Outputs: None.
# Usage:   python -m pytest python/tests/test_power_curve.py     # inside `nix develop`
# -----------------------------------------------------------------------------
import math
import os

import numpy as np
import pandas as pd

from adamant import mantel, rng
from adamant.power_curve import (CURVE_COLUMNS, curve_table, gcv_curve, minus_log10_moments,
                                 power_curve, power_run)
from adamant.ridge_cv import ridge_cv
from inputs import SMALL_OPTS, SMALL_SETTINGS, close, seeded_matrix

# Last GCV score of gcv_curve() on the small problem with seed 100, as first
# computed with R 4.6.1; tests/testthat/test-power_curve.R has the same.
STORED_GCV_LAST = 8.049197476488950898


def small_problem():
    """The small test problem: 12 observations, 3 predictors, a response in
    the column space of the predictors. Returns: (X, y)."""
    X = seeded_matrix(12, 3, seed=31)
    return X, X @ np.array([2.0, -1.0, 1.0])


def test_power_run_is_the_replication_of_seed_plus_10_times_run():
    X, y = small_problem()
    X_svd = mantel.truncated_svd(X, 0)
    lambdas = mantel.lambda_grid(X_svd["d_tilde"], 12, 3, 3, 3)
    weights = mantel.ridge_weights(X_svd["d_tilde"], lambdas)
    p = power_run(2, X_svd["u"], weights, lambdas, y, 0.5, 9, 100, exact=False)
    assert p.shape == (2, 6)
    assert ((p >= 0) & (p <= 1)).all()

    rng.set_seed(100 + 10 * 2)
    Y = y + rng.rnorm(12, 0, 0.5)
    rows = [rng.sample_int(12) - 1 for _ in range(9)]
    Z = np.column_stack([X_svd["u"].T @ Y] + [X_svd["u"][s].T @ Y for s in rows])
    statistics = mantel.mantel_statistics(Z, weights, lambdas)
    close(p[0], mantel.permutation_pvalues(statistics)[0])
    close(p[1], mantel.gamma_pvalues(statistics)[0])
    # permuting the rows of U is permuting the response
    close(Z[:, 1], X_svd["u"].T @ Y[np.argsort(rows[0])])
    close(power_run(2, X_svd["u"], weights, lambdas, y, 0.5, 9, 100, exact=True), p)


def test_minus_log10_moments_are_the_mean_and_sd_of_every_column():
    # columns of -log10(p): (1, 3, 2), (0, 2) with one missing value, a
    # single value 1, and no value at all
    p = np.array([[0.1, 1.0, np.nan, np.nan],
                  [0.001, np.nan, 0.1, np.nan],
                  [0.01, 0.01, np.nan, np.nan]])
    centre, spread = minus_log10_moments(p, 0)
    close(centre[:3], [2, 1, 1])
    close(spread[:2], [1, math.sqrt(2)])
    assert math.isnan(spread[2])                       # sd of one value
    assert math.isnan(centre[3]) and math.isnan(spread[3])
    # the shift is added before the logarithm: p = 0 gives -log10(shift)
    centre, spread = minus_log10_moments(np.zeros((4, 2)), 0.01)
    close(centre, [2, 2])
    close(spread, [0, 0])
    # p-values without a pattern, against numpy
    p = 1 / (1 + np.exp(seeded_matrix(9, 4, seed=32)))
    p[2, 1] = p[5, 1] = p[0, 3] = np.nan
    centre, spread = minus_log10_moments(p, 1e-3)
    values = -np.log10(p + 1e-3)
    close(centre, np.nanmean(values, axis=0))
    close(spread, np.nanstd(values, axis=0, ddof=1))
    # a Fortran-ordered matrix and a list give the same
    again = minus_log10_moments(np.asfortranarray(p), 1e-3)
    assert np.array_equal(again[0], centre) and np.array_equal(again[1], spread)
    close(minus_log10_moments([[0.1], [0.001]], 0)[0], [2])


def test_power_curve_has_one_row_per_penalty_and_the_columns_of_the_original_files():
    X, y = small_problem()
    curve = power_curve(X, y, 3, 3, sigma=0.01, n_perms=19, n_runs=4, alpha=0.05, seed=100, q=0)
    assert list(curve.columns) == ["lambdas", "ave_logpval", "ave_logpval_MOM", "sd_logpval",
                                   "sd_logpval_MOM", "power", "power_MOM", "power_check", "sign"]
    close(curve["lambdas"], mantel.lambda_grid(mantel.truncated_svd(X, 0)["d_tilde"], 12, 3, 3, 3))
    assert curve["sign"].tolist() == [1, 1, 2, 3, 3, 3]
    close(curve["power"], curve["power_check"])
    # a response in the column space with almost no noise: at lambda = 0 the
    # observed statistic is above every permuted one in every replication
    zero = curve[curve["lambdas"] == 0].iloc[0]
    assert zero["power"] == 1
    close(zero["ave_logpval"], -math.log10(math.sqrt(np.finfo(float).eps)))
    assert zero["sd_logpval"] == 0


def test_power_curve_has_no_power_for_pure_noise_and_does_not_depend_on_cores():
    X, _ = small_problem()
    curve = power_curve(X, np.zeros(12), 3, 3, 1, 19, 6, 0.05, 100)
    assert (curve["power"] <= 0.5).all()
    assert ((curve["power"] >= 0) & (curve["power_MOM"] >= 0) & (curve["power_MOM"] <= 1)).all()
    pd.testing.assert_frame_equal(power_curve(X, np.zeros(12), 3, 3, 1, 19, 6, 0.05, 100, cores=2),
                                  curve, check_exact=True)
    exact = power_curve(X, np.zeros(12), 3, 3, 1, 19, 6, 0.05, 100, exact=True)
    close(exact.to_numpy(), curve.to_numpy())


def test_gcv_curve_scores_the_response_plus_the_noise_of_its_seed():
    X, y = small_problem()
    curve = gcv_curve(X, y, 3, 3, sigma=0.5, seed=100)
    assert list(curve.columns) == ["lambdas", "gcv", "loo"]
    lambdas = mantel.lambda_grid(mantel.truncated_svd(X, 0.5)["d_tilde"], 12, 3, 3, 3)
    close(curve["lambdas"], lambdas)
    rng.set_seed(100)
    Y = y + rng.rnorm(12, 0, 0.5)
    cv = ridge_cv(X, Y, lambdas, variant="dec2021")
    close(curve["gcv"], cv["gcv"][0])
    close(curve["loo"], cv["loo"][0])
    close(curve["gcv"][5], STORED_GCV_LAST, rtol=1e-10)


def test_curve_table_joins_both_curves_and_stores_their_parts(tmp_path, capsys):
    X, y = small_problem()
    results = str(tmp_path)
    table = curve_table("small", X, y, SMALL_OPTS, SMALL_SETTINGS, results)
    assert list(table.columns) == CURVE_COLUMNS
    assert len(table) == 6
    assert sorted(os.listdir(results)) == ["gcv_small_grid3x3.csv",
                                           "power_small_runs2_perms9_grid3x3.csv"]
    power = power_curve(X, y, 3, 3, 0.5, 9, 2, 0.05, 8891, 1000)
    close(table["power"], power["power"])
    close(table["gcv"], gcv_curve(X, y, 3, 3, 0.5, 8891, 1000)["gcv"])
    # the second call reads the stored parts
    capsys.readouterr()
    again = curve_table("small", X, y, SMALL_OPTS, SMALL_SETTINGS, results)
    assert "reusing" in capsys.readouterr().out
    close(again.to_numpy(), table.to_numpy(), rtol=1e-13)
