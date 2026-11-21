# -----------------------------------------------------------------------------
# File:    python/tests/test_mantel.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Unit tests of python/adamant/mantel.py: every function on inputs with a known
#          answer (hand-computed values, and the closed form of the ridge
#          statistic); the cases of tests/testthat/test-mantel.R.
# Inputs:  None.
# Outputs: None.
# Usage:   python -m pytest python/tests/test_mantel.py     # inside `nix develop`
# -----------------------------------------------------------------------------
import math

import numpy as np

from adamant import mantel
from inputs import close, seeded_matrix

EPS = math.sqrt(np.finfo(float).eps)


def test_percent_rank_is_the_fraction_strictly_below():
    close(mantel.percent_rank([3, 1, 2]), [1, 0, 0.5])
    close(mantel.percent_rank([1, 1, 2]), [0, 0, 1])
    close(mantel.percent_rank([5, 1, 3, 2]), np.array([3, 0, 2, 1]) / 3)


def test_truncated_svd_zeroes_the_singular_values_below_the_quantile():
    X = np.diag([4.0, 3.0, 2.0, 1.0])
    close(mantel.truncated_svd(X, 0.5)["d"], [4, 3, 2, 1])
    close(mantel.truncated_svd(X, 0.5)["d_tilde"], [4, 3, 0, 0])
    close(mantel.truncated_svd(X, 0)["d_tilde"], [4, 3, 2, 1])
    close(mantel.truncated_svd(X, 0.25)["d_tilde"], [4, 3, 2, 0])

    X = seeded_matrix(6, 4, seed=1)
    X_svd = mantel.truncated_svd(X, 0)
    assert X_svd["u"].shape == (6, 4)
    close(X_svd["u"].T @ X_svd["u"], np.eye(4))
    close(X_svd["d"] ** 2, np.linalg.eigvalsh(X.T @ X)[::-1])


def test_lambda_grid_has_the_negative_part_zero_and_the_positive_part():
    # sqrt(p / n) = 2, smallest kept singular value 3: lower bound -(2 - 1)^2 * 3
    grid = mantel.lambda_grid(np.array([4.0, 3.0, 0.0, 0.0]), 4, 16, 4, 3)
    close(grid, [-3, (-3 - EPS) / 2, -EPS, 0, EPS, (EPS + 1000) / 2, 1000], atol=0)
    assert mantel.lambda_grid(np.array([4.0, 3.0]), 4, 16, 4, 3, upper=10).max() == 10
    # n = p: the lower bound is zero and the grid has no negative part
    close(mantel.lambda_grid(np.array([2.0, 1.0]), 5, 5, 4, 3), [0, 0, 500, 1000])


def test_ridge_weights_are_the_ratios_and_one_at_zero():
    weights = mantel.ridge_weights(np.array([2.0, 1.0, 0.0]), np.array([-0.5, 0.0, 1.0]))
    close(weights, np.column_stack([[4 / 3.5, 2, 0], [1, 1, 1], [0.8, 0.5, 0]]))
    # a penalty below -d^2 gives a negative ratio, which is set to zero
    close(mantel.ridge_weights(np.array([2.0, 1.0]), np.array([-2.0]))[:, 0], [2, 0])


def test_mantel_statistics_are_the_weighted_sums_of_squares_in_both_forms():
    Z = np.array([[1.0, 2.0], [3.0, 4.0]])            # two data sets, one per column
    weights = np.array([[1.0, 0.5], [1.0, 0.25]])     # two penalties, one per column
    expected = [[10, 2.75], [20, 6]]                  # 1 + 9, 0.5 + 2.25; 4 + 16, 2 + 4
    lambdas = np.array([0.0, 1.0])
    close(mantel.mantel_statistics(Z, weights, lambdas), expected)
    close(mantel.mantel_statistics(Z, weights, lambdas, exact=True), expected)


def test_the_statistic_is_the_closed_form_of_ridge_regression():
    X = seeded_matrix(8, 5, seed=2)
    y = seeded_matrix(8, 1, seed=3)
    lambdas = np.array([0.0, 0.5, 3.0, 1000.0])
    X_svd = mantel.truncated_svd(X, 0)
    statistics = mantel.mantel_statistics(X_svd["u"].T @ y,
                                          mantel.ridge_weights(X_svd["d_tilde"], lambdas),
                                          lambdas, exact=True)
    closed_form = [(y.T @ X @ np.linalg.solve(X.T @ X + lam * np.eye(5), X.T @ y)).item()
                   for lam in lambdas]
    close(statistics[0, :], closed_form)


def test_permutation_pvalues_are_exact_for_a_small_case():
    statistics = np.column_stack([[5, 1, 3, 2], [1, 2, 3, 4], [2, 2, 2, 2]]).astype(float)
    expected = np.column_stack([[0, 1, 1 / 3, 2 / 3], [1, 2 / 3, 1 / 3, 0], [1, 1, 1, 1]])
    close(mantel.permutation_pvalues(statistics), expected)
    close(mantel.observed_pvalues(statistics), expected[0, :])

    random = seeded_matrix(20, 6, seed=4) ** 2
    p = mantel.permutation_pvalues(random)
    assert ((p >= 0) & (p <= 1)).all()
    close(mantel.observed_pvalues(random), p[0, :])
    # without ties every multiple of 1 / 19 appears once per penalty
    close(np.sort(p[:, 0]), np.arange(20) / 19)


def test_gamma_pvalues_match_the_gamma_with_the_moments_of_the_permutations():
    # permuted statistics 1, 2, 3: mean 2, variance 1, so shape 4 and rate 2;
    # P(Gamma(4, 2) > 3) = P(Poisson(6) <= 3) = 61 exp(-6)
    statistics = np.column_stack([[3, 1, 2, 3], [6, 2, 4, 6]]).astype(float)
    p = mantel.gamma_pvalues(statistics)
    assert p.shape == (4, 2)
    close(p[0, 0], 61 * math.exp(-6))
    close(p[0, 1], 61 * math.exp(-6))          # the same after doubling every statistic
    close(p[1, 0], 19 / 3 * math.exp(-2))      # P(Poisson(2) <= 3)
    close(mantel.gamma_pvalues(statistics, rows=[0]), p[0:1, :])
    assert ((p >= 0) & (p <= 1)).all()
