# -----------------------------------------------------------------------------
# File:    python/tests/test_adamant_cv.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Unit tests of python/adamant/adamant_cv.py: the percent rank of the original runs
#          on hand-computed cases, and the test on a small seeded design
#          (structure, the statistic of the observed data in closed form, the
#          permutations of the seed, stored values); the cases of
#          tests/testthat/test-adamant_cv.R.
# Inputs:  None.
# Outputs: None.
# Usage:   python -m pytest python/tests/test_adamant_cv.py     # inside `nix develop`
# -----------------------------------------------------------------------------
import numpy as np

from adamant import rng
from adamant.adamant_cv import TESTS, adamant_cv, adaptive_pvalue, percent_rank_cumulative
from adamant.mantel import percent_rank
from inputs import close, seeded_matrix

# Values of adamant_cv(X, y, 0.5, 4, 5, 19, seed = 7) on the small problem, as
# first computed with R 4.6.1; tests/testthat/test-adamant_cv.R has the same.
STORED = {"lambda_first": -1.0926875988102316, "statistic_last": 0.32052880115284371,
          "gcv_last": 1.1434087615079229, "p_gamma_zero": 4.8335366933517551e-05}


def small_problem():
    """The small test problem: 20 observations, 5 predictors, signal in the
    first. Returns: (X, y)."""
    X = seeded_matrix(20, 5, seed=11)
    return X, X[:, 0] + 0.5 * seeded_matrix(20, 1, seed=12)[:, 0]


def ridge_statistic(X, y, lam):
    """y'X (X'X + lambda I)^-1 X'y."""
    return y @ X @ np.linalg.solve(X.T @ X + lam * np.eye(X.shape[1]), X.T @ y)


def test_percent_rank_cumulative_adds_the_counts_of_the_distinct_values():
    # distinct values 0.1, 0.2, 0.4 with counts 1, 2, 1 of m - 1 = 3
    close(percent_rank_cumulative([0.2, 0.1, 0.2, 0.4]), [1 / 3, 0, 1 / 3, 1])
    # without ties it is the usual percent rank
    close(percent_rank_cumulative([3, 1, 2]), percent_rank([3, 1, 2]))


def test_adaptive_pvalue_is_the_rank_of_the_first_element():
    close(adaptive_pvalue([0.2, 0.1, 0.2, 0.4]), 1 / 3)
    assert adaptive_pvalue([0.01, 0.5, 0.3]) == 0
    assert adaptive_pvalue([0.9, 0.5, 0.3]) == 1


def test_adamant_cv_returns_every_part_with_the_right_shape():
    X, y = small_problem()
    ada = adamant_cv(X, y, 0, 4, 5, 19, seed=7)
    assert list(ada) == ["p_values", "lambdas", "lambda_gcv", "lambda_loo", "gcv", "loo",
                         "statistics", "p_perm", "p_gamma"]
    assert list(ada["p_values"]) == TESTS
    assert all(0 <= p <= 1 for p in ada["p_values"].values())
    assert len(ada["lambdas"]) == 9
    assert (ada["lambdas"] == 0).sum() == 1
    # one row for the observed data and each of the 19 permutations
    for part in ("gcv", "loo", "statistics", "p_perm", "p_gamma"):
        assert ada[part].shape == (1 + 19, 9)
    assert len(ada["lambda_gcv"]) == 1 + 19
    assert np.isin(ada["lambda_gcv"], ada["lambdas"]).all()
    assert np.isin(ada["lambda_loo"], ada["lambdas"]).all()
    assert ((ada["p_perm"] >= 0) & (ada["p_perm"] <= 1)).all()
    assert ((ada["p_gamma"] >= 0) & (ada["p_gamma"] <= 1)).all()


def test_row_1_is_the_observed_data_and_row_2_the_first_permutation_of_the_seed():
    X, y = small_problem()
    ada = adamant_cv(X, y, 0, 4, 5, 19, seed=7)
    positive = ada["lambdas"] > 0
    close(ada["statistics"][0, positive],
          [ridge_statistic(X, y, lam) for lam in ada["lambdas"][positive]])
    rng.set_seed(7)
    permuted = y[rng.sample_int(20) - 1]
    close(ada["statistics"][1, positive],
          [ridge_statistic(X, permuted, lam) for lam in ada["lambdas"][positive]])


def test_adamant_cv_depends_on_the_seed_only_and_finds_the_signal():
    X, y = small_problem()
    ada = adamant_cv(X, y, 0, 4, 5, 19, seed=7)
    rng.runif(3)  # the state of the random stream before the call does not matter
    again = adamant_cv(X, y, 0, 4, 5, 19, seed=7)
    assert again["p_values"] == ada["p_values"]
    for part in ("statistics", "p_perm", "p_gamma", "gcv", "loo", "lambda_gcv"):
        assert np.array_equal(again[part], ada[part])
    close(ada["statistics"], adamant_cv(X, y, 0, 4, 5, 19, seed=7, exact=False)["statistics"])
    other = adamant_cv(X, y, 0, 4, 5, 19, seed=8)
    close(other["statistics"][0, :], ada["statistics"][0, :])
    assert not np.allclose(other["statistics"][1, :], ada["statistics"][1, :])
    # the observed statistic is above all 19 permuted ones at lambda = 0
    assert ada["p_perm"][0, ada["lambdas"] == 0] == 0
    assert ada["p_values"]["gcv"] == 0


def test_observed_lambda0_replaces_the_p_values_of_the_observed_data_at_zero():
    X, y = small_problem()
    ada = adamant_cv(X, y, 0, 4, 5, 19, seed=7,
                     observed_lambda0={"perms_below": 15, "gamma": 0.123})
    zero = ada["lambdas"] == 0
    close(ada["p_perm"][0, zero], 1 - 15 / 19)
    close(ada["p_gamma"][0, zero], 0.123)


def test_adamant_cv_gives_the_stored_values_of_the_small_problem():
    X, y = small_problem()
    ada = adamant_cv(X, y, 0.5, 4, 5, 19, seed=7)
    close(ada["lambdas"][0], STORED["lambda_first"], rtol=1e-10)
    close(ada["statistics"][0, 8], STORED["statistic_last"], rtol=1e-10)
    close(ada["gcv"][0, 8], STORED["gcv_last"], rtol=1e-10)
    close(ada["p_gamma"][0, 4], STORED["p_gamma_zero"], rtol=1e-8, atol=0)
