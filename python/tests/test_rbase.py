# -----------------------------------------------------------------------------
# File:    python/tests/test_rbase.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Unit tests of python/adamant/rbase.py: what base R, stats and MASS supply, on
#          inputs with a known answer: matrix products, svd, eigen and ginv
#          against numpy, sums, means and variances, scale, quantile (type
#          7), seq, which.min, pgamma and mvrnorm.
# Inputs:  None.
# Outputs: None.
# Usage:   python -m pytest python/tests/test_rbase.py     # inside `nix develop`
# -----------------------------------------------------------------------------
import math

import numpy as np
import pytest

from adamant import rbase, rng
from inputs import close, seeded_matrix


def test_as_matrix_is_a_float_matrix_in_column_major_storage():
    column = rbase.as_matrix([1, 2, 3])
    assert column.shape == (3, 1) and column.dtype == np.float64
    matrix = rbase.as_matrix(np.arange(6).reshape(2, 3))
    assert matrix.flags["F_CONTIGUOUS"] and matrix.tolist() == [[0, 1, 2], [3, 4, 5]]


def test_matprod_and_crossprod_are_the_matrix_products():
    A = seeded_matrix(4, 3, seed=51)
    B = seeded_matrix(3, 2, seed=52)
    v = seeded_matrix(3, 1, seed=53)[:, 0]
    close(rbase.matprod(A, B), A @ B)
    close(rbase.matprod(A, v), (A @ v)[:, None])            # a vector is one column
    close(rbase.matprod(A[0:1, :], B), A[0:1, :] @ B)
    close(rbase.crossprod(A, A), A.T @ A)
    close(rbase.crossprod(A, A[:, 0]), (A.T @ A[:, 0])[:, None])
    close(rbase.crossprod(A[:, 0], A), (A[:, 0] @ A)[None, :])
    close(rbase.crossprod_columns(A, A[:, :2]), A.T @ A[:, :2])
    assert rbase.matprod(A, B).flags["F_CONTIGUOUS"]


def test_matprod_carries_missing_values_as_r_does():
    A = np.array([[1.0, np.nan], [2.0, 3.0]])
    product = rbase.matprod(A, np.array([[1.0], [1.0]]))
    assert np.isnan(product[0, 0]) and product[1, 0] == 5
    assert np.isnan(rbase.crossprod(A, np.array([[1.0], [1.0]]))[1, 0])


def test_svd_returns_the_thin_decomposition():
    for n, p in [(6, 4), (4, 6)]:
        X = seeded_matrix(n, p, seed=54)
        X_svd = rbase.svd(X)
        k = min(n, p)
        assert X_svd["u"].shape == (n, k) and X_svd["v"].shape == (p, k)
        close(X_svd["d"], np.linalg.svd(X, compute_uv=False))
        close(X_svd["u"] * X_svd["d"] @ X_svd["v"].T, X)
        close(X_svd["u"].T @ X_svd["u"], np.eye(k))
    close(rbase.svd(np.diag([1.0, 3.0, 2.0]))["d"], [3, 2, 1])


def test_eigen_symmetric_returns_decreasing_values_and_their_vectors():
    S = np.array([[2.0, 1.0], [1.0, 2.0]])
    eig = rbase.eigen_symmetric(S)
    close(eig["values"], [3, 1])
    close(np.abs(eig["vectors"]), np.full((2, 2), math.sqrt(0.5)))
    close(eig["vectors"] * eig["values"] @ eig["vectors"].T, S)


def test_ginv_is_the_inverse_or_the_generalized_inverse():
    X = seeded_matrix(4, 4, seed=55)
    close(rbase.ginv(X), np.linalg.inv(X))
    singular = np.array([[1.0, 2.0], [2.0, 4.0]])     # rank one
    close(rbase.ginv(singular), singular / 25)
    close(rbase.ginv(np.zeros((2, 3))), np.zeros((3, 2)))
    tall = seeded_matrix(5, 2, seed=56)
    close(rbase.ginv(tall), np.linalg.pinv(tall))


def test_sums_means_and_variances_are_those_of_r():
    A = np.array([[1.0, 4.0], [2.0, 5.0], [3.0, 7.0]])
    close(rbase.col_sums(A), [6, 16])
    assert rbase.r_sum([0.5, 0.25, 0.25]) == 1
    assert rbase.r_mean([1, 2, 3, 6]) == 3
    assert rbase.r_var([1, 2, 3, 6]) == 14 / 3
    close(rbase.r_sd([1, 2, 3, 6]), math.sqrt(14 / 3))
    assert math.isnan(rbase.r_mean([1, np.nan]))
    assert rbase.r_mean([1, np.nan, 3], na_rm=True) == 2
    assert rbase.r_var([1, np.nan, 3], na_rm=True) == 2
    assert math.isnan(rbase.r_var([1])) and math.isnan(rbase.r_mean([]))


def test_scale_centres_and_divides_by_the_sd():
    close(rbase.scale([1, 2, 3]), [[-1], [0], [1]])
    X = seeded_matrix(8, 3, seed=57)
    scaled = rbase.scale(X)
    close(scaled, (X - X.mean(axis=0)) / X.std(axis=0, ddof=1))
    close(scaled.mean(axis=0), np.zeros(3))


@pytest.mark.parametrize("prob, expected", [(0, 1), (0.25, 1.75), (0.5, 2.5), (0.9, 3.7),
                                            (1, 4)])
def test_quantile_is_type_7(prob, expected):
    close(rbase.quantile([4, 1, 3, 2], prob), expected)


def test_seq_length_and_seq_by_are_the_sequences_of_r():
    close(rbase.seq_length(0, 1, 5), [0, 0.25, 0.5, 0.75, 1])
    close(rbase.seq_length(2, 7, 2), [2, 7])
    close(rbase.seq_length(2, 7, 1), [2])
    close(rbase.seq_length(3, 3, 4), [3, 3, 3, 3])
    assert rbase.seq_length(-3, 1000, 7)[-1] == 1000
    close(rbase.seq_by(4, -2.5, -0.5), np.arange(14) * -0.5 + 4)
    close(rbase.seq_by(0, 1, 0.3), [0, 0.3, 0.6, 0.9])


def test_which_min_rows_is_the_first_smallest_value_of_every_row():
    A = np.array([[3.0, 1.0, 1.0], [np.nan, 5.0, 4.0], [2.0, 2.0, 9.0]])
    assert rbase.which_min_rows(A).tolist() == [1, 2, 0]


def test_pgamma_upper_is_the_upper_tail_of_the_gamma():
    close(rbase.pgamma_upper(3.0, 4.0, 2.0), 61 * math.exp(-6))   # P(Poisson(6) <= 3)
    close(rbase.pgamma_upper(np.array([0.0, 1.0, 2.0]), 1.0, 1.0), np.exp([-0.0, -1.0, -2.0]))


def test_mvrnorm_draws_rows_with_the_given_covariance():
    Sigma = np.array([[1.0, 0.5], [0.5, 2.0]])
    rng.set_seed(58)
    X = rbase.mvrnorm(20000, Sigma)
    assert X.shape == (20000, 2)
    close(np.cov(X.T), Sigma, atol=0.05, rtol=0)
    close(X.mean(axis=0), [0, 0], atol=0.03)
    rng.set_seed(58)
    assert np.array_equal(rbase.mvrnorm(20000, Sigma), X)   # the draw of its seed
    # a diagonal covariance scales independent normal numbers
    rng.set_seed(59)
    D = rbase.mvrnorm(5, np.diag([4.0, 1.0]))
    rng.set_seed(59)
    Z = rng.rnorm(10).reshape((5, 2), order="F")
    close(np.abs(D), np.abs(Z * [2, 1]))
    with pytest.raises(ValueError, match="not positive definite"):
        rbase.mvrnorm(3, np.array([[1.0, 2.0], [2.0, 1.0]]))
