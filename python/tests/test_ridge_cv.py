# -----------------------------------------------------------------------------
# File:    python/tests/test_ridge_cv.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Unit tests of python/adamant/ridge_cv.py: the ridge operators and the
#          cross-validation scores against direct solves and against
#          leave-one-out refits; the cases of tests/testthat/test-ridge_cv.R.
# Inputs:  None.
# Outputs: None.
# Usage:   python -m pytest python/tests/test_ridge_cv.py     # inside `nix develop`
# -----------------------------------------------------------------------------
import numpy as np
import pytest

from adamant.ridge_cv import matrix_trace, penalty_scores, ridge_cv, ridge_operator
from inputs import close, seeded_matrix


def refit_loo(X, y, penalty):
    """Leave-one-out mean squared error of ridge regression by refitting.

    Parameters: X, y (predictors and response), penalty (number added to the
    diagonal of X'X).
    Returns: mean over i of the squared error of observation i predicted
    from the fit without it.
    """
    errors = []
    for i in range(X.shape[0]):
        keep = np.arange(X.shape[0]) != i
        beta = np.linalg.solve(X[keep].T @ X[keep] + penalty * np.eye(X.shape[1]),
                               X[keep].T @ y[keep])
        errors.append(y[i] - X[i] @ beta)
    return np.mean(np.square(errors))


def gram_of(X):
    """X'X as a Fortran-ordered matrix."""
    return np.asfortranarray(X.T @ X)


def test_matrix_trace_adds_the_diagonal():
    assert matrix_trace(np.diag([1.0, 2.0, 3.0])) == 6
    assert matrix_trace(np.array([[1.0, 3.0], [2.0, 4.0]])) == 5


def test_ridge_operator_when_n_is_above_p():
    X = seeded_matrix(8, 3, seed=5)
    op = ridge_operator(X, np.asfortranarray(X.T), gram_of(X), 0.5, zero_rule=True)
    inverse = np.linalg.inv(X.T @ X + 8 * 0.5 * np.eye(3))
    close(op["inverse"], inverse)
    close(op["W"], inverse @ X.T)
    close(op["L"], X @ inverse @ X.T)


def test_ridge_operator_when_n_is_at_most_p():
    X = seeded_matrix(4, 6, seed=6)
    gram = np.asfortranarray(X @ X.T)
    op = ridge_operator(X, np.asfortranarray(X.T), gram, 0.5, zero_rule=True)
    close(op["inverse"], np.linalg.inv(gram + 4 * 0.5 * np.eye(4)))
    close(op["W"], X.T @ op["inverse"])
    close(op["L"], gram @ op["inverse"])

    zero = ridge_operator(X, np.asfortranarray(X.T), gram, 0, zero_rule=True)
    close(zero["L"], np.eye(4) - np.linalg.inv(gram))
    plain = ridge_operator(X, np.asfortranarray(X.T), gram, 0, zero_rule=False)
    close(plain["L"], np.eye(4), atol=1e-10)   # X X' (X X')^-1


def test_penalty_scores_are_the_gcv_formula_and_the_leave_one_out_refits():
    X = seeded_matrix(10, 3, seed=7)
    y = seeded_matrix(10, 1, seed=8)
    lam = 0.3
    scores = penalty_scores(X, np.asfortranarray(X.T), gram_of(X), lam, y, feb2022=True)
    L = X @ np.linalg.inv(X.T @ X + 10 * lam * np.eye(3)) @ X.T
    residual = y - L @ y
    close(scores["gcv"], np.mean((residual / (1 - np.trace(L) / 10)) ** 2))
    close(scores["loo"], refit_loo(X, y[:, 0], 10 * lam))


def test_the_residual_of_a_later_column_is_the_observed_response_minus_its_fit():
    X = seeded_matrix(10, 3, seed=7)
    Y = seeded_matrix(10, 2, seed=9)
    scores = penalty_scores(X, np.asfortranarray(X.T), gram_of(X), 0.3, Y, feb2022=True)
    L = X @ np.linalg.inv(X.T @ X + 3 * np.eye(3)) @ X.T
    residual = Y[:, 0] - L @ Y[:, 1]
    close(scores["loo"][1], np.mean((residual / (1 - np.diag(L))) ** 2))
    assert len(scores["gcv"]) == 2


def test_ridge_cv_has_one_row_per_response_and_one_column_per_penalty():
    X = seeded_matrix(10, 3, seed=7)
    Y = seeded_matrix(10, 2, seed=9)
    lambdas = [0.1, 0.3, 2]
    cv = ridge_cv(X, Y, lambdas)
    assert cv["gcv"].shape == (2, 3) and cv["loo"].shape == (2, 3)
    close(cv["loo"][0, :], [refit_loo(X, Y[:, 0], 10 * lam) for lam in lambdas])
    one = penalty_scores(X, np.asfortranarray(X.T), gram_of(X), 0.3, Y, feb2022=True)
    close(cv["gcv"][:, 1], one["gcv"])
    # the two versions of the original code agree away from lambda = 0
    dec = ridge_cv(X, Y, lambdas, variant="dec2021")
    close(dec["gcv"], cv["gcv"])
    close(dec["loo"], cv["loo"])
    with pytest.raises(ValueError):
        ridge_cv(X, Y, lambdas, variant="other")


def test_ridge_cv_does_not_depend_on_cores():
    X = seeded_matrix(10, 3, seed=7)
    Y = seeded_matrix(10, 2, seed=9)
    one = ridge_cv(X, Y, [0.1, 0.3, 2], cores=1)
    two = ridge_cv(X, Y, [0.1, 0.3, 2], cores=2)
    assert np.array_equal(one["gcv"], two["gcv"]) and np.array_equal(one["loo"], two["loo"])


def test_at_lambda_zero_with_n_at_most_p_the_feb_2022_scores_are_finite():
    X = seeded_matrix(4, 6, seed=6)
    y = seeded_matrix(4, 1, seed=10)
    # the smoother I - (XX')^-1 keeps the denominators away from zero, and the
    # residual of an exact fit is zero up to rounding
    feb = ridge_cv(X, y, [0, 1], variant="feb2022")
    assert np.isfinite(feb["gcv"]).all() and np.isfinite(feb["loo"]).all()
    assert feb["gcv"][0, 0] < 1e-20 and feb["loo"][0, 0] < 1e-20
    assert feb["gcv"][0, 1] > 1e-3
    dec = ridge_cv(X, y, [0, 1], variant="dec2021")
    close(dec["gcv"][0, 1], feb["gcv"][0, 1])
