# -----------------------------------------------------------------------------
# File:    python/tests/test_liver_cv.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Unit tests of python/adamant/liver_cv.py: the fold assignments drawn and skipped,
#          the cross-validated curves against a ridge fit in closed form and
#          against stored values, their plot and their minima; the cases of
#          tests/testthat/test-liver_cv.R.
# Inputs:  None.
# Outputs: None.
# Usage:   python -m pytest python/tests/test_liver_cv.py     # inside `nix develop`
# -----------------------------------------------------------------------------
import numpy as np
import pandas as pd

from adamant import plots, rng
from adamant.glmnet import cv_glmnet_ridge
from adamant.liver_cv import cv_curves, cv_minima, cv_plot, fold_assignment, skip_fold_assignments
from inputs import close, seeded_matrix

LAMBDAS = np.array([10, 1, 0.1, 0.01])
# cvm of cv_curves() on the small problem with seed 3, as first computed with
# R 4.6.1 and glmnet 5.0; tests/testthat/test-liver_cv.R has the same.
STORED_CVM = [1.9565450234509685, 0.74109508987441752, 0.32795538515868683,
              0.32300771374188475, 1.1202080079946386, 0.61124977338760755,
              0.41355989194038351, 0.41167216021214931]


def small_problem():
    """The small test problem: 30 observations, 4 predictors, 2 responses.
    Returns: (X, Y)."""
    X = seeded_matrix(30, 4, seed=21)
    noise = seeded_matrix(30, 2, seed=22)
    return X, np.column_stack([X[:, 0] - X[:, 1], X[:, 2]]) + 0.5 * noise


def small_curves(seed=3):
    """cv_curves() of the small problem with the folds of a seed."""
    X, Y = small_problem()
    rng.set_seed(seed)
    return cv_curves(X, Y, LAMBDAS, 5)


def test_fold_assignment_holds_every_fold_equally_often():
    rng.set_seed(5)
    folds = fold_assignment(10, 5)
    assert sorted(folds) == [1, 1, 2, 2, 3, 3, 4, 4, 5, 5]
    assert np.bincount(fold_assignment(64, 10))[1:].tolist() == [7, 7, 7, 7, 6, 6, 6, 6, 6, 6]


def test_skip_fold_assignments_advances_the_stream_by_one_fold_draw_per_count():
    rng.set_seed(5)
    assert skip_fold_assignments(2, 10, 5) is None
    after_skipping = rng.runif(1)[0]
    rng.set_seed(5)
    first = fold_assignment(10, 5)
    second = fold_assignment(10, 5)
    assert after_skipping == rng.runif(1)[0]
    assert not np.array_equal(first, second)
    rng.set_seed(5)
    skip_fold_assignments(0, 10, 5)
    assert np.array_equal(fold_assignment(10, 5), first)


def test_cv_curves_is_a_long_table_with_one_curve_per_response():
    curves = small_curves()
    assert list(curves.columns) == ["variable", "log_lambda", "cvm", "at_minimum",
                                    "log_lambda_min"]
    assert list(curves["variable"]) == ["Variable1"] * 4 + ["Variable2"] * 4
    close(curves["log_lambda"], np.tile(np.log(LAMBDAS), 2))
    assert (curves["cvm"] > 0).all()
    for _, curve in curves.groupby("variable"):
        assert curve["at_minimum"].sum() == 1
        assert curve["cvm"][curve["at_minimum"]].item() == curve["cvm"].min()
        assert curve["log_lambda_min"].unique().tolist() == \
            curve["log_lambda"][curve["at_minimum"]].tolist()


def test_cv_curves_draws_the_folds_of_each_response_from_the_stream_in_turn():
    X, Y = small_problem()
    curves = small_curves()
    rng.set_seed(3)
    folds = [fold_assignment(30, 5) for _ in range(2)]
    for j in range(2):
        fit = cv_glmnet_ridge(X, Y[:, j], LAMBDAS, folds[j])
        assert np.array_equal(curves["cvm"][curves["variable"] == f"Variable{j + 1}"],
                              fit["cvm"])
    assert not np.allclose(small_curves(seed=4)["cvm"], curves["cvm"])


def test_the_cross_validated_error_is_that_of_ridge_fits_in_closed_form():
    # glmnet minimizes sum((y - a - X b)^2) / (2 n) + lambda * |b|^2 / 2 on
    # predictors and response scaled to variance one (divisor n): on centred
    # data b = (X'X + n lambda sd(y) I)^-1 X'y in the units of the response
    X, Y = small_problem()
    y = Y[:, 0]
    rng.set_seed(3)
    folds = fold_assignment(30, 5)
    errors = []
    for fold in range(1, 6):
        train = folds != fold
        centre, spread = X[train].mean(axis=0), X[train].std(axis=0)
        Z = (X[train] - centre) / spread
        b = np.linalg.solve(Z.T @ Z + train.sum() * LAMBDAS[1] / y[train].std() * np.eye(4),
                            Z.T @ (y[train] - y[train].mean()))
        predicted = y[train].mean() + (X[~train] - centre) / spread @ b
        errors.append(np.mean((y[~train] - predicted) ** 2))
    close(small_curves()["cvm"][1], np.mean(errors), rtol=1e-4)


def test_cv_curves_gives_the_stored_values_of_the_small_problem():
    close(small_curves()["cvm"], STORED_CVM, rtol=1e-10)


def test_cv_plot_draws_one_line_per_variable_and_a_point_at_each_minimum():
    curves = small_curves()
    plot = cv_plot(curves, y_max=10)
    assert isinstance(plot, plots.Plot)
    marks = plot.plotted_data()
    lines, points = marks[marks["layer"] == 1], marks[marks["layer"] == 2]
    assert len(lines) == 8
    assert lines["group"].nunique() == 2
    close(np.sort(points["y"]), np.sort(curves["cvm"][curves["at_minimum"]]))
    assert plot.ax.get_xlabel() == "log(lambda)"


def test_cv_minima_has_the_lowest_point_of_every_curve():
    curves = pd.DataFrame({"variable": ["Variable1"] * 3 + ["Variable2"] * 3,
                           "log_lambda": [2.0, 1.0, 0.0] * 2, "cvm": [5.0, 3, 4, 1, 2, 6],
                           "at_minimum": [False, True, False, True, False, False],
                           "log_lambda_min": [1.0] * 3 + [2.0] * 3})
    pd.testing.assert_frame_equal(
        cv_minima(curves, "A"),
        pd.DataFrame({"panel": "A", "variable": ["Variable1", "Variable2"],
                      "log_lambda_min": [1.0, 2.0], "min_cvm": [3.0, 1.0]}))
