# -----------------------------------------------------------------------------
# File:    python/adamant/liver_cv.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Figure 1.
# Purpose: K-fold cross-validated MSE of ridge regression (as glmnet computes
#          it) over a grid of penalties, for every response of a matrix, and
#          the plot of the curves.
# Inputs:  None.
# Outputs: None.
# Usage:   from adamant.liver_cv import cv_curves
#          cv_curves(X, Y, lambdas, 10)
# Source:  R/liver_cv.R.
# -----------------------------------------------------------------------------
import math

import numpy as np
import pandas as pd

from . import plots, rng
from .glmnet import cv_glmnet_ridge


def fold_assignment(n, n_folds):
    """The folds `glmnet::cv.glmnet()` draws: 1, ..., n_folds repeated to
    length n, in random order.

    Parameters: n (number of observations), n_folds (number of folds).
    Returns: int array of length n.
    """
    return rng.sample(np.resize(np.arange(1, n_folds + 1), n))


def skip_fold_assignments(count, n, n_folds):
    """Draw and discard fold assignments.

    Advances the random number stream by what `count` calls of
    `glmnet::cv.glmnet()` draw (one fold assignment each).

    Parameters: count (number of fold assignments), n (number of
    observations), n_folds (number of folds). Returns: None.
    """
    for _ in range(count):
        fold_assignment(n, n_folds)


def cv_curves(X, Y, lambdas, n_folds):
    """Cross-validated MSE of ridge regression for every response.

    The fits draw their folds from the current random number stream, one
    response after the other.

    Parameters: X (predictors), Y (responses, one column each), lambdas
    (decreasing penalties), n_folds (number of folds).
    Returns: long DataFrame with columns variable ("Variable<j>"),
    log_lambda, cvm (cross-validated MSE), at_minimum (True where cvm is
    smallest for the variable) and log_lambda_min.
    """
    log_lambda = np.array([math.log(lam) for lam in lambdas])
    curves = []
    for j in range(Y.shape[1]):
        fit = cv_glmnet_ridge(X, Y[:, j], lambdas, fold_assignment(X.shape[0], n_folds))
        curves.append(pd.DataFrame({"variable": f"Variable{j + 1}", "log_lambda": log_lambda,
                                    "cvm": fit["cvm"],
                                    "at_minimum": fit["cvm"] == fit["cvm"].min(),
                                    "log_lambda_min": math.log(fit["lambda_min"])}))
    return pd.concat(curves, ignore_index=True)


def cv_plot(curves, y_max):
    """Plot of the cross-validated MSE against log(lambda), one line per
    variable and a point at each minimum.

    Parameters: curves (output of `cv_curves()`), y_max (upper limit of the
    MSE axis).
    Returns: a figure (`plots.Plot`).
    """
    variables = list(dict.fromkeys(curves["variable"]))
    colours = plots.hue_colours(len(variables))
    plot = plots.Plot(theme="bw")
    shown = curves.copy()
    shown.loc[shown["cvm"] > y_max, "cvm"] = np.nan  # outside the axis: not drawn
    for group, variable in enumerate(variables, start=1):
        line = shown[shown["variable"] == variable].sort_values("log_lambda")
        plot.ax.plot(line["log_lambda"], line["cvm"], color=colours[group - 1], linewidth=2.0)
        plot.record(1, x=line["log_lambda"], y=line["cvm"], group=group,
                    colour=colours[group - 1])
    points = shown[shown["at_minimum"]]
    plot.ax.scatter(points["log_lambda"], points["cvm"], color="#0000FF", s=12, zorder=3)
    plot.record(2, x=points["log_lambda"], y=points["cvm"], group=-1, colour="blue")
    plot.ax.set_ylim(top=y_max)
    plot.ax.set_xlabel("log(lambda)")
    plot.ax.set_ylabel("Cross-validated MSE")
    return plot


def cv_minima(curves, panel):
    """The minimum of every curve.

    Parameters: curves (output of `cv_curves()`), panel (label of the
    panel).
    Returns: DataFrame with columns panel, variable, log_lambda_min, min_cvm.
    """
    lowest = curves[curves["at_minimum"] & (curves["log_lambda"] == curves["log_lambda_min"])]
    return pd.DataFrame({"panel": panel, "variable": lowest["variable"].to_numpy(),
                         "log_lambda_min": lowest["log_lambda_min"].to_numpy(),
                         "min_cvm": lowest["cvm"].to_numpy()})
