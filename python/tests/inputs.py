# -----------------------------------------------------------------------------
# File:    python/tests/inputs.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Inputs shared by the unit tests: small seeded matrices (the same
#          numbers as R draws), small options and settings, made-up data
#          files. tests/testthat/helper-inputs.R has the same.
# Inputs:  None.
# Outputs: None.
# Usage:   from inputs import seeded_matrix     # in a test under python/tests/
# -----------------------------------------------------------------------------
import os

import numpy as np
import pandas as pd

from adamant import rng
from adamant.setup import REPO_ROOT, write_csv

REPO = REPO_ROOT

# Settings as PAPER of adamant/setup.py with few permutations.
SMALL_SETTINGS = {"seed": 8891, "alpha": 0.05, "sigma": 0.5, "upper": 1000, "n_perms": 9,
                  "figure_size": 3}
# Options of a small rerun: 2 replications, 9 permutations, a 3 + 3 grid.
SMALL_OPTS = {"rerun": True, "runs": 2, "perms": 9, "negative": 3, "positive": 3, "cores": 1}


def seeded_matrix(n, p, seed):
    """A matrix of standard normal numbers drawn with its own seed, as
    `set.seed(seed); matrix(rnorm(n * p), n, p)` in R.

    Parameters: n, p (numbers of rows and columns), seed (seed of the draw).
    Returns: n x p Fortran-ordered array.
    """
    rng.set_seed(seed)
    return np.asfortranarray(rng.rnorm(n * p).reshape((n, p), order="F"))


def close(actual, expected, rtol=1.5e-8, atol=1e-12):
    """Assert that two arrays are equal up to rounding (the tolerance of
    testthat's `expect_equal()`).

    Parameters: actual, expected (array-likes), rtol, atol (tolerances).
    Returns: None.
    """
    np.testing.assert_allclose(np.asarray(actual, dtype=np.float64),
                               np.asarray(expected, dtype=np.float64), rtol=rtol, atol=atol)


def made_up_curves(shift=0):
    """A curve table as `curve_table()` returns it, with made-up values.

    Parameters: shift (number added to every curve, to tell tables apart).
    Returns: DataFrame with the columns CURVE_COLUMNS and 5 penalties.
    """
    return pd.DataFrame({
        "lambdas": [-1, -0.5, 0, 500, 1000],
        "ave_logpval": shift + np.array([1, 1.5, 2, 1.2, 1]),
        "ave_logpval_MOM": shift + np.array([1.1, 1.4, 2.1, 1.3, 0.9]),
        "power": [0.2, 0.4, 0.8, 0.5, 0.3], "power_MOM": [0.25, 0.45, 0.75, 0.5, 0.3],
        "gcv": shift + np.array([9, 4, 1, 2, 3]), "loo": shift + np.array([8, 5, 1.5, 2.5, 3.5])})


def write_eeg_file(eeg_dir, k, n=12, p=5):
    """Write a made-up file dat_ch<k>.csv in the layout of the EEG files: a
    row number, the coherence, then the SNPs.

    Parameters: eeg_dir (folder of the file), k (triple number; also the
    seed of the values), n, p (numbers of subjects and SNPs).
    Returns: the values written, an array with 2 + p columns.
    """
    values = np.column_stack([np.arange(1, n + 1), np.round(seeded_matrix(n, 1 + p, seed=k), 3)])
    names = ["id", "coherence"] + [f"snp{j}" for j in range(1, p + 1)]
    write_csv(pd.DataFrame(values, columns=names), os.path.join(eeg_dir, f"dat_ch{k}.csv"))
    return values


def write_liver_file(file, n=12, p=4):
    """Write a made-up file in the layout of mixOmicsData.csv: a row label,
    the 10 clinical variables, then the genes.

    Parameters: file (path of the file), n, p (numbers of animals and genes).
    Returns: the values written, an array with 10 + p columns.
    """
    values = np.round(seeded_matrix(n, 10 + p, seed=41), 3)
    names = [f"clinic{j}" for j in range(1, 11)] + [f"gene{j}" for j in range(1, p + 1)]
    table = pd.DataFrame(values, columns=names)
    table.insert(0, "id", [f"ID{i}" for i in range(1, n + 1)])
    write_csv(table, file)
    return values
