# -----------------------------------------------------------------------------
# File:    python/tests/test_eeg_tests.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Unit tests of python/adamant/eeg_tests.py: the adaptive p-values from a hand-made
#          file of p-values at the chosen penalties, from the archived run,
#          and from a rerun on made-up EEG files; the cases of
#          tests/testthat/test-eeg_tests.R.
# Inputs:  reference/, archive/real0210a/table4_pvalues.csv; made-up EEG files
#          are written to temporary folders (the EEG data are not used).
# Outputs: None.
# Usage:   python -m pytest python/tests/test_eeg_tests.py     # inside `nix develop`
# -----------------------------------------------------------------------------
import os

import numpy as np
import pandas as pd

from adamant.adamant_cv import TESTS, adamant_cv
from adamant.data import read_eeg_triple
from adamant.eeg_tests import eeg_tests, eeg_tests_archived, eeg_tests_rerun
from adamant.setup import PAPER, PATHS, REPO_ROOT, project_paths, read_csv, write_csv
from inputs import SMALL_OPTS, SMALL_SETTINGS, close, write_eeg_file

COLUMNS = TESTS + ["lambda_gcv", "lambda_loo"]


def test_eeg_tests_archived_ranks_the_observed_p_value_among_the_permutations(tmp_path):
    file = str(tmp_path / "pvalues.csv")
    write_csv(pd.DataFrame({"triple": [3] * 4 + [8] * 4, "b": [1, 2, 3, 4] * 2,
                            "p_gcv": [0.2, 0.1, 0.2, 0.4, 0.01, 0.5, 0.6, 0.7],
                            "p_loo": [0.9, 0.1, 0.2, 0.4, 0.5, 0.5, 0.5, 0.5],
                            "p_ga_gcv": [0.3, 0.1, 0.2, 0.4, 0.6, 0.5, 0.7, 0.1],
                            "p_ga_loo": [0.1, 0.1, 0.2, 0.4, 0.2, 0.3, 0.1, 0.4],
                            "lambda_gcv": [0, 1000, 5, 5, -0.5, 1, 1, 1],
                            "lambda_loo": [1000, 0, 5, 5, 2.0, 1, 1, 1]}), file)
    tests = eeg_tests_archived([8, 3], file)
    assert list(tests.columns) == COLUMNS
    close(tests.iloc[0], [0, 0, 2 / 3, 1 / 3, -0.5, 2])
    close(tests.iloc[1], [1 / 3, 1, 2 / 3, 0, 0, 1000])


def test_eeg_tests_builds_the_table_of_the_archived_run():
    table = eeg_tests([11, 248], 0.5, {"rerun": False}, PAPER, PATHS)
    assert list(table.columns) == ["triple", "channel1", "channel2", "band"] + COLUMNS
    assert table["triple"].tolist() == [11, 248]
    assert table["channel1"].tolist() == ["Fp1", "P4"]
    assert ((table[TESTS] >= 0) & (table[TESTS] <= 1)).all().all()
    # the observed p-value 0.015 of triple 11 is ranked among its 1000 permutations
    archived = read_csv(os.path.join(PATHS["archive"], "real0210a", "table4_pvalues.csv"))
    permuted = archived["p_gcv"][archived["triple"] == 11].to_numpy()[1:]
    close(table["gcv"][0], np.mean(permuted < 0.015))


def test_eeg_tests_rerun_tests_made_up_triples_with_the_saved_values_at_zero(tmp_path):
    paths = project_paths(REPO_ROOT, data=str(tmp_path))
    os.makedirs(paths["eeg"])
    for k in (11, 27):
        write_eeg_file(paths["eeg"], k)
    tests = eeg_tests_rerun([11, 27], 0.5, SMALL_OPTS, SMALL_SETTINGS, paths)
    assert tests.shape == (2, 6)
    assert list(tests.columns) == COLUMNS

    triple = read_eeg_triple(27, paths["eeg"])
    ada = adamant_cv(triple["X"], triple["y"], 0.5, 3, 3, 9, 8891, 1000,
                     observed_lambda0={"perms_below": 992, "gamma": 0.0118191573598137})
    close(tests.iloc[1], [*ada["p_values"].values(), ada["lambda_gcv"][0], ada["lambda_loo"][0]])
    close(ada["p_gamma"][0, ada["lambdas"] == 0], 0.0118191573598137)

    table = eeg_tests([11, 27], 0.5, SMALL_OPTS, SMALL_SETTINGS, paths)
    close(table[COLUMNS].to_numpy(), tests.to_numpy())
    assert table["band"].tolist() == ["delta", "delta"]
