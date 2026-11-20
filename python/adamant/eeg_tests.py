# -----------------------------------------------------------------------------
# File:    python/adamant/eeg_tests.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Table 4.
# Purpose: AdaMantCV tests of the association between SNPs and EEG coherence
#          for a set of triples, from the archived run or recomputed.
# Inputs:  archive/real0210a/table4_pvalues.csv; in a rerun
#          <data>/04_coherence_SNPs/dat_ch<k>.csv and reference/table4_cluster_lambda0.csv.
# Outputs: None.
# Usage:   from adamant.eeg_tests import eeg_tests
#          eeg_tests([11, 27], 0.5, opts, PAPER, PATHS)
# Source:  R/eeg_tests.R.
# -----------------------------------------------------------------------------
import os

import pandas as pd

from .adamant_cv import TESTS, adamant_cv, adaptive_pvalue
from .data import read_eeg_triple, triple_info
from .setup import parallel_map, read_csv


def eeg_tests_archived(triples, file):
    """Tests of EEG triples from the archived p-values at the chosen
    penalties.

    Parameters: triples (triple numbers), file (path of table4_pvalues.csv:
    per triple, the p-values at the chosen penalties, p_gcv, p_loo,
    p_ga_gcv, p_ga_loo, and the chosen penalties of the observed data,
    b = 1, and every permutation).
    Returns: DataFrame with one row per triple: the four adaptive p-values
    (named as `TESTS`) and the penalties chosen for the observed data.
    """
    archived = read_csv(file)
    rows = []
    for k in triples:
        chosen = archived[archived["triple"] == k]
        rows.append({**{test: adaptive_pvalue(chosen[f"p_{test}"].to_numpy())
                        for test in TESTS},
                     "lambda_gcv": chosen["lambda_gcv"].iloc[0],
                     "lambda_loo": chosen["lambda_loo"].iloc[0]})
    return pd.DataFrame(rows)


def eeg_tests_rerun(triples, q, opts, settings, paths):
    """Tests of EEG triples recomputed from the data.

    The p-values of the observed data at lambda = 0 are taken from the
    values saved by the original run (argument `observed_lambda0` of
    `adamant_cv()`).

    Parameters: triples (triple numbers), q (quantile below which singular
    values are set to zero), opts (options: perms, negative, positive, cores,
    parallel processes over the triples), settings (settings of the paper:
    seed, upper), paths (project paths: eeg, reference).
    Returns: DataFrame as returned by `eeg_tests_archived()`.
    """
    saved = read_csv(os.path.join(paths["reference"], "table4_cluster_lambda0.csv"))

    def test(k):
        data = read_eeg_triple(k, paths["eeg"])
        lambda0 = saved[saved["triple"] == k].iloc[0]
        ada = adamant_cv(data["X"], data["y"], q, opts["negative"], opts["positive"],
                         opts["perms"], settings["seed"], settings["upper"],
                         observed_lambda0={"perms_below": float(lambda0["perms_below"]),
                                           "gamma": float(lambda0["gamma"])})
        return {**ada["p_values"], "lambda_gcv": ada["lambda_gcv"][0],
                "lambda_loo": ada["lambda_loo"][0]}

    return pd.DataFrame(parallel_map(test, triples, opts["cores"]))


def eeg_tests(triples, q, opts, settings, paths):
    """Table of the tests of EEG triples, from the archive or recomputed.

    Parameters: triples, q, opts, settings, paths (as for
    `eeg_tests_rerun()`; `opts["rerun"]` chooses between the two).
    Returns: DataFrame with the labels of the triples (`triple_info()`), the
    four adaptive p-values and the penalties chosen for the observed data.
    """
    if opts["rerun"]:
        tests = eeg_tests_rerun(triples, q, opts, settings, paths)
    else:
        tests = eeg_tests_archived(triples, os.path.join(paths["archive"], "real0210a",
                                                         "table4_pvalues.csv"))
    return pd.concat([triple_info(triples, paths["reference"]), tests], axis=1)
