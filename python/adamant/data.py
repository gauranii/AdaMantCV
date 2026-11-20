# -----------------------------------------------------------------------------
# File:    python/adamant/data.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Table 4, Figures 1, 2 and 5 to 8.
# Purpose: Read the EEG coherence / SNP data and the liver toxicity data,
#          label the EEG triples, and return the curve tables of the EEG
#          triples and the liver variables (archived or recomputed).
# Inputs:  <data>/04_coherence_SNPs/dat_ch<k>.csv, data/liver/mixOmicsData.csv,
#          reference/triple.csv, reference/chan_names.csv,
#          archive/coh7/curves_eeg<k>.csv.gz, archive/liver/curves_liver<j>.csv.gz.
# Outputs: None.
# Usage:   from adamant.data import read_eeg_triple
#          read_eeg_triple(184, PATHS["eeg"])
# Source:  R/data.R.
# -----------------------------------------------------------------------------
import os

import numpy as np
import pandas as pd

from . import rbase
from .power_curve import curve_table
from .setup import read_csv

EEG_BANDS = ["delta", "theta", "alpha", "beta", "gamma"]
LIVER_LABELS = ["BUN", "Creatinine", "Total Protein", "Albumin", "ALT", "SDH",
                "AST", "ALP", "TBA", "Cholesterol"]
EEG_GRID = {"negative": 500, "positive": 2000}   # Dec 2021 EEG runs
LIVER_GRID = {"negative": 50, "positive": 250}   # Dec 2021 liver runs


def read_eeg_triple(k, eeg_dir):
    """EEG coherence and SNPs of one triple (channel pair and frequency band).

    Parameters: k (triple number), eeg_dir (folder with the files
    dat_ch<k>.csv).
    Returns: dict with `X` (standardized SNPs, 497 columns) and `y`
    (standardized coherence).
    """
    data = read_csv(os.path.join(eeg_dir, f"dat_ch{k}.csv")).to_numpy(dtype=np.float64)
    return {"X": rbase.scale(data[:, 2:]), "y": rbase.scale(data[:, 1])[:, 0].copy()}


def read_liver(file):
    """Liver toxicity data.

    Parameters: file (path of mixOmicsData.csv).
    Returns: dict with `X` (standardized gene expressions, 3116 columns) and
    `Y` (standardized clinical variables, 10 columns).
    """
    data = read_csv(file)
    return {"X": rbase.scale(data.iloc[:, 11:].to_numpy(dtype=np.float64)),
            "Y": rbase.scale(data.iloc[:, 1:11].to_numpy(dtype=np.float64))}


def triple_info(triples, reference_dir):
    """Channel names and band of EEG triples.

    Parameters: triples (triple numbers), reference_dir (folder with
    triple.csv and chan_names.csv).
    Returns: DataFrame with columns triple, channel1, channel2, band (for
    example Pz, PO4, alpha for triple 160).
    """
    triple = read_csv(os.path.join(reference_dir, "triple.csv"), index_col=0).loc[list(triples)]
    channels = [name.strip() for name in
                read_csv(os.path.join(reference_dir, "chan_names.csv"), header=None)[0]]
    return pd.DataFrame({"triple": list(triples),
                         "channel1": [channels[i - 1] for i in triple["chan1"]],
                         "channel2": [channels[i - 1] for i in triple["chan2"]],
                         "band": [EEG_BANDS[i - 1] for i in triple["band"]]})


def archived_curve_file(kind, k, archive_dir):
    """File of the curve table of the archived run of one response.

    Parameters: kind ("eeg" or "liver"), k (triple number for eeg, or number
    of the clinical variable for liver), archive_dir (the archive folder).
    Returns: path of the csv.gz file.
    """
    return os.path.join(archive_dir, "coh7" if kind == "eeg" else "liver",
                        f"curves_{kind}{k}.csv.gz")


def eeg_curves(triples, opts, settings, paths):
    """Curve tables of EEG triples, from the archive or recomputed.

    Parameters: triples (triple numbers), opts (options: rerun and those of
    `curve_table()`), settings (settings of the paper), paths (project
    paths).
    Returns: dict of curve tables keyed "<channel1>-<channel2>-<band>".
    """
    info = triple_info(triples, paths["reference"])
    curves = {}
    for k, channel1, channel2, band in info.itertuples(index=False):
        if not opts["rerun"]:
            table = read_csv(archived_curve_file("eeg", k, paths["archive"]))
        else:
            data = read_eeg_triple(k, paths["eeg"])
            table = curve_table(f"eeg{k}", data["X"], data["y"], opts, settings,
                                paths["results"])
        curves[f"{channel1}-{channel2}-{band}"] = table
    return curves


def liver_curves(opts, settings, paths):
    """Curve tables of the ten liver clinical variables, from the archive or
    recomputed.

    Parameters: opts, settings, paths (as for `eeg_curves()`).
    Returns: dict of curve tables keyed by clinical variable.
    """
    if opts["rerun"]:
        data = read_liver(paths["liver"])
    curves = {}
    for j, label in enumerate(LIVER_LABELS, start=1):
        if not opts["rerun"]:
            curves[label] = read_csv(archived_curve_file("liver", j, paths["archive"]))
        else:
            curves[label] = curve_table(f"liver{j}", data["X"], data["Y"][:, j - 1].copy(),
                                        opts, settings, paths["results"])
    return curves
