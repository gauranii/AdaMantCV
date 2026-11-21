# -----------------------------------------------------------------------------
# File:    python/tests/test_data.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Unit tests of python/adamant/data.py: the readers on small made-up files and on
#          the liver data, the labels of the EEG triples, and the curve
#          tables from the archived runs and from a small rerun; the cases
#          of tests/testthat/test-data.R.
# Inputs:  data/liver/mixOmicsData.csv, reference/, archive/coh7/, archive/liver/;
#          made-up EEG files are written to temporary folders (the EEG data
#          are not used).
# Outputs: None.
# Usage:   python -m pytest python/tests/test_data.py     # inside `nix develop`
# -----------------------------------------------------------------------------
import os

import numpy as np

from adamant.data import (EEG_GRID, LIVER_GRID, LIVER_LABELS, archived_curve_file, eeg_curves,
                          liver_curves, read_eeg_triple, read_liver, triple_info)
from adamant.power_curve import CURVE_COLUMNS, curve_table
from adamant.setup import PAPER, PATHS, REPO_ROOT, project_paths
from inputs import SMALL_OPTS, SMALL_SETTINGS, close, write_eeg_file, write_liver_file


def standardized(values):
    """A vector centred by its mean and divided by its sd."""
    return (values - values.mean()) / values.std(ddof=1)


def test_read_eeg_triple_standardizes_the_coherence_and_the_snps(tmp_path):
    values = write_eeg_file(str(tmp_path), k=7)
    triple = read_eeg_triple(7, str(tmp_path))
    assert triple["X"].shape == (12, 5)
    assert triple["y"].shape == (12,)
    close(triple["y"], standardized(values[:, 1]))
    close(triple["X"][:, 2], standardized(values[:, 4]))
    close(triple["X"].mean(axis=0), np.zeros(5))
    close(triple["X"].std(axis=0, ddof=1), np.ones(5))


def test_read_liver_splits_a_file_into_standardized_clinical_variables_and_genes(tmp_path):
    file = str(tmp_path / "liver.csv")
    values = write_liver_file(file)
    liver = read_liver(file)
    assert liver["Y"].shape == (12, 10)
    assert liver["X"].shape == (12, 4)
    close(liver["Y"][:, 9], standardized(values[:, 9]))
    close(liver["X"][:, 0], standardized(values[:, 10]))


def test_read_liver_reads_the_liver_toxicity_data():
    liver = read_liver(PATHS["liver"])
    assert liver["X"].shape == (64, 3116)
    assert liver["Y"].shape == (64, 10)
    close(liver["Y"].mean(axis=0), np.zeros(10))
    close(liver["X"][:, :5].std(axis=0, ddof=1), np.ones(5))
    assert not np.isnan(liver["X"]).any()


def test_triple_info_names_the_channels_and_the_band():
    info = triple_info([160, 1], PATHS["reference"])
    assert list(info.columns) == ["triple", "channel1", "channel2", "band"]
    assert info["triple"].tolist() == [160, 1]
    assert info.iloc[0, 1:].tolist() == ["Pz", "PO4", "alpha"]
    assert info.iloc[1, 1:].tolist() == ["FC3", "CP2", "delta"]


def test_archived_curve_file_points_into_the_archive():
    assert archived_curve_file("eeg", 184, "/archive") == "/archive/coh7/curves_eeg184.csv.gz"
    assert archived_curve_file("liver", 3, "/archive") == "/archive/liver/curves_liver3.csv.gz"
    assert os.path.exists(archived_curve_file("eeg", 160, PATHS["archive"]))


def test_eeg_curves_reads_the_archived_curve_tables():
    curves = eeg_curves([160, 184], {"rerun": False}, PAPER, PATHS)
    assert list(curves) == ["Pz-PO4-alpha", "F3-O2-alpha"]
    for table in curves.values():
        assert list(table.columns) == CURVE_COLUMNS
        assert len(table) == EEG_GRID["negative"] + EEG_GRID["positive"]
        assert table["lambdas"].is_monotonic_increasing
        assert ((table["power"] >= 0) & (table["power"] <= 1)).all()


def test_eeg_curves_recomputes_the_curves_of_a_made_up_triple(tmp_path):
    paths = project_paths(REPO_ROOT, data=str(tmp_path / "data"), output=str(tmp_path / "out"))
    os.makedirs(paths["eeg"])
    write_eeg_file(paths["eeg"], k=160)
    curves = eeg_curves([160], SMALL_OPTS, SMALL_SETTINGS, paths)
    assert list(curves) == ["Pz-PO4-alpha"]
    triple = read_eeg_triple(160, paths["eeg"])
    expected = curve_table("eeg160", triple["X"], triple["y"], SMALL_OPTS, SMALL_SETTINGS,
                           str(tmp_path / "other"))
    close(curves["Pz-PO4-alpha"].to_numpy(), expected.to_numpy())
    assert os.path.exists(os.path.join(paths["results"], "gcv_eeg160_grid3x3.csv"))


def test_liver_curves_reads_the_archived_curve_tables_of_the_ten_variables():
    curves = liver_curves({"rerun": False}, PAPER, PATHS)
    assert list(curves) == LIVER_LABELS
    for table in curves.values():
        assert list(table.columns) == CURVE_COLUMNS
        assert len(table) == LIVER_GRID["negative"] + LIVER_GRID["positive"]


def test_liver_curves_recomputes_the_curves_from_a_made_up_file(tmp_path):
    paths = project_paths(REPO_ROOT, output=str(tmp_path / "out"))
    paths["liver"] = str(tmp_path / "liver.csv")
    write_liver_file(paths["liver"])
    curves = liver_curves(SMALL_OPTS, SMALL_SETTINGS, paths)
    assert list(curves) == LIVER_LABELS
    liver = read_liver(paths["liver"])
    expected = curve_table("liver10", liver["X"], liver["Y"][:, 9].copy(), SMALL_OPTS,
                           SMALL_SETTINGS, str(tmp_path / "other"))
    close(curves["Cholesterol"].to_numpy(), expected.to_numpy())
    assert len(os.listdir(paths["results"])) == 20
