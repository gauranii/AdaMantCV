# -----------------------------------------------------------------------------
# File:    python/tests/test_simulation.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Unit tests of python/adamant/simulation.py: the covariance builders, the simulated
#          designs, the scenario numbers, one replication, the stored
#          replications of a small rerun, the rejection rates of a hand-made
#          table, and the tables read from the archived runs; the cases of
#          tests/testthat/test-simulation.R.
# Inputs:  archive/run0204b/ (p-values and chosen penalties of the archived runs);
#          files are written to temporary folders.
# Outputs: None.
# Usage:   python -m pytest python/tests/test_simulation.py     # inside `nix develop`
# -----------------------------------------------------------------------------
import os

import numpy as np
import pandas as pd
import pytest

from adamant import rbase, rng
from adamant.adamant_cv import TESTS, adamant_cv
from adamant.setup import PAPER, REPO_ROOT, project_paths, read_csv
from adamant.simulation import (SIMULATION_DESIGN, chosen_penalties, make_covariance,
                                rejection_rates, run_replications, simulate_design,
                                simulate_designs, simulation_run, simulation_scenarios,
                                simulation_study)
from inputs import SMALL_OPTS, SMALL_SETTINGS, close

# A design with 12 observations and 4 predictors, in place of
# SIMULATION_DESIGN.
SMALL_DESIGN = {"n": 12, "rho": 0.025, "beta_mean": 1,
                "quantiles": {"none": 0, "Q1": 0.25, "Q2": 0.5},
                "covariances": {"CS": "CSYM", "Heteroskedastic": "NCDIAG"},
                "dimensions": [4]}
P_COLUMNS = ["p_gcv", "p_loo", "p_ga_gcv", "p_ga_loo"]


def test_make_covariance_builds_the_two_structures():
    close(make_covariance("CSYM", 3, 0.025),
          [[1, 0.025, 0.025], [0.025, 1, 0.025], [0.025, 0.025, 1]])
    close(make_covariance("NCDIAG", 3, 0.025), np.diag(np.log([2, 3, 4])))
    assert make_covariance("NCDIAG", 1, 0).shape == (1, 1)
    with pytest.raises(ValueError, match="unknown covariance type AR1"):
        make_covariance("AR1", 3, 0.5)


def test_simulate_design_returns_a_standardized_design_and_the_coefficients():
    fixed = simulate_design(30, 4, "CSYM", 0.025, "fixed", mu=1, sigma=0.5, seed=8891)
    assert fixed["X"].shape == (30, 4)
    close(fixed["X"].mean(axis=0), np.zeros(4))
    close(fixed["X"].std(axis=0, ddof=1), np.ones(4))
    close(fixed["beta"], np.ones(4))
    close(simulate_design(30, 4, "CSYM", 0.025, "zero", 1, 0.5, 8891)["beta"], np.zeros(4))
    with pytest.raises(ValueError, match="unknown beta type other"):
        simulate_design(30, 4, "CSYM", 0.025, "other", 1, 0.5, 8891)


def test_simulate_design_depends_on_its_seed_only():
    fixed = simulate_design(30, 4, "CSYM", 0.025, "fixed", 1, 0.5, seed=8891)
    rng.runif(5)
    random = simulate_design(30, 4, "CSYM", 0.025, "random", 1, 0.5, seed=8891)
    assert np.array_equal(random["X"], fixed["X"])
    other = simulate_design(30, 4, "CSYM", 0.025, "fixed", 1, 0.5, 8892)
    assert not np.array_equal(other["X"], fixed["X"])
    # the random coefficients are the draws that follow the design in the stream
    rng.set_seed(8891)
    rbase.mvrnorm(30, make_covariance("CSYM", 4, 0.025))
    close(random["beta"], rng.rnorm(4, 1, 0.5))


def test_simulation_scenarios_numbers_the_12_settings_as_the_original_runs():
    null = simulation_scenarios(False, SIMULATION_DESIGN)
    assert list(null.columns) == ["scenario", "covariance", "p", "tau"]
    assert null["scenario"].tolist() == [2, 5, 8, 11, 14, 17, 3, 6, 9, 12, 15, 18]
    assert null["covariance"].tolist() == ["CS"] * 6 + ["Heteroskedastic"] * 6
    assert null["p"].tolist() == [500, 500, 500, 1000, 1000, 1000] * 2
    assert null["tau"].tolist() == ["none", "Q1", "Q2"] * 4
    fixed = simulation_scenarios(True, SIMULATION_DESIGN)
    assert fixed["scenario"].tolist() == (null["scenario"] + 18).tolist()


def test_simulation_run_draws_the_response_with_seed_10_times_run_plus_seed():
    simulated = simulate_design(20, 5, "CSYM", 0.025, "fixed", 1, 0.5, seed=8891)
    ada = simulation_run(simulated, 3, 0.25, 0.5, 3, 3, 9, seed=8891)
    rng.set_seed(10 * 3 + 8891)
    Y = (simulated["X"] @ simulated["beta"]) + rng.rnorm(20, 0, 0.5)
    direct = adamant_cv(simulated["X"], Y, 0.25, 3, 3, 9, 8891)
    assert list(ada["p_values"]) == TESTS
    assert ada["p_values"] == direct["p_values"]
    close(ada["statistics"], direct["statistics"])
    other = simulation_run(simulated, 4, 0.25, 0.5, 3, 3, 9, seed=8891)
    assert not np.array_equal(other["statistics"], ada["statistics"])


def test_simulate_designs_makes_one_design_per_covariance_and_dimension():
    design = {**SMALL_DESIGN, "dimensions": [4, 6]}
    scenarios = simulation_scenarios(False, design)
    simulated = simulate_designs(scenarios, "zero", SMALL_SETTINGS, design)
    assert list(simulated) == ["CS 4", "CS 6", "Heteroskedastic 4", "Heteroskedastic 6"]
    assert simulated["Heteroskedastic 6"]["X"].shape == (12, 6)
    assert np.array_equal(simulated["CS 4"]["X"],
                          simulate_design(12, 4, "CSYM", 0.025, "zero", 1, 0.5, 8891)["X"])


def test_run_replications_stores_one_file_per_replication_and_reuses_them(tmp_path, capsys):
    run_dir = str(tmp_path / "runs")
    scenarios = simulation_scenarios(False, SMALL_DESIGN).iloc[:2]
    opts = {**SMALL_OPTS, "runs": "2"}
    pvalues = run_replications(run_dir, "zero", scenarios, opts, SMALL_SETTINGS, SMALL_DESIGN)
    assert "4 replications to run" in capsys.readouterr().out
    assert list(pvalues.columns) == ["scenario", "run"] + P_COLUMNS
    assert pvalues["scenario"].tolist() == [2, 2, 5, 5]
    assert pvalues["run"].tolist() == [1, 2, 1, 2]
    assert sorted(os.listdir(run_dir)) == ["scenario02_run0001.csv", "scenario02_run0002.csv",
                                           "scenario05_run0001.csv", "scenario05_run0002.csv"]
    assert ((pvalues[P_COLUMNS] >= 0) & (pvalues[P_COLUMNS] <= 1)).all().all()

    # a replication is the test of simulation_run() on the design of its setting
    simulated = simulate_design(12, 4, "CSYM", 0.025, "zero", 1, 0.5, 8891)
    ada = simulation_run(simulated, 2, 0.25, 0.5, 3, 3, 9, 8891)
    stored = pvalues[(pvalues["scenario"] == 5) & (pvalues["run"] == 2)].iloc[0]
    close(stored[P_COLUMNS], list(ada["p_values"].values()))

    # stored replications are not run again; new ones are added
    opts["runs"] = "3"
    more = run_replications(run_dir, "zero", scenarios, opts, SMALL_SETTINGS, SMALL_DESIGN)
    assert "2 replications to run" in capsys.readouterr().out
    assert len(more) == 6
    pd.testing.assert_frame_equal(more[more["run"] <= 2].reset_index(drop=True), pvalues)


def test_rejection_rates_counts_the_p_values_below_alpha_per_setting():
    scenarios = pd.DataFrame({"scenario": [5, 2], "covariance": "CS", "p": [4, 4],
                              "tau": ["Q1", "none"]})
    pvalues = pd.DataFrame({"scenario": [2, 2, 2, 2, 5, 5], "run": [1, 2, 3, 4, 1, 2],
                            "p_gcv": [0.01, 0.04, 0.05, 0.9, 0.5, 0.6],
                            "p_loo": [0.01, 0.2, 0.3, 0.9, 0.01, 0.6],
                            "p_ga_gcv": [0.01, 0.02, 0.03, 0.04, 0.5, 0.6],
                            "p_ga_loo": [0.5, 0.5, 0.5, 0.5, 0.01, 0.02]})
    rates = rejection_rates(scenarios, pvalues, alpha=0.05)
    assert rates["scenario"].tolist() == [5, 2]      # the order of `scenarios`
    assert rates["runs"].tolist() == [2, 4]
    assert rates["gcv"].tolist() == [0, 0.5]         # 0.05 itself is not below alpha
    assert rates["loocv"].tolist() == [0.5, 0.25]
    assert rates["ga_gcv"].tolist() == [0, 1]
    assert rates["ga_loocv"].tolist() == [1, 0]
    assert rates["tau"].tolist() == ["Q1", "none"]


def test_simulation_study_builds_table_1_from_the_archived_p_values(tmp_path, capsys):
    paths = project_paths(REPO_ROOT, output=str(tmp_path))
    opts = {"rerun": False, "scenarios": "all"}
    table = simulation_study("Table1_type1-error", "zero", opts, PAPER, paths)
    assert table.shape == (12, 9)
    assert table["scenario"].tolist() == \
        simulation_scenarios(False, SIMULATION_DESIGN)["scenario"].tolist()
    assert (table["runs"] == 500).all()
    # under the null every rejection rate is near the level 0.05
    rates = table[["gcv", "loocv", "ga_gcv", "ga_loocv"]].to_numpy()
    assert ((rates >= 0.01) & (rates <= 0.1)).all()
    written = read_csv(os.path.join(paths["results"], "Table1_type1-error.csv"))
    pd.testing.assert_frame_equal(written, table)
    assert "Heteroskedastic" in capsys.readouterr().out

    opts["scenarios"] = "2,17"
    chosen = simulation_study("Table1_type1-error", "zero", opts, PAPER, paths)
    assert chosen["scenario"].tolist() == [2, 17]


def test_simulation_study_reruns_a_small_design(tmp_path):
    paths = project_paths(REPO_ROOT, output=str(tmp_path))
    opts = {**SMALL_OPTS, "runs": "2", "scenarios": "20,23"}
    table = simulation_study("Table2_power_fixed-effects", "fixed", opts, SMALL_SETTINGS, paths,
                             design=SMALL_DESIGN)
    assert table["scenario"].tolist() == [20, 23]
    assert table["runs"].tolist() == [2, 2]
    assert os.path.isdir(os.path.join(paths["results"],
                                      "Table2_power_fixed-effects_perms9_grid3x3"))


def test_chosen_penalties_reads_the_archive_or_chooses_the_penalties_again(tmp_path):
    panels = pd.DataFrame({"panel": ["a", "b"], "beta_type": ["zero", "fixed"],
                           "tau": ["Q1", "Q2"]})
    paths = project_paths(REPO_ROOT, output=str(tmp_path))
    archived = chosen_penalties(panels, 1, 500, "CSYM", {"rerun": False}, PAPER, paths)
    assert list(archived.columns) == ["panel", "b", "GCV", "LOOCV"]
    assert archived["panel"].value_counts(sort=False).to_dict() == \
        {"a": 1001, "b": 1001, "c": 1001, "d": 1001}

    chosen = chosen_penalties(panels, 1, 5, "CSYM", SMALL_OPTS, SMALL_SETTINGS, paths,
                              design=SMALL_DESIGN)
    assert list(chosen.columns) == ["panel", "beta_type", "tau", "b", "GCV", "LOOCV"]
    assert chosen["panel"].tolist() == ["a"] * 10 + ["b"] * 10
    assert chosen["b"].tolist() == list(range(1, 11)) * 2
    simulated = simulate_design(12, 5, "CSYM", 0.025, "fixed", 1, 0.5, 8891)
    ada = simulation_run(simulated, 1, 0.5, 0.5, 3, 3, 9, 8891)
    close(chosen["GCV"][chosen["panel"] == "b"], ada["lambda_gcv"])
    assert os.path.exists(os.path.join(paths["results"], "fig3_lambdas_perms9_grid3x3.csv"))
