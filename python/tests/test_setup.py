# -----------------------------------------------------------------------------
# File:    python/tests/test_setup.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Unit tests of python/adamant/setup.py: the settings, the project paths, the
#          command line parser, the index and label helpers, the csv reader
#          and writer (numbers as R writes them), the csv cache and the
#          parallel map; the cases of tests/testthat/test-setup.R and those
#          of the helpers R itself supplies.
# Inputs:  None (files are written to temporary folders).
# Outputs: None.
# Usage:   python -m pytest python/tests/test_setup.py     # inside `nix develop`
# -----------------------------------------------------------------------------
import os
import re

import pandas as pd
import pytest

from adamant import adamant_cv, data, eeg_tests, liver_cv, plots, power_curve, ridge_cv
from adamant import rng, setup, simulation
from adamant.setup import PAPER, PATHS, REPO_ROOT


def test_setup_finds_the_repository_and_the_settings_of_the_paper():
    assert os.path.exists(os.path.join(REPO_ROOT, "flake.nix"))
    for module, name in [(adamant_cv, "adamant_cv"), (ridge_cv, "ridge_cv"),
                         (simulation, "simulation_study"), (power_curve, "power_curve"),
                         (data, "read_liver"), (eeg_tests, "eeg_tests"),
                         (liver_cv, "cv_curves"), (plots, "heat_map")]:
        assert callable(getattr(module, name))
    assert {key: PAPER[key] for key in ("seed", "alpha", "sigma", "upper", "n_perms")} == \
        {"seed": 8891, "alpha": 0.05, "sigma": 0.5, "upper": 1000, "n_perms": 1000}


def test_project_paths_puts_data_and_outputs_where_it_is_told():
    paths = setup.project_paths("/repo", data="/elsewhere", output="/out")
    assert paths == {"eeg": "/elsewhere/04_coherence_SNPs",
                     "liver": "/repo/data/liver/mixOmicsData.csv",
                     "reference": "/repo/reference", "archive": "/repo/archive",
                     "results": "/out/results", "figures": "/out/figures"}
    assert os.path.exists(PATHS["liver"])


def test_project_paths_reads_the_environment(monkeypatch):
    monkeypatch.setenv("ADAMANT_DATA", "/data")
    monkeypatch.setenv("ADAMANT_OUTPUT", "/work")
    paths = setup.project_paths("/repo")
    assert paths["eeg"] == "/data/04_coherence_SNPs" and paths["results"] == "/work/results"


def test_output_file_creates_the_folder_and_returns_the_path(tmp_path):
    folder = str(tmp_path / "a" / "b")
    assert setup.output_file(folder, "x.csv") == os.path.join(folder, "x.csv")
    assert os.path.isdir(folder)


def test_parse_args_converts_values_to_the_type_of_their_default():
    defaults = {"runs": "500", "perms": 1000.0, "cores": 4, "fast": False}
    assert setup.parse_args(defaults, []) == {"rerun": False, **defaults}
    opts = setup.parse_args(defaults, ["--rerun", "--runs=1:100", "--perms=50", "--cores=2",
                                       "--fast"])
    assert opts == {"rerun": True, "runs": "1:100", "perms": 50.0, "cores": 2, "fast": True}
    assert isinstance(opts["cores"], int) and isinstance(opts["perms"], float)


@pytest.mark.parametrize("argument, message", [
    ("--other=1", "unknown argument --other"), ("perms=5", "must look like"),
    ("--perms", "--perms needs a value"), ("--fast=1", "--fast takes no value")])
def test_parse_args_refuses_what_it_does_not_know(argument, message):
    with pytest.raises(SystemExit, match=message):
        setup.parse_args({"perms": 1000, "fast": False}, [argument])


def test_curve_options_and_test_options_are_the_defaults_of_the_scripts():
    assert setup.curve_options({"negative": 50, "positive": 250}, PAPER) == \
        {"runs": 1000, "perms": 1000, "negative": 50, "positive": 250, "cores": 4}
    assert setup.test_options(PAPER) == {"perms": 1000, "negative": 50, "positive": 200,
                                         "cores": 4}
    assert setup.test_options(PAPER, cores=5)["cores"] == 5


def test_parse_index_set_reads_ranges_lists_and_counts():
    assert setup.parse_index_set("3:7") == [3, 4, 5, 6, 7]
    assert setup.parse_index_set("2,5,9") == [2, 5, 9]
    assert setup.parse_index_set("4") == [4]
    assert setup.parse_index_set("4", count=True) == [1, 2, 3, 4]
    assert setup.parse_index_set("2,5", count=True) == [2, 5]


def test_settings_tag_names_the_settings_of_a_rerun():
    opts = {"perms": 1000, "negative": 50, "positive": 200}
    assert setup.settings_tag(opts) == "perms1000_grid50x200"
    assert setup.settings_tag(opts, n_runs=25) == "runs25_perms1000_grid50x200"


@pytest.mark.parametrize("number, text", [
    (0.0, "0"), (1.0, "1"), (-2.5, "-2.5"), (1000.0, "1000"), (0.1, "0.1"),
    (1 / 3, "0.333333333333333"), (2 / 3, "0.666666666666667"), (123456.789, "123456.789"),
    (1e5, "1e+05"), (100000.5, "100000.5"), (1e-5, "1e-05"), (0.0001, "1e-04"),
    (1e15, "1e+15"), (1.5e-10, "1.5e-10"), (0.993000000000001, "0.993000000000001"),
    (1.4901161193847656e-08, "1.49011611938477e-08"), (float("nan"), "NA"),
    (float("inf"), "Inf"), (float("-inf"), "-Inf")])
def test_format_number_writes_a_number_as_r_does(number, text):
    # the texts are what R's write.csv() writes (15 significant digits)
    assert setup.format_number(number) == text


def test_write_csv_and_read_csv_keep_every_type(tmp_path):
    file = str(tmp_path / "table.csv")
    table = pd.DataFrame({"count": [1, 2], "value": [0.5, 1 / 3], "flag": [True, False],
                          "label": ["a", 'b "quoted"']})
    setup.write_csv(table, file)
    with open(file) as handle:
        assert handle.read() == ('"count","value","flag","label"\n'
                                 '1,0.5,TRUE,"a"\n'
                                 '2,0.333333333333333,FALSE,"b ""quoted"""\n')
    back = setup.read_csv(file)
    assert list(back["count"]) == [1, 2] and list(back["label"]) == ["a", 'b "quoted"']
    assert back["value"][1] == 0.333333333333333


def test_read_csv_reads_every_digit(tmp_path):
    file = str(tmp_path / "digits.csv")
    with open(file, "w") as handle:
        handle.write("x\n0.1\n4.8335366933517551e-05\n")
    assert list(setup.read_csv(file)["x"]) == [0.1, 4.8335366933517551e-05]


def test_cached_csv_computes_once_and_then_reads_the_file(tmp_path, capsys):
    file = str(tmp_path / "cache.csv")
    calls = []

    def compute():
        calls.append(1)
        return pd.DataFrame({"x": [1.0, 2.5], "label": ["a", "b"]})

    first = setup.cached_csv(file, compute)
    assert os.path.exists(file)
    second = setup.cached_csv(file, compute)
    assert len(calls) == 1
    assert "reusing" in capsys.readouterr().out
    pd.testing.assert_frame_equal(second, first)


def _square(x):
    return x * x


def test_parallel_map_keeps_the_order_for_any_number_of_cores():
    assert setup.parallel_map(_square, [3, 1, 2], cores=1) == [9, 1, 4]
    assert setup.parallel_map(_square, range(5), cores=2) == [0, 1, 4, 9, 16]
    assert setup.parallel_map(lambda x: x + 1, [1, 2, 3], cores=3) == [2, 3, 4]
    assert setup.parallel_map(_square, [], cores=2) == []


def test_parallel_map_leaves_the_stream_of_this_process_alone():
    # a forked process starts with a copy of the stream and may run several
    # jobs, so a job that draws must seed itself, as every job of adamant does
    rng.set_seed(1)
    expected = rng.runif(2).tolist()
    rng.set_seed(1)
    draws = setup.parallel_map(lambda i: float(rng.runif(1)[0]), [1, 2], cores=2)
    assert set(draws) <= set(expected)
    assert float(rng.runif(1)[0]) == expected[0]

    def seeded(i):
        rng.set_seed(i)
        return float(rng.runif(1)[0])

    assert setup.parallel_map(seeded, [1, 1, 2], cores=2) == \
        setup.parallel_map(seeded, [1, 1, 2], cores=1)


def test_log_message_prints_the_time_and_the_formatted_message(capsys):
    assert setup.log_message("%d of %s", 3, "five") is None
    assert re.match(r"[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9:]{8} - 3 of five", capsys.readouterr().out)
