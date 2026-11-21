# -----------------------------------------------------------------------------
# File:    python/tests/test_plots.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Unit tests of python/adamant/plots.py: every plot function returns a figure whose
#          recorded marks are the expected ones, the save functions write
#          the files they name, and the colours and bins are those of
#          ggplot2; the cases of tests/testthat/test-plots.R and those of
#          the helpers ggplot2 supplies in R.
# Inputs:  None (figures are written to temporary folders).
# Outputs: None.
# Usage:   python -m pytest python/tests/test_plots.py     # inside `nix develop`
# -----------------------------------------------------------------------------
import os

import numpy as np
import pandas as pd

from adamant import plots
from adamant.setup import read_csv
from inputs import close, made_up_curves

CURVES = {"first": made_up_curves(), "second": made_up_curves(shift=1)}


def test_plot_keeps_the_numbers_of_its_layers():
    plot = plots.Plot(theme="bw")
    plot.record(1, x=[1.0, 2.0], y=[3.0, 4.0], group=1)
    plot.record(2, xintercept=[5.0])
    marks = plot.plotted_data()
    assert list(marks.columns) == plots.PLOT_DATA_COLUMNS
    assert marks["layer"].tolist() == [1, 1, 2]
    assert marks["x"].tolist()[:2] == [1.0, 2.0] and np.isnan(marks["x"][2])
    assert marks["xintercept"][2] == 5.0


def test_save_figure_writes_a_square_pdf_and_the_plotted_numbers(tmp_path, monkeypatch, capsys):
    figures = str(tmp_path / "figures")
    monkeypatch.delenv("ADAMANT_PLOT_DATA", raising=False)
    plots.save_figure(plots.sign_scatter(made_up_curves(), "power", "power"), "scatter.pdf",
                      figures, 3)
    assert "wrote figures/scatter.pdf" in capsys.readouterr().out
    with open(os.path.join(figures, "scatter.pdf"), "rb") as handle:
        assert handle.read(4) == b"%PDF"
    assert os.listdir(figures) == ["scatter.pdf"]

    monkeypatch.setenv("ADAMANT_PLOT_DATA", str(tmp_path / "numbers"))
    plot = plots.sign_scatter(made_up_curves(), "power", "power")
    plots.save_figure(plot, "scatter.pdf", figures, 3)
    assert plot.figure.get_size_inches().tolist() == [3, 3]
    written = read_csv(str(tmp_path / "numbers" / "scatter.csv"))
    assert list(written.columns) == plots.PLOT_DATA_COLUMNS
    close(written["y"], made_up_curves()["power"])


def test_hue_colours_are_the_default_colours_of_ggplot2():
    assert plots.hue_colours(1) == ["#F8766D"]
    assert plots.hue_colours(2) == ["#F8766D", "#00BFC4"]
    assert plots.hue_colours(3) == ["#F8766D", "#00BA38", "#619CFF"]


def test_diverging_colours_run_from_navy_through_white_to_red():
    colours = plots.diverging_colours()
    assert colours.N == 511
    close(colours(0.0)[:3], [0, 0, 128 / 255], atol=1e-5)
    close(colours(0.5)[:3], [1, 1, 1], atol=1e-5)
    close(colours(1.0)[:3], [1, 0, 0], atol=1e-5)
    blue_side, red_side = colours(0.25), colours(0.75)
    assert blue_side[2] > blue_side[0] and red_side[0] > red_side[2]


def test_r_round_rounds_to_decimals():
    close(plots.r_round([1.2345, -0.5004, 1000.0], 2), [1.23, -0.5, 1000])
    close(plots.r_round([0.577, 2.0], 0), [1, 2])


def test_heat_map_has_one_tile_per_table_and_penalty_and_the_line_at_zero():
    plot = plots.heat_map(CURVES, "power", 0.4, 3, 3, 10)
    assert isinstance(plot, plots.Plot)
    marks = plot.plotted_data()
    tiles = marks[marks["layer"] == 1]
    assert len(tiles) == 10
    assert tiles["y"].value_counts().to_dict() == {1: 5, 2: 5}   # one row of tiles per table
    assert sorted(tiles["x"].unique()) == [1, 2, 3, 4, 5]
    close(tiles["value"], np.tile(made_up_curves()["power"], 2))
    assert marks[marks["layer"] == 2]["xintercept"].tolist() == [3]
    assert plot.ax.get_xlabel() == "lambda"
    assert [label.get_text() for label in plot.ax.get_yticklabels()] == ["first", "second"]


def test_heat_map_rounds_the_penalties_so_close_penalties_share_a_column():
    curves = {"only": made_up_curves().assign(lambdas=[-1, -0.5004, -0.5001, 0, 1000])}
    marks = plots.heat_map(curves, "gcv", 2, 2, 3, 10).plotted_data()
    assert sorted(marks[marks["layer"] == 1]["x"].unique()) == [1, 2, 3, 4]


def test_save_heat_maps_writes_panels_a_and_b_with_their_names(tmp_path):
    plots.save_heat_maps("Fig9", "EEG", "_test", CURVES, "power_MOM", {"logp": 1.5, "power": 0.4},
                         3, 3, 10, str(tmp_path), 3)
    assert sorted(os.listdir(tmp_path)) == ["Fig9A_EEG_heatmap_mean-log10-pvalue_test.pdf",
                                            "Fig9B_EEG_heatmap_power_test.pdf"]


def test_sign_scatter_has_one_point_per_penalty_coloured_by_its_sign():
    table = made_up_curves()
    plot = plots.sign_scatter(table, "gcv", "Generalized Cross Validation")
    assert isinstance(plot, plots.Plot)
    points = plot.plotted_data()
    close(points["x"], table["lambdas"])
    close(points["y"], table["gcv"])
    colours = points["colour"].tolist()
    assert len(set(colours)) == 3
    assert colours[0] == colours[1] and colours[3] == colours[4] and colours[2] != colours[0]
    assert plot.ax.get_ylabel() == "Generalized Cross Validation"


def test_save_power_gcv_panels_writes_a_power_and_a_gcv_panel_per_triple(tmp_path):
    windows = pd.DataFrame({"triple": [1, 2], "first": [2, 1], "last": [4, 99],
                            "y_max": [3, 100]})
    plots.save_power_gcv_panels("Fig9", "rising", windows, CURVES, str(tmp_path), 3)
    assert sorted(os.listdir(tmp_path)) == [
        "Fig9A_EEG_power-rising_first.pdf", "Fig9B_EEG_power-rising_second.pdf",
        "Fig9C_EEG_GCV_first.pdf", "Fig9D_EEG_GCV_second.pdf"]


def test_histogram_bins_are_centred_on_the_multiples_of_the_width():
    bins = plots.histogram_bins(np.array([0.1, 0.9, 1.1, 2.0, 2.0, 5.2]), (0.1, 5.2), 1)
    close(bins["x"], [0, 1, 2, 3, 4, 5])
    close(bins["xmin"], np.arange(6) - 0.5)
    close(bins["xmax"], np.arange(6) + 0.5)
    assert bins["count"].tolist() == [1, 2, 2, 0, 0, 1]
    # values on the two edges of a bin both fall in it (as in ggplot2)
    edges = plots.histogram_bins(np.array([0.5, 1.5]), (0.5, 1.5), 1)
    assert edges["x"].tolist() == [1] and edges["count"].tolist() == [2]


def test_penalty_histogram_counts_every_chosen_penalty_once_per_method():
    chosen = pd.DataFrame({"GCV": [0.1, 0.9, 1.1, 2.0, 2.0, 5.2],
                           "LOOCV": [0.1, 0.1, 0.1, 5, 5, 5]})
    plot = plots.penalty_histogram(chosen, "Q1", 1)
    assert isinstance(plot, plots.Plot)
    marks = plot.plotted_data()
    assert sorted(marks["layer"].unique()) == [1, 2]
    # bins of width 1 centred on the integers
    gcv = marks[(marks["layer"] == 1) & (marks["count"] > 0)]
    assert gcv["x"].tolist() == [0, 1, 2, 5] and gcv["count"].tolist() == [1, 2, 2, 1]
    loocv = marks[(marks["layer"] == 2) & (marks["count"] > 0)]
    assert loocv["x"].tolist() == [0, 5] and loocv["count"].tolist() == [3, 3]
    labels = [text.get_text() for text in plot.figure.legends[0].get_texts()]
    assert labels == ["variable", "GCV tau = Q1", "LOOCV tau = Q1"]
