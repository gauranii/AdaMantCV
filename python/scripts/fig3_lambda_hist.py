# -----------------------------------------------------------------------------
# File:    python/scripts/fig3_lambda_hist.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Figure 3.
# Purpose: Histograms of the ridge penalties chosen by GCV and LOOCV over the
#          permutations of one simulated data set (replication 500, compound
#          symmetric covariance, p = 1000), under beta = 0 (A, B) and beta = 1
#          (C, D), with tau = Q1 (A, C) and tau = Q2 (B, D).
# Inputs:  archive/run0204b/fig3_lambdas.csv (chosen penalties of the archived
#          run); none with --rerun (simulated data).
# Outputs: figures/Fig3{A,B}_simulation_lambda-histogram_null_threshold-{Q1,median}.pdf,
#          figures/Fig3{C,D}_simulation_lambda-histogram_alternative_threshold-{Q1,median}.pdf;
#          with --rerun also results/fig3_lambdas_<settings>.csv.
# Usage:   python python/scripts/fig3_lambda_hist.py [--rerun] [--perms=1000] [--negative=50]
#          [--positive=200] [--cores=4]
#          --rerun: simulate the data and choose the penalties here.
# Source:  scripts/fig3_lambda_hist.R.
# -----------------------------------------------------------------------------
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import pandas as pd  # noqa: E402

from adamant.plots import penalty_histogram, save_figure  # noqa: E402
from adamant.setup import PAPER, PATHS, parse_args, test_options  # noqa: E402
from adamant.simulation import chosen_penalties  # noqa: E402

REPLICATION = 500     # the replication shown
DIMENSION = 1000
COVARIANCE = "CSYM"
BIN_WIDTH = 0.25
PANELS = pd.DataFrame({"panel": ["a", "b", "c", "d"],
                       "beta_type": ["zero", "zero", "fixed", "fixed"],
                       "tau": ["Q1", "Q2", "Q1", "Q2"]})
# Words of the file names
HYPOTHESIS = {"zero": "null", "fixed": "alternative"}
THRESHOLD = {"Q1": "Q1", "Q2": "median"}

opts = parse_args(test_options(PAPER))

lambdas = chosen_penalties(PANELS, REPLICATION, DIMENSION, COVARIANCE, opts, PAPER, PATHS)

# The histograms show the permutations (b = 1 is the observed data).
lambdas = lambdas[lambdas["b"] >= 2]

for panel, beta_type, tau in PANELS.itertuples(index=False):
    hist = penalty_histogram(lambdas[lambdas["panel"] == panel], tau, BIN_WIDTH)
    name = (f"Fig3{panel.upper()}_simulation_lambda-histogram_{HYPOTHESIS[beta_type]}"
            f"_threshold-{THRESHOLD[tau]}.pdf")
    save_figure(hist, name, PATHS["figures"], PAPER["figure_size"])
