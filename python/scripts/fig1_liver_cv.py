# -----------------------------------------------------------------------------
# File:    python/scripts/fig1_liver_cv.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Figure 1.
# Purpose: Ten-fold cross-validated MSE of ridge regression against log(lambda)
#          for the ten clinical variables of the liver toxicity data, with 50
#          randomly chosen genes (A) and with all 3116 genes (B).
# Inputs:  data/liver/mixOmicsData.csv.
# Outputs: figures/Fig1A_liver_ridge-CV-error_p50-random-genes.pdf,
#          figures/Fig1B_liver_ridge-CV-error_p3116-all-genes.pdf,
#          results/fig1_minima.csv.
# Usage:   python python/scripts/fig1_liver_cv.py
# Source:  scripts/fig1_liver_cv.R.
# -----------------------------------------------------------------------------
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import math  # noqa: E402

import pandas as pd  # noqa: E402

from adamant import rbase, rng  # noqa: E402
from adamant.data import read_liver  # noqa: E402
from adamant.liver_cv import cv_curves, cv_minima, cv_plot, skip_fold_assignments  # noqa: E402
from adamant.plots import save_figure  # noqa: E402
from adamant.setup import PAPER, PATHS, output_file, parse_args, write_csv  # noqa: E402

SUBSET_SEED = 889191 + 1  # seed of panel A
SUBSET_SIZE = 50
SUBSET_LAMBDAS = [math.pow(10.0, e) for e in rbase.seq_by(4, -2.5, -0.001)]
# Panel B continues the random number stream of the original session: seed
# 889191 + 2, one sample of 50 genes, then eleven fold assignments.
FULL_SEED = 889191 + 2
FULL_SKIPPED_FOLDS = 11
FULL_LAMBDAS = [math.pow(10.0, e) for e in rbase.seq_by(9, -5, -0.1)]
N_FOLDS = 10
MSE_AXIS_MAX = 1.2

opts = parse_args({})  # no options; Figure 1 is always computed from the data

liver = read_liver(PATHS["liver"])
n_genes = liver["X"].shape[1]

rng.set_seed(SUBSET_SEED)
genes = rng.sample_int(n_genes, SUBSET_SIZE) - 1
subset_curves = cv_curves(liver["X"][:, genes], liver["Y"], SUBSET_LAMBDAS, N_FOLDS)

rng.set_seed(FULL_SEED)
rng.sample_int(n_genes, SUBSET_SIZE)
skip_fold_assignments(FULL_SKIPPED_FOLDS, liver["X"].shape[0], N_FOLDS)
full_curves = cv_curves(liver["X"], liver["Y"], FULL_LAMBDAS, N_FOLDS)

save_figure(cv_plot(subset_curves, MSE_AXIS_MAX), "Fig1A_liver_ridge-CV-error_p50-random-genes.pdf",
            PATHS["figures"], PAPER["figure_size"])
save_figure(cv_plot(full_curves, MSE_AXIS_MAX), "Fig1B_liver_ridge-CV-error_p3116-all-genes.pdf",
            PATHS["figures"], PAPER["figure_size"])
write_csv(pd.concat([cv_minima(subset_curves, "A"), cv_minima(full_curves, "B")],
                    ignore_index=True),
          output_file(PATHS["results"], "fig1_minima.csv"))
