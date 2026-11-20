# -----------------------------------------------------------------------------
# File:    python/scripts/fig6_eeg_heat.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Figure 6.
# Purpose: Heat maps of the mean -log10 p-value (A) and the power (B) of the
#          ridge-penalized Mantel test over the ridge penalties, for the 12
#          EEG triples with weak association.
# Inputs:  archive/coh7/curves_eeg<k>.csv.gz for the 12 triples, reference/*.csv;
#          with --rerun <data>/04_coherence_SNPs/dat_ch<k>.csv instead of the archive.
# Outputs: figures/Fig6A_EEG_heatmap_mean-log10-pvalue_other-triples.pdf,
#          figures/Fig6B_EEG_heatmap_power_other-triples.pdf;
#          with --rerun also results/power_eeg<k>_<settings>.csv and
#          results/gcv_eeg<k>_<grid>.csv (shared by Figures 5 to 8).
# Usage:   python python/scripts/fig6_eeg_heat.py [--rerun] [--runs=1000] [--perms=1000]
#          [--negative=500] [--positive=2000] [--cores=4]
#          --rerun: recompute the curves from the data (ADAMANT_DATA).
# Source:  scripts/fig6_eeg_heat.R.
# -----------------------------------------------------------------------------
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from adamant.data import EEG_GRID, eeg_curves  # noqa: E402
from adamant.plots import save_heat_maps  # noqa: E402
from adamant.setup import PAPER, PATHS, curve_options, parse_args  # noqa: E402

# Triples and midpoints of the published figure (Dec 2021 version of R_Heat2.R).
TRIPLES = [181, 231, 200, 144, 159, 136, 223, 205, 25, 179, 219, 123]  # bottom to top
MIDPOINTS = {"logp": 0.75, "power": 0.4}  # white points of the colour scales
LAMBDA_DIGITS = 3
AXIS_TEXT_SIZE = 8

opts = parse_args(curve_options(EEG_GRID, PAPER))

save_heat_maps("Fig6", "EEG", "_other-triples", eeg_curves(TRIPLES, opts, PAPER, PATHS), "power",
               MIDPOINTS, LAMBDA_DIGITS, opts["negative"], AXIS_TEXT_SIZE,
               PATHS["figures"], PAPER["figure_size"])
