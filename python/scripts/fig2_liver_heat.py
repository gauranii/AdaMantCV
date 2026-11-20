# -----------------------------------------------------------------------------
# File:    python/scripts/fig2_liver_heat.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Figure 2.
# Purpose: Heat maps of the mean -log10 p-value (A) and the power (B) of the
#          ridge-penalized Mantel test over the ridge penalties, for the ten
#          clinical variables of the liver toxicity data (3116 genes).
# Inputs:  archive/liver/curves_liver<j>.csv.gz;
#          with --rerun data/liver/mixOmicsData.csv instead.
# Outputs: figures/Fig2A_liver_heatmap_mean-log10-pvalue.pdf,
#          figures/Fig2B_liver_heatmap_power.pdf;
#          with --rerun also results/power_liver<j>_<settings>.csv and
#          results/gcv_liver<j>_<grid>.csv.
# Usage:   python python/scripts/fig2_liver_heat.py [--rerun] [--runs=1000] [--perms=1000]
#          [--negative=50] [--positive=250] [--cores=4]
#          --rerun: recompute the curves from the data.
# Source:  scripts/fig2_liver_heat.R.
# -----------------------------------------------------------------------------
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from adamant.data import LIVER_GRID, liver_curves  # noqa: E402
from adamant.plots import save_heat_maps  # noqa: E402
from adamant.setup import PAPER, PATHS, curve_options, parse_args  # noqa: E402

MIDPOINTS = {"logp": 4, "power": 0.5}  # white points of the colour scales
LAMBDA_DIGITS = 1
AXIS_TEXT_SIZE = 8

opts = parse_args(curve_options(LIVER_GRID, PAPER))

# Panel B shows the power of the gamma approximation (column power_MOM), as
# drawn by R_Heat3.R.
save_heat_maps("Fig2", "liver", "", liver_curves(opts, PAPER, PATHS), "power_MOM",
               MIDPOINTS, LAMBDA_DIGITS, opts["negative"] - 1, AXIS_TEXT_SIZE,
               PATHS["figures"], PAPER["figure_size"])
