# -----------------------------------------------------------------------------
# File:    scripts/fig2_liver_heat.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Figure 2.
# Purpose: Heat maps of the mean -log10 p-value (A) and the power (B) of the
#          ridge-penalized Mantel test over the ridge penalties, for the ten
#          clinical variables of the liver toxicity data (3116 genes).
# Inputs:  archive/liver/curves_liver<j>.csv.gz;
#          with --rerun data/liver/mixOmicsData.csv instead.
# Outputs: figures/fig2{a_logp,b_power}_liver_heat.pdf;
#          with --rerun also results/power_liver<j>_<settings>.csv and
#          results/gcv_liver<j>_<grid>.csv.
# Usage:   Rscript scripts/fig2_liver_heat.R [--rerun] [--runs=1000] [--perms=1000]
#          [--negative=50] [--positive=250] [--cores=4]
#          --rerun: recompute the curves from the data.
# Source:  laptop_coh10: RPC_read.R, RPC_init.R, RPC_fcns.R, R_Heat3.R.
# -----------------------------------------------------------------------------

source("R/setup.R")

MIDPOINTS      <- c(logp = 4, power = 0.5)  # white points of the colour scales
LAMBDA_DIGITS  <- 1
AXIS_TEXT_SIZE <- 8

opts <- parse_args(curve_options(LIVER_GRID, PAPER))

# Panel B shows the power of the gamma approximation (column power_MOM), as
# drawn by R_Heat3.R.
save_heat_maps("fig2", "liver_heat.pdf", liver_curves(opts, PAPER, PATHS), "power_MOM",
               MIDPOINTS, LAMBDA_DIGITS, opts$negative - 1, AXIS_TEXT_SIZE,
               PATHS$figures, PAPER$figure_size)
