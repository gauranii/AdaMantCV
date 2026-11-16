# -----------------------------------------------------------------------------
# File:    scripts/fig2_liver_heat.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Figure 2.
# Purpose: Heat maps of the mean -log10 p-value (A) and the power (B) of the
#          ridge-penalized Mantel test over the ridge penalties, for the ten
#          clinical variables of the liver toxicity data (3116 genes).
# Inputs:  <data>/liver/mixOmicsData.csv.
# Outputs: results/power_liver<j>_<settings>.csv, results/gcv_liver<j>_<grid>.csv,
#          figures/fig2{a_logp,b_power}_liver_heat.pdf.
# Usage:   Rscript scripts/fig2_liver_heat.R [--runs=1000] [--perms=1000] [--negative=50]
#          [--positive=250] [--cores=4]
# Source:  laptop_coh10: RPC_read.R, RPC_init.R, RPC_fcns.R, R_Heat3.R.
# -----------------------------------------------------------------------------

source("R/setup.R")

MIDPOINT_LOGP  <- 4    # white point of the colour scale, panel A
MIDPOINT_POWER <- 0.5  # white point of the colour scale, panel B
LAMBDA_DIGITS  <- 1
AXIS_TEXT_SIZE <- 8

opts <- parse_args(list(runs = 1000L, perms = N_PERMS, negative = LIVER_GRID$negative,
                        positive = LIVER_GRID$positive, cores = 4))

curves <- liver_curves(opts)
save_figure(heat_map(curves, "ave_logpval", MIDPOINT_LOGP, LAMBDA_DIGITS,
                     opts$negative - 1, AXIS_TEXT_SIZE), "fig2a_logp_liver_heat.pdf")
# Panel B shows the power of the gamma approximation (power_MOM): in R_Heat3.R
# the file LiverHeatPower<p>.pdf was drawn from that column, while the
# permutation power went to LiverMHeatPower<p>.pdf.
save_figure(heat_map(curves, "power_MOM", MIDPOINT_POWER, LAMBDA_DIGITS,
                     opts$negative - 1, AXIS_TEXT_SIZE), "fig2b_power_liver_heat.pdf")
