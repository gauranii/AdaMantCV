# -----------------------------------------------------------------------------
# File:    scripts/fig5_eeg_heat.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Figure 5.
# Purpose: Heat maps of the mean -log10 p-value (A) and the power (B) of the
#          ridge-penalized Mantel test over the ridge penalties, for the 12
#          EEG triples with the strongest association.
# Inputs:  <data>/eeg/dat_ch<k>.csv for the 12 triples; reference/*.csv.
# Outputs: results/power_eeg<k>_<settings>.csv, results/gcv_eeg<k>_<grid>.csv
#          (shared with Figs 5 to 8), figures/fig5{a_logp,b_power}_eeg_heat.pdf.
# Usage:   Rscript scripts/fig5_eeg_heat.R [--runs=1000] [--perms=1000] [--negative=500]
#          [--positive=2000] [--cores=4]
# Source:  R_Heat1.R (heat1, tile1); EEG_t184.R, RPC_fcns.R (the curves).
# -----------------------------------------------------------------------------

source("R/setup.R")

TRIPLES <- c(160, 248, 191, 183, 184, 220, 148, 177, 22, 157, 31, 151)  # bottom to top
MIDPOINT_LOGP  <- 1.25  # white point of the colour scale, panel A
MIDPOINT_POWER <- 0.4  # white point of the colour scale, panel B
LAMBDA_DIGITS  <- 3
AXIS_TEXT_SIZE <- 10

opts <- parse_args(list(runs = 1000L, perms = N_PERMS, negative = EEG_GRID$negative,
                        positive = EEG_GRID$positive, cores = 4))

curves <- eeg_curves(TRIPLES, opts)
save_figure(heat_map(curves, "ave_logpval", MIDPOINT_LOGP, LAMBDA_DIGITS,
                     opts$negative, AXIS_TEXT_SIZE), "fig5a_logp_eeg_heat.pdf")
save_figure(heat_map(curves, "power", MIDPOINT_POWER, LAMBDA_DIGITS,
                     opts$negative, AXIS_TEXT_SIZE), "fig5b_power_eeg_heat.pdf")
