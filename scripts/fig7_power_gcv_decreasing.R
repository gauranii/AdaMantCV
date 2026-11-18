# -----------------------------------------------------------------------------
# File:    scripts/fig7_power_gcv_decreasing.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Figure 7.
# Purpose: Power of the ridge-penalized Mantel test (top) and GCV score (bottom)
#          against the ridge penalty for three EEG triples whose power
#          decreases with the penalty.
# Inputs:  archive/coh7/curves_eeg<k>.csv.gz for the three triples, reference/*.csv;
#          with --rerun <data>/eeg/dat_ch<k>.csv instead of the archive.
# Outputs: figures/fig7_{power,gcv}_<k>.pdf;
#          with --rerun also results/power_eeg<k>_<settings>.csv and
#          results/gcv_eeg<k>_<grid>.csv (shared by Figures 5 to 8).
# Usage:   Rscript scripts/fig7_power_gcv_decreasing.R [--rerun] [--runs=1000] [--perms=1000]
#          [--negative=500] [--positive=2000] [--cores=4]
#          --rerun: recompute the curves from the data (ADAMANT_DATA).
# Source:  R_Output_plots.R (point1), GCV_plots.R / Output_plots.R (GCV);
#          EEG_t184.R, RPC_fcns.R (the curves).
# -----------------------------------------------------------------------------

source("R/setup.R")

# Rows of the curve table shown in each GCV panel and the largest GCV shown,
# as in the published panels.
GCV_WINDOWS <- data.frame(triple = c(223, 22, 25),
                          first  = c(59, 329, 79),
                          last   = c(750, 650, 750),
                          y_max  = c(25, 25, 18))

opts <- parse_args(curve_options(EEG_GRID, PAPER))

save_power_gcv_panels("fig7", GCV_WINDOWS, eeg_curves(GCV_WINDOWS$triple, opts, PAPER, PATHS),
                      PATHS$figures, PAPER$figure_size)
