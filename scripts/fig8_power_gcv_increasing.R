# -----------------------------------------------------------------------------
# File:    scripts/fig8_power_gcv_increasing.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Figure 8.
# Purpose: Power of the ridge-penalized Mantel test (top) and GCV score (bottom)
#          against the ridge penalty for three EEG triples whose power
#          increases with the penalty.
# Inputs:  archive/coh7/curves_eeg<k>.csv.gz for the three triples, reference/*.csv;
#          with --rerun <data>/eeg/dat_ch<k>.csv instead of the archive.
# Outputs: figures/fig8_{power,gcv}_<k>.pdf;
#          with --rerun also results/power_eeg<k>_<settings>.csv and
#          results/gcv_eeg<k>_<grid>.csv (shared by Figures 5 to 8).
# Usage:   Rscript scripts/fig8_power_gcv_increasing.R [--rerun] [--runs=1000] [--perms=1000]
#          [--negative=500] [--positive=2000] [--cores=4]
#          --rerun: recompute the curves from the data (ADAMANT_DATA).
# Source:  R_Output_plots.R (point1), GCV_plots.R / Output_plots.R (GCV);
#          EEG_t184.R, RPC_fcns.R (the curves).
# -----------------------------------------------------------------------------

source("R/setup.R")

# Rows of the curve table shown in each GCV panel and the largest GCV shown,
# as in the published panels.
GCV_WINDOWS <- data.frame(triple = c(191, 184, 160),
                          first  = c(413, 453, 31),
                          last   = c(650, 650, 650),
                          y_max  = c(16, 15.3, 21.6625))

opts <- parse_args(list(runs = 1000L, perms = N_PERMS, negative = EEG_GRID$negative,
                        positive = EEG_GRID$positive, cores = 4))

save_power_gcv_panels("fig8", GCV_WINDOWS, opts)
