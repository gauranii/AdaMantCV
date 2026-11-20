# -----------------------------------------------------------------------------
# File:    scripts/fig8_power_gcv_increasing.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Figure 8.
# Purpose: Power of the ridge-penalized Mantel test (top) and GCV score (bottom)
#          against the ridge penalty for three EEG triples whose power
#          increases with the penalty.
# Inputs:  archive/coh7/curves_eeg<k>.csv.gz for the three triples, reference/*.csv;
#          with --rerun <data>/04_coherence_SNPs/dat_ch<k>.csv instead of the archive.
# Outputs: figures/Fig8{A,B,C}_EEG_power-increasing_<channel>-<channel>-<band>.pdf,
#          figures/Fig8{D,E,F}_EEG_GCV_<channel>-<channel>-<band>.pdf;
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

opts <- parse_args(curve_options(EEG_GRID, PAPER))

save_power_gcv_panels("Fig8", "increasing", GCV_WINDOWS, eeg_curves(GCV_WINDOWS$triple, opts, PAPER, PATHS),
                      PATHS$figures, PAPER$figure_size)
