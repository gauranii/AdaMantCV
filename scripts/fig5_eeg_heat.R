# -----------------------------------------------------------------------------
# File:    scripts/fig5_eeg_heat.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Figure 5.
# Purpose: Heat maps of the mean -log10 p-value (A) and the power (B) of the
#          ridge-penalized Mantel test over the ridge penalties, for the 12
#          EEG triples with the strongest association.
# Inputs:  archive/coh7/curves_eeg<k>.csv.gz for the 12 triples, reference/*.csv;
#          with --rerun <data>/04_coherence_SNPs/dat_ch<k>.csv instead of the archive.
# Outputs: figures/Fig5A_EEG_heatmap_mean-log10-pvalue_top-triples.pdf,
#          figures/Fig5B_EEG_heatmap_power_top-triples.pdf;
#          with --rerun also results/power_eeg<k>_<settings>.csv and
#          results/gcv_eeg<k>_<grid>.csv (shared by Figures 5 to 8).
# Usage:   Rscript scripts/fig5_eeg_heat.R [--rerun] [--runs=1000] [--perms=1000]
#          [--negative=500] [--positive=2000] [--cores=4]
#          --rerun: recompute the curves from the data (ADAMANT_DATA).
# Source:  R_Heat1.R (heat1, tile1); EEG_t184.R, RPC_fcns.R (the curves).
# -----------------------------------------------------------------------------

source("R/setup.R")

TRIPLES <- c(160, 248, 191, 183, 184, 220, 148, 177, 22, 157, 31, 151)  # bottom to top
MIDPOINTS      <- c(logp = 1.25, power = 0.4)  # white points of the colour scales
LAMBDA_DIGITS  <- 3
AXIS_TEXT_SIZE <- 10

opts <- parse_args(curve_options(EEG_GRID, PAPER))

save_heat_maps("Fig5", "EEG", "_top-triples", eeg_curves(TRIPLES, opts, PAPER, PATHS), "power",
               MIDPOINTS, LAMBDA_DIGITS, opts$negative, AXIS_TEXT_SIZE,
               PATHS$figures, PAPER$figure_size)
