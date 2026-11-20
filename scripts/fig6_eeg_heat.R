# -----------------------------------------------------------------------------
# File:    scripts/fig6_eeg_heat.R
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
# Usage:   Rscript scripts/fig6_eeg_heat.R [--rerun] [--runs=1000] [--perms=1000]
#          [--negative=500] [--positive=2000] [--cores=4]
#          --rerun: recompute the curves from the data (ADAMANT_DATA).
# Source:  R_Heat2.R (heat1, tile1; Dec 2021 version);
#          EEG_t184.R, RPC_fcns.R (the curves).
# -----------------------------------------------------------------------------

source("R/setup.R")

# Triples and midpoints of the published figure (Dec 2021 version of R_Heat2.R).
TRIPLES <- c(181, 231, 200, 144, 159, 136, 223, 205, 25, 179, 219, 123)  # bottom to top
MIDPOINTS      <- c(logp = 0.75, power = 0.4)  # white points of the colour scales
LAMBDA_DIGITS  <- 3
AXIS_TEXT_SIZE <- 8

opts <- parse_args(curve_options(EEG_GRID, PAPER))

save_heat_maps("Fig6", "EEG", "_other-triples", eeg_curves(TRIPLES, opts, PAPER, PATHS), "power",
               MIDPOINTS, LAMBDA_DIGITS, opts$negative, AXIS_TEXT_SIZE,
               PATHS$figures, PAPER$figure_size)
