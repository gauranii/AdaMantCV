# -----------------------------------------------------------------------------
# File:    scripts/table4_eeg_pvalues.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Table 4.
# Purpose: AdaMantCV p-values (GCV, LOOCV and gamma approximations) of the
#          association between SNPs and EEG coherence for five triples.
# Inputs:  archive/real0210a/table4_pvalues.csv (p-values at the chosen
#          penalties of the observed data and every permutation, archived
#          run), reference/*.csv;
#          with --rerun <data>/04_coherence_SNPs/dat_ch<k>.csv for k = 11, 27, 70, 178, 248
#          and reference/table4_cluster_lambda0.csv (observed test statistics
#          saved by the original run) instead of the archive.
# Outputs: results/Table4_EEG_pvalues.csv.
# Usage:   Rscript scripts/table4_eeg_pvalues.R [--rerun] [--perms=1000]
#          [--negative=50] [--positive=200] [--cores=5]
#          --rerun: recompute the tests from the data (ADAMANT_DATA).
# Source:  real0210a: R_ch11.R, R_driver.R; PctRk_<k>q_0.5.csv (original run).
# -----------------------------------------------------------------------------

source("R/setup.R")

TRIPLES  <- c(11, 27, 70, 178, 248)
QUANTILE <- 0.5  # tau = Q2

opts <- parse_args(test_options(PAPER, cores = 5))

table <- eeg_tests(TRIPLES, QUANTILE, opts, PAPER, PATHS)
write.csv(table, output_file(PATHS$results, "Table4_EEG_pvalues.csv"), row.names = FALSE)
print(table, row.names = FALSE)
