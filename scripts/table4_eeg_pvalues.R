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
#          with --rerun <data>/eeg/dat_ch<k>.csv for k = 11, 27, 70, 178, 248
#          and reference/table4_cluster_lambda0.csv (observed test statistics
#          saved by the original run) instead of the archive.
# Outputs: results/table4_eeg_pvalues.csv.
# Usage:   Rscript scripts/table4_eeg_pvalues.R [--rerun] [--perms=1000]
#          [--negative=50] [--positive=200] [--cores=5]
#          --rerun: recompute the tests from the data (ADAMANT_DATA).
# Source:  real0210a: R_ch11.R, R_driver.R; PctRk_<k>q_0.5.csv (original run).
# -----------------------------------------------------------------------------

source("R/setup.R")

QUANTILE <- 0.5  # tau = Q2

opts <- parse_args(list(perms = N_PERMS, negative = 50, positive = 200, cores = 5))

# Adaptive p-values and chosen penalties of the observed data, one row per
# triple, from the archived p-values at the chosen penalties.
from_archive <- function() {
  archived <- read.csv(archive_path("real0210a", "table4_pvalues.csv"))
  t(vapply(TABLE4_TRIPLES, function(k) {
    chosen <- archived[archived$triple == k, ]
    c(gcv = adaptive_pvalue(chosen$p_gcv), loo = adaptive_pvalue(chosen$p_loo),
      ga_gcv = adaptive_pvalue(chosen$p_ga_gcv), ga_loo = adaptive_pvalue(chosen$p_ga_loo),
      lambda_gcv = chosen$lambda_gcv[1], lambda_loo = chosen$lambda_loo[1])
  }, numeric(6)))
}

# The same, recomputed from the data.
rerun <- function() {
  saved <- read.csv(reference_path("table4_cluster_lambda0.csv"))
  runs <- parallel::mclapply(TABLE4_TRIPLES, function(k) {
    data <- read_eeg_triple(k)
    ada <- adamant_cv(data$X, data$y, QUANTILE, opts$negative, opts$positive,
                      opts$perms, SEED, UPPER,
                      observed_lambda0 = unlist(saved[saved$triple == k,
                                                      c("perms_below", "gamma")]))
    c(ada$p_values, lambda_gcv = ada$lambda_gcv[1], lambda_loo = ada$lambda_loo[1])
  }, mc.cores = opts$cores)
  do.call(rbind, runs)
}

table <- cbind(triple_info(TABLE4_TRIPLES), if (opts$rerun) rerun() else from_archive())
write.csv(table, results_path("table4_eeg_pvalues.csv"), row.names = FALSE)
print(table, row.names = FALSE)
