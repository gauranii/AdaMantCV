# -----------------------------------------------------------------------------
# File:    scripts/table4_eeg_pvalues.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Table 4.
# Purpose: AdaMantCV p-values (GCV, LOOCV and gamma approximations) of the
#          association between SNPs and EEG coherence for five triples.
# Inputs:  <data>/eeg/dat_ch<k>.csv for k = 11, 27, 70, 178, 248; reference/*.csv.
# Outputs: results/table4_eeg_pvalues_<version>_<settings>.csv.
# Usage:   Rscript scripts/table4_eeg_pvalues.R [--version=best_fit] [--perms=1000]
#          [--negative=50] [--positive=200] [--cores=5]
#          --version: best_fit (default) the original code, with the observed
#                       data's lambda = 0 p-values replayed from the 2022
#                       cluster run (reference/table4_cluster_lambda0.csv);
#                     original  the original code, rounding of this machine;
#                     paper     the method as the paper states it (ties
#                       counted, CV on each permutation's own response).
# Source:  real0210a: R_ch11.R, R_driver.R; PctRk_<k>q_0.5.csv (cluster run).
# -----------------------------------------------------------------------------

source("R/setup.R")

TRIPLES  <- c(11, 27, 70, 178, 248)
QUANTILE <- 0.5  # tau = Q2

opts <- parse_args(list(version = "best_fit", perms = N_PERMS, negative = 50,
                        positive = 200, cores = 5))
stopifnot(opts$version %in% c("best_fit", "original", "paper"))
cluster <- read.csv(reference_path("table4_cluster_lambda0.csv"))

runs <- parallel::mclapply(TRIPLES, function(k) {
  data <- read_eeg_triple(k)
  replay <- if (opts$version == "best_fit") {
    unlist(cluster[cluster$triple == k, c("perms_below", "gamma")])
  }
  ada <- adamant_cv(data$X, data$y, QUANTILE, opts$negative, opts$positive,
                    opts$perms, SEED, UPPER,
                    method = if (opts$version == "paper") "paper" else "original",
                    replay_lambda0 = replay)
  c(ada$p_values, lambda_gcv = ada$lambda_gcv[1], lambda_loo = ada$lambda_loo[1],
    perms_at_zero_gcv = sum(ada$lambda_gcv[-1] == 0),
    perms_at_zero_loo = sum(ada$lambda_loo[-1] == 0))
}, mc.cores = opts$cores)

table <- cbind(triple_info(TRIPLES), do.call(rbind, runs))
write.csv(table, results_path(sprintf("table4_eeg_pvalues_%s_%s.csv", opts$version,
                                      settings_tag(opts))), row.names = FALSE)
print(table, row.names = FALSE)
