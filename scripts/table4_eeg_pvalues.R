# -----------------------------------------------------------------------------
# File:    scripts/table4_eeg_pvalues.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Table 4.
# Purpose: AdaMantCV p-values (GCV, LOOCV and gamma approximations) of the
#          association between SNPs and EEG coherence for five triples.
# Inputs:  <data>/eeg/dat_ch<k>.csv for k = 11, 27, 70, 178, 248; reference/*.csv.
# Outputs: results/table4_eeg_pvalues_<settings>.csv.
# Usage:   Rscript scripts/table4_eeg_pvalues.R [--perms=1000] [--negative=50]
#          [--positive=200] [--cores=5]
# Source:  real0210a: R_ch11.R, R_driver.R.
# -----------------------------------------------------------------------------

source("R/setup.R")

TRIPLES  <- c(11, 27, 70, 178, 248)
QUANTILE <- 0.5  # tau = Q2

opts <- parse_args(list(perms = N_PERMS, negative = 50, positive = 200, cores = 5))

pvalues <- parallel::mclapply(TRIPLES, function(k) {
  data <- read_eeg_triple(k)
  adamant_cv(data$X, data$y, QUANTILE, opts$negative, opts$positive, opts$perms,
             SEED, UPPER)$p_values
}, mc.cores = opts$cores)

table <- cbind(triple_info(TRIPLES), do.call(rbind, pvalues))
write.csv(table, results_path(sprintf("table4_eeg_pvalues_%s.csv", settings_tag(opts))),
          row.names = FALSE)
print(table, row.names = FALSE)
