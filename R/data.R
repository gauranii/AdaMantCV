# -----------------------------------------------------------------------------
# File:    R/data.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Figures 1, 2, 5 to 8, Table 4.
# Purpose: Read the EEG coherence / SNP data and the liver toxicity data,
#          label the EEG triples, and get their power and cross-validation
#          curves, from the archived runs or recomputed.
# Inputs:  <data>/eeg/dat_ch<k>.csv, data/liver/mixOmicsData.csv,
#          reference/triple.csv, reference/chan_names.csv,
#          archive/coh7/curves_eeg<k>.csv.gz, archive/liver/curves_liver<j>.csv.gz.
# Outputs: None.
# Usage:   source("R/setup.R"); read_eeg_triple(184)
# Source:  EEG_t184.R, RPC_read.R, real0210a/R_driver.R.
# -----------------------------------------------------------------------------

EEG_BANDS <- c("delta", "theta", "alpha", "beta", "gamma")
LIVER_LABELS <- c("BUN", "Creatinine", "Total Protein", "Albumin", "ALT", "SDH",
                  "AST", "ALP", "TBA", "Cholesterol")
EEG_GRID   <- list(negative = 500L, positive = 2000L)  # Dec 2021 EEG runs
LIVER_GRID <- list(negative = 50L, positive = 250L)    # Dec 2021 liver runs

TABLE4_TRIPLES <- c(11, 27, 70, 178, 248)  # the five triples of Table 4

# Columns of a curve table: the penalties, the mean -log10 p-value and the
# power of the permutation test and of its gamma approximation (MOM), and the
# GCV and LOOCV scores.
CURVE_COLUMNS <- c("lambdas", "ave_logpval", "ave_logpval_MOM", "power", "power_MOM",
                   "gcv", "loo")

# Standardized SNPs (X, 497 columns) and coherence (y) of triple k
# (channel pair and frequency band).
read_eeg_triple <- function(k) {
  data <- as.matrix(read.csv(data_path("eeg", sprintf("dat_ch%d.csv", k))))
  list(X = scale(data[, 3:ncol(data)]), y = scale(data[, 2])[, 1])
}

# Standardized gene expressions (X, 3116 columns) and clinical variables (Y,
# 10 columns) of the liver toxicity data.
read_liver <- function() {
  data <- read.csv(file.path(REPO_ROOT, "data", "liver", "mixOmicsData.csv"))
  list(X = scale(as.matrix(data[, 12:ncol(data)])),
       Y = scale(as.matrix(data[, 2:11])))
}

# Channel names and band of each triple, e.g. Pz, PO4, alpha for triple 160.
triple_info <- function(triples) {
  triple <- read.csv(reference_path("triple.csv"), row.names = 1)[triples, ]
  channels <- trimws(read.csv(reference_path("chan_names.csv"), header = FALSE)[[1]])
  data.frame(triple = triples, channel1 = channels[triple$chan1],
             channel2 = channels[triple$chan2], band = EEG_BANDS[triple$band])
}

# Curve table of the archived run of one response: kind "eeg" (triple k) or
# "liver" (clinical variable k).
archived_curve_file <- function(kind, k) {
  archive_path(if (kind == "eeg") "coh7" else "liver", sprintf("curves_%s%d.csv.gz", kind, k))
}

# Curve tables of EEG triples, named "<channel1>-<channel2>-<band>".
eeg_curves <- function(triples, opts) {
  info <- triple_info(triples)
  curves <- lapply(triples, function(k) {
    if (!opts$rerun) return(read.csv(archived_curve_file("eeg", k)))
    data <- read_eeg_triple(k)
    curve_table(sprintf("eeg%d", k), data$X, data$y, opts)
  })
  setNames(curves, paste(info$channel1, info$channel2, info$band, sep = "-"))
}

# Curve tables of the ten liver clinical variables, named by variable.
liver_curves <- function(opts) {
  if (opts$rerun) data <- read_liver()
  curves <- lapply(seq_along(LIVER_LABELS), function(j) {
    if (!opts$rerun) return(read.csv(archived_curve_file("liver", j)))
    curve_table(sprintf("liver%d", j), data$X, data$Y[, j], opts)
  })
  setNames(curves, LIVER_LABELS)
}
