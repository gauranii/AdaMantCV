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
# Usage:   source("R/setup.R"); read_eeg_triple(184, PATHS$eeg)
# Source:  EEG_t184.R, RPC_read.R, real0210a/R_driver.R.
# -----------------------------------------------------------------------------

EEG_BANDS <- c("delta", "theta", "alpha", "beta", "gamma")
LIVER_LABELS <- c("BUN", "Creatinine", "Total Protein", "Albumin", "ALT", "SDH",
                  "AST", "ALP", "TBA", "Cholesterol")
EEG_GRID   <- list(negative = 500L, positive = 2000L)  # Dec 2021 EEG runs
LIVER_GRID <- list(negative = 50L, positive = 250L)    # Dec 2021 liver runs

#' EEG coherence and SNPs of one triple (channel pair and frequency band).
#'
#' @param k Triple number.
#' @param eeg_dir Folder with the files dat_ch<k>.csv.
#' @return List with `X` (standardized SNPs, 497 columns) and `y`
#'   (standardized coherence).
read_eeg_triple <- function(k, eeg_dir) {
  data <- as.matrix(read.csv(file.path(eeg_dir, sprintf("dat_ch%d.csv", k))))
  list(X = scale(data[, 3:ncol(data)]), y = scale(data[, 2])[, 1])
}

#' Liver toxicity data.
#'
#' @param file Path of mixOmicsData.csv.
#' @return List with `X` (standardized gene expressions, 3116 columns) and
#'   `Y` (standardized clinical variables, 10 columns).
read_liver <- function(file) {
  data <- read.csv(file)
  list(X = scale(as.matrix(data[, 12:ncol(data)])),
       Y = scale(as.matrix(data[, 2:11])))
}

#' Channel names and band of EEG triples.
#'
#' @param triples Triple numbers.
#' @param reference_dir Folder with triple.csv and chan_names.csv.
#' @return Data frame with columns `triple`, `channel1`, `channel2`, `band`
#'   (for example Pz, PO4, alpha for triple 160).
triple_info <- function(triples, reference_dir) {
  triple <- read.csv(file.path(reference_dir, "triple.csv"), row.names = 1)[triples, ]
  channels <- trimws(read.csv(file.path(reference_dir, "chan_names.csv"), header = FALSE)[[1]])
  data.frame(triple = triples, channel1 = channels[triple$chan1],
             channel2 = channels[triple$chan2], band = EEG_BANDS[triple$band])
}

#' File of the curve table of the archived run of one response.
#'
#' @param kind "eeg" or "liver".
#' @param k Triple number (eeg) or number of the clinical variable (liver).
#' @param archive_dir The archive folder.
#' @return Path of the csv.gz file.
archived_curve_file <- function(kind, k, archive_dir) {
  file.path(archive_dir, if (kind == "eeg") "coh7" else "liver",
            sprintf("curves_%s%d.csv.gz", kind, k))
}

#' Curve tables of EEG triples, from the archive or recomputed.
#'
#' @param triples Triple numbers.
#' @param opts Options: `rerun` and those of `curve_table()`.
#' @param settings Settings of the paper.
#' @param paths Project paths.
#' @return List of curve tables named "<channel1>-<channel2>-<band>".
eeg_curves <- function(triples, opts, settings, paths) {
  info <- triple_info(triples, paths$reference)
  curves <- lapply(triples, function(k) {
    if (!opts$rerun) return(read.csv(archived_curve_file("eeg", k, paths$archive)))
    data <- read_eeg_triple(k, paths$eeg)
    curve_table(sprintf("eeg%d", k), data$X, data$y, opts, settings, paths$results)
  })
  setNames(curves, paste(info$channel1, info$channel2, info$band, sep = "-"))
}

#' Curve tables of the ten liver clinical variables, from the archive or
#' recomputed.
#'
#' @param opts,settings,paths As for `eeg_curves()`.
#' @return List of curve tables named by clinical variable.
liver_curves <- function(opts, settings, paths) {
  if (opts$rerun) data <- read_liver(paths$liver)
  curves <- lapply(seq_along(LIVER_LABELS), function(j) {
    if (!opts$rerun) return(read.csv(archived_curve_file("liver", j, paths$archive)))
    curve_table(sprintf("liver%d", j), data$X, data$Y[, j], opts, settings, paths$results)
  })
  setNames(curves, LIVER_LABELS)
}
