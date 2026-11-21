# -----------------------------------------------------------------------------
# File:    tests/testthat/helper-inputs.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Inputs shared by the unit tests (testthat loads this file before
#          the tests): small seeded matrices, temporary folders, small
#          options and settings. python/tests/conftest.py has the same.
# Inputs:  None.
# Outputs: None.
# Usage:   loaded by testthat; see tests/testthat.R
# -----------------------------------------------------------------------------

REPO <- normalizePath(file.path("..", ".."))

#' A matrix of standard normal numbers drawn with its own seed.
#'
#' @param n,p Numbers of rows and columns.
#' @param seed Seed of the draw.
#' @return n x p matrix, filled by column from `rnorm(n * p)`.
seeded_matrix <- function(n, p, seed) {
  set.seed(seed)
  matrix(rnorm(n * p), n, p)
}

#' A new empty temporary folder.
#'
#' @return Its path.
new_folder <- function() {
  folder <- tempfile("adamant-test-")
  dir.create(folder)
  folder
}

#' Run code without what it prints.
#'
#' @param code Expression that prints time-stamped messages.
#' @return The value of `code`.
quietly <- function(code) {
  capture.output(value <- code)
  value
}

# Settings as PAPER of R/setup.R with few permutations.
SMALL_SETTINGS <- list(seed = 8891, alpha = 0.05, sigma = 0.5, upper = 1000, n_perms = 9,
                       figure_size = 3)
# Options of a small rerun: 2 replications, 9 permutations, a 3 + 3 grid.
SMALL_OPTS <- list(rerun = TRUE, runs = 2L, perms = 9, negative = 3, positive = 3, cores = 1)

#' A curve table as `curve_table()` returns it, with made-up values.
#'
#' @param shift Number added to every curve, to tell tables apart.
#' @return Data frame with the columns `CURVE_COLUMNS` and 5 penalties.
made_up_curves <- function(shift = 0) {
  lambdas <- c(-1, -0.5, 0, 500, 1000)
  data.frame(lambdas = lambdas, ave_logpval = shift + c(1, 1.5, 2, 1.2, 1),
             ave_logpval_MOM = shift + c(1.1, 1.4, 2.1, 1.3, 0.9),
             power = c(0.2, 0.4, 0.8, 0.5, 0.3), power_MOM = c(0.25, 0.45, 0.75, 0.5, 0.3),
             gcv = shift + c(9, 4, 1, 2, 3), loo = shift + c(8, 5, 1.5, 2.5, 3.5))
}

#' Write a made-up file dat_ch<k>.csv in the layout of the EEG files: a row
#' number, the coherence, then the SNPs.
#'
#' @param eeg_dir Folder of the file.
#' @param k Triple number; also the seed of the values.
#' @param n,p Numbers of subjects and SNPs.
#' @return The values written, as a matrix with 2 + p columns.
write_eeg_file <- function(eeg_dir, k, n = 12, p = 5) {
  values <- cbind(seq_len(n), round(seeded_matrix(n, 1 + p, seed = k), 3))
  colnames(values) <- c("id", "coherence", paste0("snp", seq_len(p)))
  write.csv(values, file.path(eeg_dir, sprintf("dat_ch%d.csv", k)), row.names = FALSE)
  values
}

#' Write a made-up file in the layout of mixOmicsData.csv: a row label, the
#' 10 clinical variables, then the genes.
#'
#' @param file Path of the file.
#' @param n,p Numbers of animals and genes.
#' @return The values written, as a matrix with 10 + p columns.
write_liver_file <- function(file, n = 12, p = 4) {
  values <- round(seeded_matrix(n, 10 + p, seed = 41), 3)
  colnames(values) <- c(paste0("clinic", 1:10), paste0("gene", seq_len(p)))
  write.csv(data.frame(id = paste0("ID", seq_len(n)), values), file, row.names = FALSE)
  values
}
