# -----------------------------------------------------------------------------
# File:    R/eeg_tests.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Table 4.
# Purpose: AdaMantCV tests of the association between SNPs and EEG coherence
#          for a set of triples, from the archived run or recomputed.
# Inputs:  archive/real0210a/table4_pvalues.csv; in a rerun
#          <data>/04_coherence_SNPs/dat_ch<k>.csv and reference/table4_cluster_lambda0.csv.
# Outputs: None.
# Usage:   source("R/setup.R"); eeg_tests(c(11, 27), 0.5, opts, PAPER, PATHS)
# Source:  real0210a: R_ch11.R, R_driver.R; PctRk_<k>q_0.5.csv (original run).
# -----------------------------------------------------------------------------

#' Tests of EEG triples from the archived p-values at the chosen penalties.
#'
#' @param triples Triple numbers.
#' @param file Path of table4_pvalues.csv: per triple, the p-values at the
#'   chosen penalties (`p_gcv`, `p_loo`, `p_ga_gcv`, `p_ga_loo`) and the
#'   chosen penalties of the observed data (b = 1) and every permutation.
#' @return Matrix with one row per triple: the four adaptive p-values (named
#'   as `TESTS`) and the penalties chosen for the observed data.
eeg_tests_archived <- function(triples, file) {
  archived <- read.csv(file)
  t(vapply(triples, function(k) {
    chosen <- archived[archived$triple == k, ]
    c(vapply(setNames(paste0("p_", TESTS), TESTS),
             function(column) adaptive_pvalue(chosen[[column]]), numeric(1)),
      lambda_gcv = chosen$lambda_gcv[1], lambda_loo = chosen$lambda_loo[1])
  }, numeric(length(TESTS) + 2)))
}

#' Tests of EEG triples recomputed from the data.
#'
#' The p-values of the observed data at lambda = 0 are taken from the values
#' saved by the original run (argument `observed_lambda0` of `adamant_cv()`).
#'
#' @param triples Triple numbers.
#' @param q Quantile below which singular values are set to zero.
#' @param opts Options: `perms`, `negative`, `positive`, `cores` (parallel
#'   processes over the triples).
#' @param settings Settings of the paper (`seed`, `upper`).
#' @param paths Project paths (`eeg`, `reference`).
#' @return Matrix as returned by `eeg_tests_archived()`.
eeg_tests_rerun <- function(triples, q, opts, settings, paths) {
  saved <- read.csv(file.path(paths$reference, "table4_cluster_lambda0.csv"))
  runs <- parallel::mclapply(triples, function(k) {
    data <- read_eeg_triple(k, paths$eeg)
    ada <- adamant_cv(data$X, data$y, q, opts$negative, opts$positive, opts$perms,
                      settings$seed, settings$upper,
                      observed_lambda0 = unlist(saved[saved$triple == k,
                                                      c("perms_below", "gamma")]))
    c(ada$p_values, lambda_gcv = ada$lambda_gcv[1], lambda_loo = ada$lambda_loo[1])
  }, mc.cores = opts$cores)
  do.call(rbind, runs)
}

#' Table of the tests of EEG triples, from the archive or recomputed.
#'
#' @param triples,q,opts,settings,paths As for `eeg_tests_rerun()`;
#'   `opts$rerun` chooses between the two.
#' @return Data frame with the labels of the triples (`triple_info()`), the
#'   four adaptive p-values and the penalties chosen for the observed data.
eeg_tests <- function(triples, q, opts, settings, paths) {
  tests <- if (opts$rerun) {
    eeg_tests_rerun(triples, q, opts, settings, paths)
  } else {
    eeg_tests_archived(triples, file.path(paths$archive, "real0210a", "table4_pvalues.csv"))
  }
  cbind(triple_info(triples, paths$reference), tests)
}
