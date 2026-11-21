# -----------------------------------------------------------------------------
# File:    tests/testthat/test-eeg_tests.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Unit tests of R/eeg_tests.R: the adaptive p-values from a hand-made file of
#          p-values at the chosen penalties, from the archived run, and from
#          a rerun on made-up EEG files.
# Inputs:  reference/, archive/real0210a/table4_pvalues.csv; made-up EEG files
#          are written to temporary folders (the EEG data are not used).
# Outputs: None.
# Usage:   Rscript tests/testthat.R      # all test files
# -----------------------------------------------------------------------------

source("../../R/setup.R")

test_that("eeg_tests_archived ranks the observed p-value among the permutations", {
  file <- file.path(new_folder(), "pvalues.csv")
  write.csv(data.frame(triple = rep(c(3, 8), each = 4), b = rep(1:4, 2),
                       p_gcv = c(0.2, 0.1, 0.2, 0.4, 0.01, 0.5, 0.6, 0.7),
                       p_loo = c(0.9, 0.1, 0.2, 0.4, 0.5, 0.5, 0.5, 0.5),
                       p_ga_gcv = c(0.3, 0.1, 0.2, 0.4, 0.6, 0.5, 0.7, 0.1),
                       p_ga_loo = c(0.1, 0.1, 0.2, 0.4, 0.2, 0.3, 0.1, 0.4),
                       lambda_gcv = c(0, 1000, 5, 5, -0.5, 1, 1, 1),
                       lambda_loo = c(1000, 0, 5, 5, 2, 1, 1, 1)), file, row.names = FALSE)
  tests <- eeg_tests_archived(c(8, 3), file)
  expect_equal(colnames(tests), c(TESTS, "lambda_gcv", "lambda_loo"))
  expect_equal(tests[1, ], c(gcv = 0, loo = 0, ga_gcv = 2 / 3, ga_loo = 1 / 3,
                             lambda_gcv = -0.5, lambda_loo = 2))
  expect_equal(tests[2, ], c(gcv = 1 / 3, loo = 1, ga_gcv = 2 / 3, ga_loo = 0,
                             lambda_gcv = 0, lambda_loo = 1000))
})

test_that("eeg_tests builds the table of the archived run", {
  table <- eeg_tests(c(11, 248), q = 0.5, list(rerun = FALSE), PAPER, PATHS)
  expect_named(table, c("triple", "channel1", "channel2", "band", TESTS,
                        "lambda_gcv", "lambda_loo"))
  expect_equal(table$triple, c(11, 248))
  expect_equal(table$channel1, c("Fp1", "P4"))
  expect_true(all(unlist(table[TESTS]) >= 0 & unlist(table[TESTS]) <= 1))
  # the observed p-value 0.015 of triple 11 is ranked among its 1000 permutations
  archived <- read.csv(file.path(PATHS$archive, "real0210a", "table4_pvalues.csv"))
  expect_equal(table$gcv[1], mean(archived$p_gcv[archived$triple == 11][-1] < 0.015))
})

test_that("eeg_tests_rerun tests made-up triples with the saved values at lambda = 0", {
  paths <- project_paths(REPO_ROOT, data = new_folder())
  dir.create(paths$eeg)
  for (k in c(11, 27)) write_eeg_file(paths$eeg, k)
  tests <- eeg_tests_rerun(c(11, 27), q = 0.5, SMALL_OPTS, SMALL_SETTINGS, paths)
  expect_equal(dim(tests), c(2, 6))
  expect_equal(colnames(tests), c(TESTS, "lambda_gcv", "lambda_loo"))

  triple <- read_eeg_triple(27, paths$eeg)
  ada <- adamant_cv(triple$X, triple$y, 0.5, 3, 3, 9, 8891, 1000,
                    observed_lambda0 = c(perms_below = 992, gamma = 0.0118191573598137))
  expect_equal(tests[2, ], c(ada$p_values, lambda_gcv = ada$lambda_gcv[1],
                             lambda_loo = ada$lambda_loo[1]))
  expect_equal(ada$p_gamma[1, ada$lambdas == 0], 0.0118191573598137)

  table <- eeg_tests(c(11, 27), 0.5, SMALL_OPTS, SMALL_SETTINGS, paths)
  expect_equal(unname(as.matrix(table[5:10])), unname(tests))
  expect_equal(table$band, c("delta", "delta"))
})
