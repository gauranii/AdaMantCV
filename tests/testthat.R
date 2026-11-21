# -----------------------------------------------------------------------------
# File:    tests/testthat.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. All figures and tables.
# Purpose: Run the unit tests of the functions in R/: every file
#          tests/testthat/test-*.R, on small synthetic inputs.
# Inputs:  tests/testthat/test-*.R; the tests read reference/, archive/ and
#          data/liver/ and write only to temporary folders.
# Outputs: prints one line per test file; exit status 1 if a test fails.
# Usage:   Rscript tests/testthat.R      # from anywhere, inside `nix develop`
#          nix run .#test                # this and the Python tests
# -----------------------------------------------------------------------------

library(testthat)

runner <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
tests_dir <- file.path(dirname(normalizePath(runner)), "testthat")

test_dir(tests_dir, reporter = "summary", stop_on_failure = TRUE, stop_on_warning = FALSE)
