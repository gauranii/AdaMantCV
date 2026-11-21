# -----------------------------------------------------------------------------
# File:    tests/testthat/test-setup.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Unit tests of R/setup.R: the settings, the project paths, the command line
#          parser, the index and label helpers and the csv cache.
# Inputs:  None (files are written to temporary folders).
# Outputs: None.
# Usage:   Rscript tests/testthat.R      # all test files
# -----------------------------------------------------------------------------

source("../../R/setup.R")

test_that("sourced_file is the file of the innermost source() call", {
  folder <- new_folder()
  writeLines("inner <- basename(sourced_file())", file.path(folder, "inner.R"))
  writeLines(c(sprintf("source('%s', local = TRUE)", file.path(folder, "inner.R")),
               "outer <- basename(sourced_file())"), file.path(folder, "outer.R"))
  seen <- new.env()
  source(file.path(folder, "outer.R"), local = seen)
  expect_equal(c(seen$inner, seen$outer), c("inner.R", "outer.R"))
  expect_error(sourced_file(), "must be loaded with source")
})

test_that("setup.R finds the repository and loads every module", {
  expect_true(file.exists(file.path(REPO_ROOT, "flake.nix")))
  for (name in c("adamant_cv", "ridge_cv", "simulation_study", "power_curve", "read_liver",
                 "eeg_tests", "cv_curves", "heat_map")) {
    expect_true(is.function(get(name)), info = name)
  }
  expect_equal(PAPER[c("seed", "alpha", "sigma", "upper", "n_perms")],
               list(seed = 8891, alpha = 0.05, sigma = 0.5, upper = 1000, n_perms = 1000))
})

test_that("project_paths puts data and outputs where it is told", {
  paths <- project_paths("/repo", data = "/elsewhere", output = "/out")
  expect_named(paths, c("eeg", "liver", "reference", "archive", "results", "figures"))
  expect_equal(paths$eeg, "/elsewhere/04_coherence_SNPs")
  expect_equal(paths$liver, "/repo/data/liver/mixOmicsData.csv")
  expect_equal(paths$reference, "/repo/reference")
  expect_equal(paths$archive, "/repo/archive")
  expect_equal(paths$results, "/out/results")
  expect_equal(paths$figures, "/out/figures")
  expect_true(file.exists(PATHS$liver))
})

test_that("output_file creates the folder and returns the path", {
  folder <- file.path(new_folder(), "a", "b")
  expect_equal(output_file(folder, "x.csv"), file.path(folder, "x.csv"))
  expect_true(dir.exists(folder))
})

test_that("parse_args converts values to the type of their default", {
  defaults <- list(runs = "500", perms = 1000, cores = 4L, fast = FALSE)
  expect_equal(parse_args(defaults, character(0)), c(list(rerun = FALSE), defaults))
  opts <- parse_args(defaults, c("--rerun", "--runs=1:100", "--perms=50", "--cores=2", "--fast"))
  expect_equal(opts, list(rerun = TRUE, runs = "1:100", perms = 50, cores = 2L, fast = TRUE))
  expect_type(opts$cores, "integer")
})

test_that("parse_args refuses what it does not know", {
  defaults <- list(perms = 1000, fast = FALSE)
  expect_error(parse_args(defaults, "--other=1"), "unknown argument --other")
  expect_error(parse_args(defaults, "perms=5"), "must look like")
  expect_error(parse_args(defaults, "--perms"), "--perms needs a value")
  expect_error(parse_args(defaults, "--fast=1"), "--fast takes no value")
})

test_that("curve_options and test_options are the defaults of the scripts", {
  expect_equal(curve_options(list(negative = 50L, positive = 250L), PAPER),
               list(runs = 1000L, perms = 1000, negative = 50L, positive = 250L, cores = 4))
  expect_equal(test_options(PAPER), list(perms = 1000, negative = 50, positive = 200, cores = 4))
  expect_equal(test_options(PAPER, cores = 5)$cores, 5)
})

test_that("parse_index_set reads ranges, lists and counts", {
  expect_equal(parse_index_set("3:7"), 3:7)
  expect_equal(parse_index_set("2,5,9"), c(2L, 5L, 9L))
  expect_equal(parse_index_set("4"), 4L)
  expect_equal(parse_index_set("4", count = TRUE), 1:4)
  expect_equal(parse_index_set("2,5", count = TRUE), c(2L, 5L))
})

test_that("settings_tag names the settings of a rerun", {
  opts <- list(perms = 1000, negative = 50, positive = 200)
  expect_equal(settings_tag(opts), "perms1000_grid50x200")
  expect_equal(settings_tag(opts, n_runs = 25L), "runs25_perms1000_grid50x200")
})

test_that("cached_csv computes once and then reads the file", {
  file <- file.path(new_folder(), "cache.csv")
  calls <- 0
  compute <- function() {
    calls <<- calls + 1
    data.frame(x = c(1, 2.5), label = c("a", "b"))
  }
  first <- cached_csv(file, compute)
  expect_true(file.exists(file))
  expect_output(second <- cached_csv(file, compute), "reusing")
  expect_equal(calls, 1)
  expect_equal(second, first)
})

test_that("log_message prints the time and the formatted message", {
  expect_output(log_message("%d of %s", 3L, "five"),
                "^[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9:]{8} - 3 of five")
  expect_null(quietly(log_message("x")))
})
