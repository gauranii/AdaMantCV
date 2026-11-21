# -----------------------------------------------------------------------------
# File:    tests/testthat/test-data.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Unit tests of R/data.R: the readers on small made-up files and on the liver
#          data, the labels of the EEG triples, and the curve tables from the
#          archived runs and from a small rerun.
# Inputs:  data/liver/mixOmicsData.csv, reference/, archive/coh7/, archive/liver/;
#          made-up EEG files are written to temporary folders (the EEG data
#          are not used).
# Outputs: None.
# Usage:   Rscript tests/testthat.R      # all test files
# -----------------------------------------------------------------------------

source("../../R/setup.R")

test_that("read_eeg_triple standardizes the coherence and the SNPs", {
  eeg_dir <- new_folder()
  values <- write_eeg_file(eeg_dir, k = 7)
  triple <- read_eeg_triple(7, eeg_dir)
  expect_equal(dim(triple$X), c(12, 5))
  expect_length(triple$y, 12)
  expect_equal(triple$y, (values[, 2] - mean(values[, 2])) / sd(values[, 2]))
  expect_equal(unname(triple$X[, 3]), (values[, 5] - mean(values[, 5])) / sd(values[, 5]))
  expect_equal(unname(colMeans(triple$X)), rep(0, 5))
  expect_equal(unname(apply(triple$X, 2, sd)), rep(1, 5))
})

test_that("read_liver splits a file into standardized clinical variables and genes", {
  file <- file.path(new_folder(), "liver.csv")
  values <- write_liver_file(file)
  liver <- read_liver(file)
  expect_equal(dim(liver$Y), c(12, 10))
  expect_equal(dim(liver$X), c(12, 4))
  expect_equal(unname(liver$Y[, 10]), (values[, 10] - mean(values[, 10])) / sd(values[, 10]))
  expect_equal(unname(liver$X[, 1]), (values[, 11] - mean(values[, 11])) / sd(values[, 11]))
})

test_that("read_liver reads the liver toxicity data", {
  liver <- read_liver(PATHS$liver)
  expect_equal(dim(liver$X), c(64, 3116))
  expect_equal(dim(liver$Y), c(64, 10))
  expect_equal(unname(colMeans(liver$Y)), rep(0, 10))
  expect_equal(unname(apply(liver$X[, 1:5], 2, sd)), rep(1, 5))
  expect_false(anyNA(liver$X))
})

test_that("triple_info names the channels and the band", {
  info <- triple_info(c(160, 1), PATHS$reference)
  expect_named(info, c("triple", "channel1", "channel2", "band"))
  expect_equal(info$triple, c(160, 1))
  expect_equal(unlist(info[1, 2:4], use.names = FALSE), c("Pz", "PO4", "alpha"))
  expect_equal(unlist(info[2, 2:4], use.names = FALSE), c("FC3", "CP2", "delta"))
})

test_that("archived_curve_file points into the archive", {
  expect_equal(archived_curve_file("eeg", 184, "/archive"), "/archive/coh7/curves_eeg184.csv.gz")
  expect_equal(archived_curve_file("liver", 3, "/archive"), "/archive/liver/curves_liver3.csv.gz")
  expect_true(file.exists(archived_curve_file("eeg", 160, PATHS$archive)))
})

test_that("eeg_curves reads the archived curve tables", {
  curves <- eeg_curves(c(160, 184), list(rerun = FALSE), PAPER, PATHS)
  expect_named(curves, c("Pz-PO4-alpha", "F3-O2-alpha"))
  for (table in curves) {
    expect_named(table, CURVE_COLUMNS)
    expect_equal(nrow(table), EEG_GRID$negative + EEG_GRID$positive)
    expect_false(is.unsorted(table$lambdas))
    expect_true(all(table$power >= 0 & table$power <= 1))
  }
})

test_that("eeg_curves recomputes the curves of a made-up triple", {
  paths <- project_paths(REPO_ROOT, data = new_folder(), output = new_folder())
  dir.create(paths$eeg)
  write_eeg_file(paths$eeg, k = 160)
  curves <- eeg_curves(160, SMALL_OPTS, SMALL_SETTINGS, paths)
  expect_named(curves, "Pz-PO4-alpha")
  triple <- read_eeg_triple(160, paths$eeg)
  expected <- quietly(curve_table("eeg160", triple$X, triple$y, SMALL_OPTS, SMALL_SETTINGS,
                                  new_folder()))
  expect_equal(curves[[1]], expected)
  expect_true(file.exists(file.path(paths$results, "gcv_eeg160_grid3x3.csv")))
})

test_that("liver_curves reads the archived curve tables of the ten variables", {
  curves <- liver_curves(list(rerun = FALSE), PAPER, PATHS)
  expect_named(curves, LIVER_LABELS)
  for (table in curves) {
    expect_named(table, CURVE_COLUMNS)
    expect_equal(nrow(table), LIVER_GRID$negative + LIVER_GRID$positive)
  }
})

test_that("liver_curves recomputes the curves from a made-up file", {
  paths <- project_paths(REPO_ROOT, output = new_folder())
  paths$liver <- file.path(new_folder(), "liver.csv")
  write_liver_file(paths$liver)
  curves <- liver_curves(SMALL_OPTS, SMALL_SETTINGS, paths)
  expect_named(curves, LIVER_LABELS)
  liver <- read_liver(paths$liver)
  expected <- quietly(curve_table("liver10", liver$X, liver$Y[, 10], SMALL_OPTS, SMALL_SETTINGS,
                                  new_folder()))
  expect_equal(curves$Cholesterol, expected)
  expect_length(list.files(paths$results), 20)
})
