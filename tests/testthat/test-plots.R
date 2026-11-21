# -----------------------------------------------------------------------------
# File:    tests/testthat/test-plots.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Unit tests of R/plots.R: every plot function returns a ggplot whose layers
#          hold the expected marks, and the save functions write the files
#          they name.
# Inputs:  None (figures are written to temporary folders).
# Outputs: None.
# Usage:   Rscript tests/testthat.R      # all test files
# -----------------------------------------------------------------------------

source("../../R/setup.R")

CURVES <- list(first = made_up_curves(), second = made_up_curves(shift = 1))

test_that("save_figure writes a square pdf", {
  figures <- file.path(new_folder(), "figures")
  plot <- sign_scatter(made_up_curves(), "power", "power")
  expect_output(save_figure(plot, "scatter.pdf", figures, size = 3), "wrote figures/scatter.pdf")
  file <- file.path(figures, "scatter.pdf")
  expect_gt(file.size(file), 1000)
  expect_equal(readBin(file, "raw", 4), charToRaw("%PDF"))
})

test_that("heat_map has one tile per table and penalty and the line at zero", {
  plot <- heat_map(CURVES, "power", midpoint = 0.4, digits = 3, zero_position = 3,
                   text_size = 10)
  expect_s3_class(plot, "ggplot")
  layers <- ggplot_build(plot)$data
  tiles <- layers[[1]]
  expect_equal(nrow(tiles), 10)
  expect_equal(as.vector(table(tiles$y)), c(5, 5))     # one row of tiles per table
  expect_equal(sort(unique(as.numeric(tiles$x))), 1:5)
  expect_equal(as.numeric(layers[[2]]$xintercept), 3)
  # the midpoint is white, above it red, below it blue
  rgb <- grDevices::col2rgb(tiles$fill)
  expect_true(all(rgb["red", CURVES$first$power > 0.4] > rgb["blue", CURVES$first$power > 0.4]))
  expect_true(all(rgb["red", CURVES$first$power < 0.4] < rgb["blue", CURVES$first$power < 0.4]))
  expect_equal(plot$labels$x, "lambda")
})

test_that("heat_map rounds the penalties, so close penalties share a column", {
  curves <- list(only = transform(made_up_curves(), lambdas = c(-1, -0.5004, -0.5001, 0, 1000)))
  tiles <- ggplot_build(heat_map(curves, "gcv", 2, digits = 2, 3, 10))$data[[1]]
  expect_equal(sort(unique(as.numeric(tiles$x))), 1:4)
})

test_that("save_heat_maps writes panels A and B with their names", {
  figures <- new_folder()
  quietly(save_heat_maps("Fig9", "EEG", "_test", CURVES, "power_MOM", c(logp = 1.5, power = 0.4),
                         digits = 3, zero_position = 3, text_size = 10, figures, size = 3))
  expect_equal(sort(list.files(figures)),
               c("Fig9A_EEG_heatmap_mean-log10-pvalue_test.pdf",
                 "Fig9B_EEG_heatmap_power_test.pdf"))
})

test_that("sign_scatter has one point per penalty, coloured by its sign", {
  table <- made_up_curves()
  plot <- sign_scatter(table, "gcv", "Generalized Cross Validation")
  expect_s3_class(plot, "ggplot")
  points <- ggplot_build(plot)$data[[1]]
  expect_equal(points$x, table$lambdas)
  expect_equal(points$y, table$gcv)
  expect_equal(length(unique(points$colour)), 3)
  expect_equal(match(points$colour, unique(points$colour)), c(1, 1, 2, 3, 3))
  expect_equal(plot$labels$y, "Generalized Cross Validation")
})

test_that("save_power_gcv_panels writes a power and a GCV panel per triple", {
  figures <- new_folder()
  windows <- data.frame(triple = c(1, 2), first = c(2, 1), last = c(4, 99), y_max = c(3, 100))
  quietly(save_power_gcv_panels("Fig9", "rising", windows, CURVES, figures, size = 3))
  expect_equal(sort(list.files(figures)),
               c("Fig9A_EEG_power-rising_first.pdf", "Fig9B_EEG_power-rising_second.pdf",
                 "Fig9C_EEG_GCV_first.pdf", "Fig9D_EEG_GCV_second.pdf"))
})

test_that("penalty_histogram counts every chosen penalty once per method", {
  chosen <- data.frame(GCV = c(0.1, 0.9, 1.1, 2.0, 2.0, 5.2), LOOCV = c(0.1, 0.1, 0.1, 5, 5, 5))
  plot <- penalty_histogram(chosen, tau = "Q1", bin_width = 1)
  expect_s3_class(plot, "ggplot")
  layers <- ggplot_build(plot)$data
  expect_length(layers, 2)
  # bins of width 1 centred on the integers
  gcv <- layers[[1]][layers[[1]]$count > 0, ]
  expect_equal(gcv$x, c(0, 1, 2, 5))
  expect_equal(gcv$count, c(1, 2, 2, 1))
  loocv <- layers[[2]][layers[[2]]$count > 0, ]
  expect_equal(loocv$x, c(0, 5))
  expect_equal(loocv$count, c(3, 3))
  expect_equal(levels(factor(plot$data$variable)), c("GCV tau = Q1", "LOOCV tau = Q1"))
})
