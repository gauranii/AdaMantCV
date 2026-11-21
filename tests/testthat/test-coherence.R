# -----------------------------------------------------------------------------
# File:    tests/testthat/test-coherence.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Unit tests of R/coherence.R: agreement with the original functions
#          of coherence/coherence_functions.R and with the loops of
#          coherence/compute_band_coherence.R and build_channel_band_files.R
#          on made-up signals, the bands of the frequencies, the properties
#          of the coherence, and the files of steps 3 and 4.
# Inputs:  coherence/coherence_functions.R (its functions only; the libraries
#          it loads are not needed); made-up signals and files in temporary
#          folders (the EEG and SNP data are not used).
# Outputs: None.
# Usage:   Rscript tests/testthat.R      # all test files
# -----------------------------------------------------------------------------

source("../../R/coherence.R")

# The original functions, without the library() calls of their file.
original <- new.env()
for (expr in parse("../../coherence/coherence_functions.R")) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-"))) eval(expr, original)
}

#' A made-up recording: white noise in every channel, and one sinusoid shared
#' by channels 1 and 2 with a phase that changes from trial to trial.
#'
#' @param n_chans,n_time,n_trials Size of the recording.
#' @param cycles Cycles per trial of the shared sinusoid.
#' @param seed Seed of the draw.
#' @return Array, channels x time points x trials.
made_up_eeg <- function(n_chans, n_time, n_trials, cycles, seed) {
  set.seed(seed)
  eeg <- array(rnorm(n_chans * n_time * n_trials), c(n_chans, n_time, n_trials))
  for (trial in 1:n_trials) {
    wave <- 30 * sin(2 * pi * cycles * (1:n_time) / n_time + runif(1, 0, 2 * pi))
    eeg[1, , trial] <- eeg[1, , trial] + wave
    eeg[2, , trial] <- eeg[2, , trial] + wave
  }
  eeg
}

#' The loop of coherence/compute_band_coherence.R up to the spectral matrix.
#'
#' @param eeg Array, channels x time points x trials.
#' @param freq_cutoff Its `freq_cutoff`.
#' @return The 5 x channels x channels array `spectmat` of the script.
original_script_spectmat <- function(eeg, freq_cutoff) {
  get_pgram <- original$get_pgram
  freq <- spec.pgram(eeg[1, , 1], plot = FALSE)$freq
  n_chans <- dim(eeg)[1]
  n_time_points <- dim(eeg)[2]
  n_trials <- dim(eeg)[3]
  n_bands <- 5
  freq_bands <- na.omit(cut(freq, breaks = (c(0, 4, 8, 16, 32, 45) / n_time_points)))
  trial_spectmats <- array(dim = c(n_trials, n_bands, n_chans, n_chans))
  for (trial in 1:n_trials) {
    I <- get_pgram(eeg[, , trial])[freq < freq_cutoff, , ]
    for (band in 1:n_bands) {
      inds <- which(freq_bands == levels(freq_bands)[band])
      trial_spectmats[trial, band, , ] <- apply(I[inds, , ], c(2, 3), mean)
    }
  }
  apply(trial_spectmats, c(2, 3, 4), mean)
}

test_that("get_pgram is the original periodogram matrix", {
  trial <- made_up_eeg(4, 60, 1, cycles = 5, seed = 1)[, , 1]
  pgram <- get_pgram(trial)
  expect_identical(pgram, original$get_pgram(trial))
  expect_equal(dim(pgram), c(60, 4, 4))
  expect_equal(pgram[7, 2, 3], fft(trial[2, ])[7] * Conj(fft(trial[3, ])[7]))
  expect_equal(pgram[7, , ], Conj(t(pgram[7, , ])))
})

test_that("pgram_frequencies are those of spec.pgram, padded to nextn", {
  for (n_time in c(60, 127, 200, 256, 1000)) {
    expect_identical(pgram_frequencies(n_time), spec.pgram(rnorm(n_time), plot = FALSE)$freq)
  }
  expect_length(pgram_frequencies(127), 64)
})

test_that("frequencies are assigned to the five bands by their right-closed edges", {
  n_time <- 256
  cycles <- c(0, 0.5, 4, 4.5, 8, 9, 16, 17, 32, 33, 45, 46)
  expect_identical(band_of_frequency(cycles / n_time, n_time),
                   c(NA, 1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 5L, NA))
  expect_identical(tabulate(band_of_frequency(pgram_frequencies(n_time), n_time)),
                   c(4L, 4L, 8L, 16L, 13L))
})

test_that("band_rows lie one row below the band edges, the delta band from frequency zero", {
  rows <- band_rows(256, freq_cutoff = 0.4)
  expect_identical(rows, list(1:4, 5:8, 9:16, 17:32, 33:45))
  expect_identical(band_rows(256, freq_cutoff = 50 / 256), rows)
  # a series padded from 101 to 108 points has 54 frequencies, the gamma band
  # the 35th to the 48th; those below the cutoff are the first 43, which
  # recycled along the 101 rows keep rows 1 to 43 and 55 to 97
  expect_identical(band_rows(101, freq_cutoff = 0.4)[[5]], c(35:43, 55:59))
  # a cutoff inside the bands moves the later rows, as the original subsetting
  # does: rows 1 to 40 and, the frequencies being recycled, 129 to 168 are kept
  expect_identical(band_rows(256, freq_cutoff = 40.5 / 256)[[5]], c(33:40, 129:133))
  expect_error(band_rows(256, freq_cutoff = 3.5 / 256))
})

test_that("band_spectrum is the spectral matrix of the original function and script", {
  for (n_time in c(256, 200, 127, 101)) {
    eeg <- made_up_eeg(3, n_time, 4, cycles = 10, seed = n_time)
    expect_identical(band_spectrum(eeg, 50 / n_time), original$get_spectmat(eeg, freq_cutoff = 50))
    expect_identical(band_spectrum(eeg, 0.4), original_script_spectmat(eeg, 0.4))
  }
  eeg <- made_up_eeg(3, 256, 2, cycles = 10, seed = 7)
  expect_identical(band_spectrum(eeg, 40.5 / 256), original$get_spectmat(eeg, freq_cutoff = 40.5))
})

test_that("band_coherence is the original coherence", {
  eeg <- made_up_eeg(5, 256, 6, cycles = 10, seed = 2)
  expect_identical(band_coherence(eeg, 50 / 256), original$coherence(eeg, freq_cutoff = 50))
  expect_identical(band_coherence(eeg, 0.4),
                   coherence_matrix(original_script_spectmat(eeg, 0.4)))
})

test_that("channels sharing a sinusoid are coherent in its band only", {
  # 10 cycles per trial is row 11 of the periodogram: the alpha band (rows 9 to 16)
  coh <- band_coherence(made_up_eeg(5, 256, 40, cycles = 10, seed = 3), 0.4)
  expect_gt(coh[3, 1, 2], 0.99)
  expect_lt(max(coh[-3, 1, 2]), 0.2)
  expect_lt(max(coh[3, 1:2, 3:5]), 0.2)
  expect_lt(max(coh[3, 3:5, 3:5][upper.tri(diag(3))]), 0.2)
})

test_that("the coherence is symmetric, between 0 and 1, and 1 on the diagonal", {
  coh <- band_coherence(made_up_eeg(6, 200, 5, cycles = 20, seed = 4), 0.4)
  expect_equal(dim(coh), c(5, 6, 6))
  expect_true(all(coh >= 0 & coh <= 1 + 1e-12))
  for (band in 1:5) {
    expect_equal(coh[band, , ], t(coh[band, , ]))
    expect_equal(diag(coh[band, , ]), rep(1, 6))
  }
})

test_that("one trial has coherence 1 where a band is a single frequency", {
  # a single outer product d d^H has |S_jk|^2 = S_jj S_kk
  S <- get_pgram(made_up_eeg(4, 64, 1, cycles = 3, seed = 5)[, , 1])[3:4, , ]
  expect_equal(coherence_matrix(S), array(1, c(2, 4, 4)))
})

#' A made-up coherence array, SNP table and SNP list for steps 3 and 4.
#'
#' @param folder Folder the list of SNP names is written to.
#' @return List with `coherences` (6 subjects x 4 x 4 channels x 5 bands),
#'   `genotypes` (data frame of 499 SNPs without two of the listed ones),
#'   `file` (the list of 500 names).
made_up_inputs <- function(folder) {
  set.seed(6)
  names <- sprintf("rs%d", 1000 + 1:500)
  file <- file.path(folder, "snps.txt")
  writeLines(names, file)
  genotyped <- sample(names[-c(165, 358)])
  genotypes <- as.data.frame(matrix(rnorm(6 * 498), 6, dimnames = list(NULL, genotyped)))
  list(coherences = array(runif(6 * 4 * 4 * 5), c(6, 4, 4, 5)), genotypes = genotypes,
       file = file, names = names)
}

test_that("the 497 candidate SNPs are those of the original loop", {
  inputs <- made_up_inputs(new_folder())
  snp_names <- read_snp_names(inputs$file)
  expect_identical(snp_names, inputs$names[-c(1, 165, 358)])
  snps <- candidate_snps(inputs$genotypes, snp_names)
  # build_channel_band_files.R
  dat_snp_all <- inputs$genotypes
  X_sub <- read.delim(inputs$file)
  idx <- c(1:163, 165:356, 358:499)
  X_snp <- matrix(NA, nrow = 6, ncol = length(idx))
  for (i in seq_along(idx)) {
    j <- idx[i]
    X_snp[, i] <- unlist(dat_snp_all[, which(colnames(dat_snp_all) == X_sub[j, 1])])
  }
  expect_identical(snps, X_snp)
  expect_error(candidate_snps(inputs$genotypes, inputs$names[165]))
})

test_that("the files of step 3 are those of the original loop", {
  folder <- new_folder()
  inputs <- made_up_inputs(folder)
  dat_array <- inputs$coherences
  X_snp <- candidate_snps(inputs$genotypes, read_snp_names(inputs$file))
  write_channel_band_files(dat_array, X_snp, band = 5, folder)
  # build_channel_band_files.R, written to another folder
  expected <- new_folder()
  band <- 5
  for (ch in 1:4) {
    Y_coh <- as.matrix(dat_array[, , ch, band])
    df <- as.data.frame(cbind(Y_coh, X_snp))
    write.csv(df, file.path(expected, paste("dat_ch", ch, "b", band, ".csv", sep = "")))
  }
  for (ch in 1:4) {
    name <- sprintf("dat_ch%db5.csv", ch)
    expect_identical(readLines(file.path(folder, name)), readLines(file.path(expected, name)))
  }
  table <- read.csv(file.path(folder, "dat_ch3b5.csv"), row.names = 1)
  expect_equal(dim(table), c(6, 4 + 497))
  expect_equal(table$V2, dat_array[, 2, 3, 5])
})

test_that("the file of a triple is a row id, its coherence, then the SNPs", {
  step3 <- new_folder()
  step4 <- new_folder()
  inputs <- made_up_inputs(step3)
  snps <- candidate_snps(inputs$genotypes, read_snp_names(inputs$file))
  for (band in 1:5) write_channel_band_files(inputs$coherences, snps, band, step3)
  triples <- data.frame(chan1 = c(1, 2), chan2 = c(4, 3), band = c(2, 5), row.names = c(7, 9))
  write_triple_files(triples, n_chans = 4, step3, step4)
  expect_setequal(list.files(step4), c("dat_ch7.csv", "dat_ch9.csv"))
  lines <- readLines(file.path(step4, "dat_ch9.csv"))
  expect_identical(lines[1], paste(c('""', sprintf('"V%d"', 1:498)), collapse = ","))
  table <- read.csv(file.path(step4, "dat_ch9.csv"))
  expect_equal(dim(table), c(6, 499))
  expect_identical(table[[1]], 1:6)
  expect_equal(table[[2]], inputs$coherences[, 2, 3, 5])
  expect_equal(unname(as.matrix(table[, 3:499])), snps)
  # the same file written from the table in memory
  direct <- file.path(step4, "direct.csv")
  write.csv(triple_table(channel_band_table(inputs$coherences, snps, 3, 5), 2, n_chans = 4), direct)
  expect_identical(readLines(direct), lines)
})
