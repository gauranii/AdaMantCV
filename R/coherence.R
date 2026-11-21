# -----------------------------------------------------------------------------
# File:    R/coherence.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Input of Table 4 and Figures 5 to 8.
# Purpose: How the EEG input of the paper was made, as functions: the
#          periodogram matrix of the channels of one trial, its averages
#          within the five frequency bands and over the trials, the band
#          coherence, the files per channel and band (coherences and the 497
#          candidate SNPs) and the files per channel-channel-band triple.
# Inputs:  None (the data are private; see coherence/README.md).
# Outputs: None.
# Usage:   source("R/coherence.R"); band_coherence(eeg, freq_cutoff = 0.4)
# Source:  coherence/coherence_functions.R (get_pgram, get_spectmat,
#          coherence), coherence/compute_band_coherence.R,
#          coherence/build_channel_band_files.R. The script that wrote the
#          files per triple was not found; write_triple_files() is
#          reconstructed from the files (checks/check_coherence.sh).
# -----------------------------------------------------------------------------

# Edges of the bands delta, theta, alpha, beta, gamma, in cycles per trial
# (Hz for trials of one second).
COHERENCE_BAND_EDGES <- c(0, 4, 8, 16, 32, 45)

# Rows of the table of candidate SNPs that are used: its 499 rows without the
# two SNPs (rows 164 and 357) that are not among the genotyped SNPs.
SNP_ROWS <- c(1:163, 165:356, 358:499)

# ----- Step 1: periodogram, band averages, coherence -----

#' Periodogram matrix of the channels of one trial.
#'
#' @param eeg_trial Matrix, channels x time points.
#' @return Complex array, time points x channels x channels: at Fourier
#'   frequency (k - 1) / time points, the outer product d_k d_k^H of the
#'   discrete Fourier transforms of the channels (not scaled).
get_pgram <- function(eeg_trial) {
  d <- t(mvfft(t(eeg_trial)))
  n_chans <- dim(eeg_trial)[1]
  n_freq <- dim(d)[2]
  I <- array(dim = c(n_freq, n_chans, n_chans))
  for (k in 1:n_freq) I[k, , ] <- d[, k] %*% Conj(t(d[, k]))
  I
}

#' Frequencies `spec.pgram()` reports for a series.
#'
#' The original code took its frequencies from `spec.pgram(x, plot = FALSE)$freq`.
#' Only the frequencies were used, so the taper and the detrending of
#' `spec.pgram()` play no part; its padding does: the series is padded to
#' `nextn(n_time)` points.
#'
#' @param n_time Number of time points.
#' @return The frequencies 1/N, 2/N, ..., floor(N/2)/N in cycles per time
#'   point, with N = `nextn(n_time)`.
pgram_frequencies <- function(n_time) {
  N <- nextn(n_time)
  seq.int(from = 1 / N, by = 1 / N, length.out = floor(N / 2))
}

#' Band of each frequency.
#'
#' @param freq Frequencies in cycles per time point.
#' @param n_time Number of time points of a trial.
#' @return Integer vector: band b (1 to 5) when
#'   edge[b] / n_time < freq <= edge[b + 1] / n_time, otherwise NA.
band_of_frequency <- function(freq, n_time) {
  as.integer(cut(freq, breaks = COHERENCE_BAND_EDGES / n_time))
}

#' Rows of the periodogram matrix that are averaged in each band.
#'
#' As in the original code: the rows of `get_pgram()` are first restricted
#' to those whose frequency is below `freq_cutoff`, the frequencies being
#' recycled along the rows, and the rows of band b are then those at the
#' positions of the frequencies of band b. Row i of `get_pgram()` is the
#' Fourier frequency (i - 1) / n_time while `pgram_frequencies()` starts at
#' 1 / N, so every band lies one row below its edges, and the delta band
#' includes frequency zero.
#'
#' @param n_time Number of time points of a trial.
#' @param freq_cutoff Frequencies at or above it are dropped, in cycles per
#'   time point (0.4 in compute_band_coherence.R, 50 / n_time in
#'   `get_spectmat()`).
#' @return List of five integer vectors of row numbers.
band_rows <- function(n_time, freq_cutoff) {
  freq <- pgram_frequencies(n_time)
  kept <- which(rep_len(freq < freq_cutoff, n_time))
  band <- band_of_frequency(freq, n_time)
  rows <- lapply(seq_len(length(COHERENCE_BAND_EDGES) - 1), function(b) kept[which(band == b)])
  stopifnot(!anyNA(unlist(rows)))
  rows
}

# ----- Step 2: band coherence of one subject -----

#' Spectral matrix of one subject in the five bands.
#'
#' @param eeg Array, channels x time points x trials.
#' @param freq_cutoff See `band_rows()`.
#' @return Complex array, 5 x channels x channels: the periodogram matrix
#'   averaged over the rows of each band, then over the trials.
band_spectrum <- function(eeg, freq_cutoff) {
  stopifnot(length(dim(eeg)) == 3)
  n_chans <- dim(eeg)[1]
  n_trials <- dim(eeg)[3]
  rows <- band_rows(dim(eeg)[2], freq_cutoff)
  trial_spectra <- array(dim = c(n_trials, length(rows), n_chans, n_chans))
  for (trial in 1:n_trials) {
    I <- get_pgram(eeg[, , trial])
    for (band in seq_along(rows)) {
      trial_spectra[trial, band, , ] <- apply(I[rows[[band]], , , drop = FALSE], c(2, 3), mean)
    }
  }
  apply(trial_spectra, c(2, 3, 4), mean)
}

#' Coherence from a spectral matrix.
#'
#' @param spectmat Complex array, bands x channels x channels.
#' @return Array of the same size: |S_jk|^2 / (|S_jj| |S_kk|) in each band.
coherence_matrix <- function(spectmat) {
  cohmat <- array(dim = dim(spectmat))
  for (band in seq_len(dim(spectmat)[1])) {
    S <- spectmat[band, , ]
    power <- Mod(diag(S))
    cohmat[band, , ] <- Mod(S)^2 / (power %o% power)
  }
  cohmat
}

#' Band coherence of one subject.
#'
#' @param eeg Array, channels x time points x trials.
#' @param freq_cutoff See `band_rows()`.
#' @return Array, 5 x channels x channels.
band_coherence <- function(eeg, freq_cutoff) {
  coherence_matrix(band_spectrum(eeg, freq_cutoff))
}

# ----- Step 3: one file per channel and band -----

#' Names of the candidate SNPs that are used.
#'
#' @param file Path of cups_top_500_SNPs.txt, one SNP per line. As in the
#'   original code its first line is read as a header, which leaves 499.
#' @return Character vector of the 497 names in the rows `SNP_ROWS`.
read_snp_names <- function(file) {
  as.character(read.delim(file)[SNP_ROWS, 1])
}

#' Genotypes of the candidate SNPs.
#'
#' @param genotypes Data frame, subjects x SNPs, with the SNP names as
#'   column names (`dat_snp_all`).
#' @param snp_names Names of the SNPs to take.
#' @return Matrix, subjects x SNPs in the order of `snp_names`, no names.
candidate_snps <- function(genotypes, snp_names) {
  stopifnot(snp_names %in% colnames(genotypes))
  unname(vapply(snp_names, function(name) unlist(genotypes[[name]]),
                numeric(nrow(genotypes))))
}

#' Coherences of one channel in one band, and the SNPs.
#'
#' @param coherences Array, subjects x channels x channels x bands
#'   (`dat_array`).
#' @param snps Matrix of `candidate_snps()`.
#' @param channel,band Channel and band numbers.
#' @return Data frame with columns V1, V2, ...: the coherence of every
#'   channel with `channel`, `coherences[, , channel, band]`, then the SNPs.
channel_band_table <- function(coherences, snps, channel, band) {
  as.data.frame(cbind(as.matrix(coherences[, , channel, band]), snps))
}

#' Write the files dat_ch<channel>b<band>.csv of one band.
#'
#' @param coherences,snps,band As for `channel_band_table()`.
#' @param out_dir Folder the files are written to.
#' @return None.
write_channel_band_files <- function(coherences, snps, band, out_dir) {
  for (channel in seq_len(dim(coherences)[3])) {
    write.csv(channel_band_table(coherences, snps, channel, band),
              file.path(out_dir, sprintf("dat_ch%db%d.csv", channel, band)))
  }
}

# ----- Step 4: one file per triple -----

#' Coherence of one channel pair and the SNPs, from a table of step 3.
#'
#' @param channel_band Table of `channel_band_table()` (or the file written
#'   from it, read with its row names) of one channel and band.
#' @param channel Number of the other channel.
#' @param n_chans Number of channels: the columns before the SNPs.
#' @return Data frame with columns V1, V2, ...: the coherence in column
#'   `channel` of the table, then its SNPs.
triple_table <- function(channel_band, channel, n_chans) {
  as.data.frame(unname(cbind(channel_band[, channel],
                             as.matrix(channel_band[, -seq_len(n_chans)]))))
}

#' Write the files dat_ch<k>.csv of a set of triples.
#'
#' Triple k (chan1, chan2, band) is cut from the file of step 3 of channel
#' chan2 and that band, column chan1: `coherences[, chan1, chan2, band]`.
#' (The coherence array is symmetric only up to rounding; this is the
#' orientation of the files of the paper.)
#'
#' @param triples Data frame with columns `chan1`, `chan2`, `band`; its row
#'   names are the triple numbers k (`triple.csv` read with `row.names = 1`).
#' @param n_chans Number of channels.
#' @param channel_band_dir Folder with the files of step 3.
#' @param out_dir Folder the files are written to.
#' @return None.
write_triple_files <- function(triples, n_chans, channel_band_dir, out_dir) {
  for (i in seq_len(nrow(triples))) {
    channel_band <- read.csv(file.path(channel_band_dir,
                                       sprintf("dat_ch%db%d.csv", triples$chan2[i], triples$band[i])),
                             row.names = 1)
    write.csv(triple_table(channel_band, triples$chan1[i], n_chans),
              file.path(out_dir, sprintf("dat_ch%s.csv", rownames(triples)[i])))
  }
}
