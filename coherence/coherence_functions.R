## Libraries ####
library(tidyverse)
library(stringr)

## Functions ####

get_pgram <- function(eeg_trial) {
  d <- t(mvfft(t(eeg_trial)))
  n_chans <- dim(eeg_trial)[1]
  n_freq <- dim(d)[2]
  I <- array(dim=c(n_freq, n_chans, n_chans))
  for (k in 1:n_freq)
    I[k, , ] <- d[, k]%*%Conj(t(d[, k]))
  return(I)
}

get_pgram_builtin <- function(eeg, chan, trials, limits=c(), bandwidth=c(0,0.5), span_width=3) {
  if (length(limits) == 0)
    limits <- 1:length(eeg[1, , 1])
  freq <- spec.pgram(eeg[chan, limits, 1], spans=span_width, fast=F, plot=F)$freq
  pgram <- matrix(nrow=length(trials), ncol=length(freq))
  for (trial in trials)
    pgram[trial, ] <- log(spec.pgram(eeg[chan, limits, trial], spans=span_width, fast=F, plot=F)$spec) - digamma(1)
  pgram <- pgram[, freq > bandwidth[1] & freq < bandwidth[2]]
  freq <- freq[freq > bandwidth[1] & freq < bandwidth[2]]
  return(list(pgram=pgram, freq=freq))
}

## GET_SPECTMAT ####
# Calculates spectral matrix "manually" from fft of channels.
# eeg is an 3-dim array of Channels x Time pts x Trials
#   the data for one subject

get_spectmat <- function(eeg, freq_cutoff = 50) {
  if (length(dim(eeg)) != 3) {
    warning("Incorrect dimensions for EEG Object.  3-D array required.")
    return(NaN)
  }
  ## SET UP VARS
  freq <- spec.pgram(eeg[1, , 1], plot=F)$freq
  n_chans <- dim(eeg)[1]
  n_time_points <- dim(eeg)[2]
  n_freq <- dim(eeg)[2]/2
  n_trials <- dim(eeg)[3]
  n_bands <- 5
  freq_cutoff <- freq_cutoff / n_time_points
  
  freq_bands <- na.omit(cut(freq, breaks=(c(0, 4, 8, 16, 32, 45)/n_time_points)))
  
  trial_spectmats <- array(dim=c(n_trials, n_bands, n_chans, n_chans))
  for (trial in 1:n_trials) {
    I <- get_pgram(eeg[, , trial])[freq < freq_cutoff, , ]
    
    for (band in 1:n_bands) {
      inds <- which(freq_bands == levels(freq_bands)[band])
      trial_spectmats[trial, band, , ] <- apply(I[inds, , ], c(2, 3), mean)
    }
  }
  spectmat <- apply(trial_spectmats, c(2, 3, 4), mean)
  return(spectmat)
}

coherence <- function(eeg, freq_cutoff = 50) {
  spectmat <- get_spectmat(eeg, freq_cutoff)
  q <- dim(spectmat)[2]
  cohmat <- array(dim = c(5, q, q))
  for (band in 1:5) {
    for (j in 1:q) {
      for (k in 1:q)
        cohmat[band, j, k] <- Mod(spectmat[band, j, k])^2/(Mod(spectmat[band, j, j])*Mod(spectmat[band, k, k]))
    }
    cohmat[upper.tri(cohmat)] <- cohmat[lower.tri(cohmat)]
  }
  return(cohmat)
}

plot_eeg <- function(eeg, chan, trials, grid=FALSE, plot_mean=FALSE) {
  if (plot_mean)
    plot(colMeans(eeg[chan, , ]), ty="l")
  else {
    plot(eeg[chan, , trials[1]], ty="l",
         ylim=c(min(eeg[chan, , trials]), max(eeg[chan, , trials])))
    for (trial in trials) {
      lines(eeg[chan, , trial], col=rainbow(length(trials))[trial])
    }
  }
}

plot_pgram <- function(pgram, trials, bandwidth=c(0, 0.5), ylim=NA,
                       grid=FALSE, plot_mean=FALSE) {
  if (plot_mean)
    plot(pgram$freq, colMeans(pgram$pgram), xlim=bandwidth, ty="l")
  else {
    if(is.na(ylim))
      ylim <- c(min(pgram$pgram), max(pgram$pgram))
    plot(pgram$freq, pgram$pgram[trials[1], ], xlim=bandwidth,
         ylim=ylim, ty="l")
    for (trial in trials[2:length(trials)])
      lines(pgram$freq, pgram$pgram[trial, ], col=rainbow(length(trials))[trial])
  }
}

