
eeg <- dat_subj
freq <- spec.pgram(eeg[1, , 1], plot=F)$freq
n_chans <- dim(eeg)[1]
n_time_points <- dim(eeg)[2]
n_freq <- dim(eeg)[2]/2
n_trials <- dim(eeg)[3]
n_bands <- 5
freq_cutoff <- 0.4

freq_bands <- na.omit(cut(freq, breaks=(c(0, 4, 8, 16, 32, 45)/n_time_points)))

trial_spectmats <- array(dim=c(n_trials, n_bands, n_chans, n_chans))
for (trial in 1:n_trials) {
  I <- get_pgram(eeg[, , trial])[freq < freq_cutoff, , ]
  
  for (band in 1:n_bands) {
    inds <- which(freq_bands == levels(freq_bands)[band])
    trial_spectmats[trial, band, , ] <- apply(I[inds, , ], c(2, 3), mean)
  }
}
dim(trial_spectmats)
spectmat <- apply(trial_spectmats, c(2, 3, 4), mean)
dim(spectmat)
spectmat[1, , ]
coherence(spectmat[5, , ])
coherence_array <- apply(spectmat, 1, coherence)
