rm(list=ls())
setwd("~/Downloads/FRONTIERS/Coherence")
load("~/Downloads/FRONTIERS/Coherence/vwm_coh_clean.RData")
load("~/Downloads/FRONTIERS/Coherence/eeg_dat_snp_all.RData")
#source("RPC_init.R")
#source("RPC_fcns.R")

X_sub <- read.delim("cups_top_500_SNPs.txt")

idx <- c(1:163,165:356,358:499)
X_snp <- matrix(NA, nrow = 350, ncol = length(idx))
for (i in seq_along(idx)) {
  j <- idx[i]
  X_snp[,i] <- unlist(dat_snp_all[,which(colnames(dat_snp_all) == X_sub[j,1])])
}

band <- 5
for(ch in 1:64) {
  Y_coh <- as.matrix(dat_array[,,ch,band])
  df <- as.data.frame(cbind(Y_coh, X_snp))
  
  write.csv(df, paste("dat_ch", ch, "b", band,".csv", sep = ""))
}

