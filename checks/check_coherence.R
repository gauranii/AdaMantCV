# -----------------------------------------------------------------------------
# File:    checks/check_coherence.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Rebuild steps 3 and 4 of the coherence pipeline with R/coherence.R
#          from the private data and compare the files of the triples with
#          those the paper used; hand the inputs to checks/check_coherence.py
#          as binary files (Python has no reader of .RData files).
# Inputs:  <data>/01_raw/eeg_dat_snp_all.RData, cups_top_500_SNPs.txt,
#          <data>/02_coherence/vwm_coh_clean.RData,
#          <data>/03_triples_channels/triple.csv,
#          <data>/04_coherence_SNPs/dat_ch<k>.csv.
# Outputs: In the work folder: genotypes.bin and genotyped.txt (the genotypes
#          of the listed SNPs and their names), coherences.bin (the coherence
#          array), r/step4/dat_ch<k>.csv, r/step3.md5 (checksums of the
#          files of step 3, which are deleted band by band) and r/summary.csv.
#          Prints counts and the largest difference, never data.
# Usage:   Rscript checks/check_coherence.R <data folder> <work folder>
#          (called by checks/check_coherence.sh from the repository root)
# -----------------------------------------------------------------------------

source("R/coherence.R")
args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]
work <- args[2]
snp_file <- file.path(data_dir, "01_raw", "cups_top_500_SNPs.txt")

load(file.path(data_dir, "01_raw", "eeg_dat_snp_all.RData"))  # dat_snp_all
snps <- candidate_snps(dat_snp_all, read_snp_names(snp_file))
genotyped <- intersect(readLines(snp_file), colnames(dat_snp_all))
writeBin(as.vector(candidate_snps(dat_snp_all, genotyped)), file.path(work, "genotypes.bin"))
writeLines(genotyped, file.path(work, "genotyped.txt"))
rm(dat_snp_all)
invisible(gc())

load(file.path(data_dir, "02_coherence", "vwm_coh_clean.RData"))  # dat_array
writeBin(as.vector(dat_array), file.path(work, "coherences.bin"))
n_chans <- dim(dat_array)[2]
triples <- read.csv(file.path(data_dir, "03_triples_channels", "triple.csv"), row.names = 1)

step3 <- file.path(work, "r", "step3")
step4 <- file.path(work, "r", "step4")
dir.create(step3, recursive = TRUE)
dir.create(step4)
checksums <- character(0)
for (band in seq_len(dim(dat_array)[4])) {
  write_channel_band_files(dat_array, snps, band, step3)
  checksums <- c(checksums, tools::md5sum(list.files(step3, full.names = TRUE)))
  write_triple_files(triples[triples$band == band, ], n_chans, step3, step4)
  unlink(list.files(step3, full.names = TRUE))
}
writeLines(paste(checksums, basename(names(checksums))), file.path(work, "r", "step3.md5"))

#' Compare the rebuilt file of one triple with the file of the paper.
#'
#' @param k Triple number.
#' @return Named vector: whether the two files have the same bytes, whether
#'   the file of the paper has the bytes of the rebuilt numbers written
#'   correctly rounded to 15 significant digits, the number of fields that
#'   then differ, and the largest absolute difference of the numbers.
compare_triple <- function(k) {
  name <- sprintf("dat_ch%d.csv", k)
  paper <- readLines(file.path(data_dir, "04_coherence_SNPs", name))
  rebuilt <- triple_table(channel_band_table(dat_array, snps, triples$chan2[k], triples$band[k]),
                          triples$chan1[k], n_chans)
  rounded <- matrix(sprintf("%.15g", as.matrix(rebuilt)), nrow(rebuilt))
  fields <- do.call(rbind, strsplit(paper[-1], ",", fixed = TRUE))
  stopifnot(identical(dim(fields), dim(rounded) + c(0L, 1L)),
            identical(paper[1], readLines(file.path(step4, name), n = 1)),
            identical(fields[, 1], sprintf('"%d"', seq_len(nrow(rebuilt)))))
  values <- matrix(as.numeric(fields[, -1]), nrow(fields))
  c(same_bytes = identical(paper, readLines(file.path(step4, name))),
    same_bytes_rounded = identical(fields[, -1], rounded),
    fields_differing = sum(fields[, -1] != rounded),
    max_difference = max(abs(values - as.matrix(rebuilt))))
}

summary <- as.data.frame(t(vapply(seq_len(nrow(triples)), compare_triple, numeric(4))))
write.csv(cbind(triple = seq_len(nrow(triples)), summary), file.path(work, "r", "summary.csv"),
          row.names = FALSE)
cat(sprintf(paste0("R:      %d files of step 3, %d of step 4; files with the bytes of the paper: ",
                   "%d as written, %d written correctly rounded (%d fields differ); ",
                   "largest difference %.3g\n"),
            length(checksums), nrow(triples), sum(summary$same_bytes),
            sum(summary$same_bytes_rounded), sum(summary$fields_differing),
            max(summary$max_difference)))
