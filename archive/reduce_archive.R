# -----------------------------------------------------------------------------
# File:    archive/reduce_archive.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Tables 1 to 4, Figures 2, 3, 5 to 8.
# Purpose: Copy from the output folders of the original runs the values the
#          table and figure scripts read, into the files kept in archive/.
#          Numbers are copied as text, so no digit changes.
# Inputs:  <runs>/run0204b/{T1,P}_pvals_<s>.csv, {T1,P}_summary_<s>.csv
#          <runs>/real0210a/PctRk_<k>q_0.5.csv
#          <runs>/coh7/CV_ch<k>.csv, <runs>/liver/liver_CV_ch<j>.csv
# Outputs: archive/run0204b/table{1,2,3}_pvalues.csv, fig3_lambdas.csv,
#          archive/real0210a/table4_pvalues.csv,
#          archive/coh7/curves_eeg<k>.csv.gz, archive/liver/curves_liver<j>.csv.gz
# Usage:   Rscript archive/reduce_archive.R /path/to/original/run/folders
# -----------------------------------------------------------------------------

source("R/setup.R")

runs_dir <- commandArgs(trailingOnly = TRUE)[1]
read_text <- function(...) {
  read.csv(file.path(runs_dir, ...), colClasses = "character", check.names = FALSE)[-1]
}
write_text <- function(table, file) {
  dir.create(dirname(file), showWarnings = FALSE)
  write.csv(table, if (endsWith(file, ".gz")) gzfile(file) else file,
            row.names = FALSE, quote = FALSE)
}
P_COLUMNS <- c("p_gcv", "p_loo", "p_ga_gcv", "p_ga_loo")

# Tables 1 to 3: the four p-values of every replication of the 12 settings.
tables <- list(table1 = list(prefix = "T1", beta_fixed = FALSE),
               table2 = list(prefix = "P", beta_fixed = TRUE),
               table3 = list(prefix = "P", beta_fixed = FALSE))
for (name in names(tables)) {
  scenarios <- simulation_scenarios(tables[[name]]$beta_fixed)$scenario
  pvalues <- do.call(rbind, lapply(scenarios, function(s) {
    p <- t(read_text("run0204b", sprintf("%s_pvals_%d.csv", tables[[name]]$prefix, s)))
    data.frame(scenario = s, run = seq_len(nrow(p)), setNames(as.data.frame(p), P_COLUMNS))
  }))
  write_text(pvalues, archive_path("run0204b", paste0(name, "_pvalues.csv")))
}

# Figure 3: penalties chosen by GCV and LOOCV for the observed data (b = 1)
# and every permutation of replication 500.
panels <- data.frame(panel = c("a", "b", "c", "d"),
                     file = c("T1_summary_14", "T1_summary_17", "P_summary_32", "P_summary_35"))
lambdas <- do.call(rbind, lapply(seq_len(nrow(panels)), function(i) {
  summary <- read_text("run0204b", paste0(panels$file[i], ".csv"))
  data.frame(panel = panels$panel[i], b = summary[[1]],
             GCV = summary[[ncol(summary) - 1]], LOOCV = summary[[ncol(summary)]])
}))
write_text(lambdas, archive_path("run0204b", "fig3_lambdas.csv"))

# Table 4: p-values at the chosen penalties of the observed data (b = 1) and
# every permutation, and the chosen penalties.
pvalues <- do.call(rbind, lapply(TABLE4_TRIPLES, function(k) {
  ranks <- read_text("real0210a", sprintf("PctRk_%dq_0.5.csv", k))
  data.frame(triple = k, b = ranks$b,
             setNames(ranks[c("Pval_b_gcv", "Pval_b_loo", "Pval_g_gcv", "Pval_g_loo")], P_COLUMNS),
             lambda_gcv = ranks$l_gcv, lambda_loo = ranks$l_loo)
}))
write_text(pvalues, archive_path("real0210a", "table4_pvalues.csv"))

# Figures 2 and 5 to 8: power and cross-validation curves over the penalties.
curves <- function(source, target) {
  table <- read_text(source)
  names(table) <- sub("l_opt$", "", names(table), fixed = TRUE)
  write_text(table[CURVE_COLUMNS], target)
}
EEG_FIGURE_TRIPLES <- c(160, 248, 191, 183, 184, 220, 148, 177, 22, 157, 31, 151,   # Figure 5
                        181, 231, 200, 144, 159, 136, 223, 205, 25, 179, 219, 123)  # Figure 6
for (k in EEG_FIGURE_TRIPLES) {
  curves(file.path("coh7", sprintf("CV_ch%d.csv", k)), archived_curve_file("eeg", k))
}
for (j in seq_along(LIVER_LABELS)) {
  curves(file.path("liver", sprintf("liver_CV_ch%d.csv", j)), archived_curve_file("liver", j))
}
