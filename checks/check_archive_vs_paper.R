# -----------------------------------------------------------------------------
# File:    checks/check_archive_vs_paper.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Tie the archived 2022 run outputs to the numbers printed in the paper
#          (Tables 1 to 4), so later checks can compare reruns with the archive.
# Inputs:  checks/paper_values.csv; archive (env ADAMANT_ARCHIVE) with
#          run0204b/{T1_type1,P_power,T1_pvals,P_pvals}_<s>.csv and
#          real0210a/AdaPval_<k>q_0.5.csv
# Outputs: results/checks/archive_vs_paper.csv; prints a summary line
# Usage:   ADAMANT_ARCHIVE=/path/to/archive Rscript checks/check_archive_vs_paper.R
# -----------------------------------------------------------------------------

archive <- Sys.getenv("ADAMANT_ARCHIVE", "archive")
paper   <- read.csv("checks/paper_values.csv", stringsAsFactors = FALSE)
methods <- c("gcv", "loo", "ga_gcv", "ga_loo")

# Scenario number of a paper row in Run0204b (see R_initial.R settings grid).
scenario <- function(table, covariance, p, tau) {
  18 * (table == 2) + 9 * (p == 1000) + 3 * match(tau, c("none", "Q1", "Q2")) - 3 +
    match(covariance, c("CS", "Het")) + 1
}

read_col <- function(file) read.csv(file)[[2]]

archive_value <- function(r) {
  m <- match(r$method, methods)
  if (r$table == 4) {
    k <- c(11, 27, 70, 178, 248)[r$row]
    return(read_col(file.path(archive, "real0210a", sprintf("AdaPval_%dq_0.5.csv", k)))[m])
  }
  s <- scenario(r$table, r$covariance, as.numeric(r$p), r$tau)
  prefix <- if (r$table == 1) "T1" else "P"
  agg <- if (r$table == 1) "type1" else "power"
  stored <- read_col(file.path(archive, "run0204b", sprintf("%s_%s_%d.csv", prefix, agg, s)))[m]
  pvals <- as.matrix(read.csv(file.path(archive, "run0204b", sprintf("%s_pvals_%d.csv", prefix, s)))[, -1])
  c(stored, mean(pvals[m, ] < 0.05), mean(pvals[m, ] <= 0.05))
}

rows <- lapply(seq_len(nrow(paper)), function(i) {
  r <- paper[i, ]
  v <- archive_value(r)
  data.frame(r, archive = v[1],
             from_runs_lt = if (length(v) > 1) v[2] else NA,
             from_runs_le = if (length(v) > 1) v[3] else NA)
})
out <- do.call(rbind, rows)
out$match <- abs(round(out$archive, 3) - out$value) < 1e-9
dir.create("results/checks", recursive = TRUE, showWarnings = FALSE)
write.csv(out, "results/checks/archive_vs_paper.csv", row.names = FALSE)

sim <- out$table < 4
cat(sprintf("archive vs paper: %d/%d cells equal (3 decimals)\n", sum(out$match), nrow(out)))
cat(sprintf("Tables 1-3 from the 500 per-run p-values: %d/%d with p < 0.05, %d/%d with p <= 0.05\n",
            sum(abs(out$from_runs_lt[sim] - out$value[sim]) < 1e-9), sum(sim),
            sum(abs(out$from_runs_le[sim] - out$value[sim]) < 1e-9), sum(sim)))
