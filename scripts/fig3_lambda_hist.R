# -----------------------------------------------------------------------------
# File:    scripts/fig3_lambda_hist.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Figure 3.
# Purpose: Histograms of the ridge penalties chosen by GCV and LOOCV over the
#          permutations of one simulated data set (replication 500, compound
#          symmetric covariance, p = 1000), under beta = 0 (A, B) and beta = 1
#          (C, D), with tau = Q1 (A, C) and tau = Q2 (B, D).
# Inputs:  archive/run0204b/fig3_lambdas.csv (chosen penalties of the archived
#          run); none with --rerun (simulated data).
# Outputs: figures/Fig3{A,B}_simulation_lambda-histogram_null_threshold-{Q1,median}.pdf,
#          figures/Fig3{C,D}_simulation_lambda-histogram_alternative_threshold-{Q1,median}.pdf;
#          with --rerun also results/fig3_lambdas_<settings>.csv.
# Usage:   Rscript scripts/fig3_lambda_hist.R [--rerun] [--perms=1000] [--negative=50]
#          [--positive=200] [--cores=4]
#          --rerun: simulate the data and choose the penalties here.
# Source:  run0204b (T1_summary_*, P_summary_* of scenarios 11, 14, 17 and
#          29, 32, 35), R_Lambda.R (hist2, hist4).
# -----------------------------------------------------------------------------

source("R/setup.R")

REPLICATION <- 500     # the replication shown
DIMENSION   <- 1000
COVARIANCE  <- "CSYM"
BIN_WIDTH   <- 0.25
PANELS <- data.frame(panel = c("a", "b", "c", "d"),
                     beta_type = c("zero", "zero", "fixed", "fixed"),
                     tau = c("Q1", "Q2", "Q1", "Q2"))
# Words of the file names
HYPOTHESIS <- c(zero = "null", fixed = "alternative")
THRESHOLD  <- c(Q1 = "Q1", Q2 = "median")

opts <- parse_args(test_options(PAPER))

lambdas <- chosen_penalties(PANELS, REPLICATION, DIMENSION, COVARIANCE, opts, PAPER, PATHS)

# The histograms show the permutations (b = 1 is the observed data).
lambdas <- lambdas[lambdas$b >= 2, ]

for (i in seq_len(nrow(PANELS))) {
  hist <- penalty_histogram(lambdas[lambdas$panel == PANELS$panel[i], ], PANELS$tau[i],
                            BIN_WIDTH)
  name <- sprintf("Fig3%s_simulation_lambda-histogram_%s_threshold-%s.pdf",
                  toupper(PANELS$panel[i]), HYPOTHESIS[[PANELS$beta_type[i]]],
                  THRESHOLD[[PANELS$tau[i]]])
  save_figure(hist, name, PATHS$figures, PAPER$figure_size)
}
