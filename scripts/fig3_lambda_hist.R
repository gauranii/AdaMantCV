# -----------------------------------------------------------------------------
# File:    scripts/fig3_lambda_hist.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Figure 3.
# Purpose: Histograms of the ridge penalties chosen by GCV and LOOCV over the
#          permutations of one simulated data set (replication 500, compound
#          symmetric covariance, p = 1000), under beta = 0 (A, B) and beta = 1
#          (C, D), with tau = Q1 (A, C) and tau = Q2 (B, D).
# Inputs:  None (simulated data).
# Outputs: results/fig3_lambdas_<settings>.csv, figures/fig3{a,b,c,d}_lambda_hist.pdf.
# Usage:   Rscript scripts/fig3_lambda_hist.R [--perms=1000] [--negative=50]
#          [--positive=200] [--cores=4]
#          --method=paper: the method as the paper states it (README).
# Source:  run0204b (T1_summary_*, P_summary_* of scenarios 11, 14, 17 and
#          29, 32, 35), R_Lambda.R (hist2, hist4).
# -----------------------------------------------------------------------------

source("R/setup.R")

REPLICATION <- 500     # the last replication, the one whose lambdas were saved
DIMENSION   <- 1000
COVARIANCE  <- "CSYM"
BIN_WIDTH   <- 0.25
PANELS <- data.frame(panel = c("a", "b", "c", "d"),
                     beta_type = c("zero", "zero", "fixed", "fixed"),
                     tau = c("Q1", "Q2", "Q1", "Q2"))

opts <- parse_args(list(method = "original",
                        perms = N_PERMS, negative = 50, positive = 200, cores = 4))

# Penalty chosen by GCV and LOOCV for the observed data (b = 1) and every
# permutation of one replication, for each row of PANELS.
chosen_lambdas <- function() {
  designs <- lapply(c(zero = "zero", fixed = "fixed"), function(beta_type) {
    simulate_design(SIM_N, DIMENSION, COVARIANCE, SIM_RHO, beta_type, SIM_BETA_MEAN,
                    SIGMA, SEED)
  })
  chosen <- parallel::mclapply(seq_len(nrow(PANELS)), function(i) {
    ada <- simulation_run(designs[[PANELS$beta_type[i]]], REPLICATION,
                          SIM_QUANTILES[[PANELS$tau[i]]], SIGMA, opts$negative,
                          opts$positive, opts$perms, SEED, UPPER)
    data.frame(PANELS[i, ], b = seq_along(ada$lambda_gcv), GCV = ada$lambda_gcv,
               LOOCV = ada$lambda_loo, row.names = NULL)
  }, mc.cores = opts$cores)
  do.call(rbind, chosen)
}

lambdas <- cached_csv(results_path(sprintf("fig3_lambdas_%s.csv", settings_tag(opts))),
                      chosen_lambdas)

# The observed data (b = 1) are left out, as in the published figure.
lambdas <- lambdas[lambdas$b >= 2, ]

for (i in seq_len(nrow(PANELS))) {
  shown <- lambdas[lambdas$panel == PANELS$panel[i], ]
  long <- data.frame(variable = rep(paste(c("GCV", "LOOCV"), "tau =", PANELS$tau[i]),
                                    each = nrow(shown)),
                     value = c(shown$GCV, shown$LOOCV))
  is_gcv <- startsWith(long$variable, "GCV")
  hist <- ggplot(long, aes(x = value, color = variable)) +
    geom_histogram(data = long[is_gcv, ], fill = "red", binwidth = BIN_WIDTH,
                   alpha = 0.2) +
    geom_histogram(data = long[!is_gcv, ], fill = "white", binwidth = BIN_WIDTH,
                   alpha = 0.2) +
    xlab("optimal ridge penalty estimates") + ylab("") + theme_bw() +
    theme(legend.position = "bottom")
  save_figure(hist, sprintf("fig3%s_lambda_hist.pdf", PANELS$panel[i]))
}
