# -----------------------------------------------------------------------------
# File:    scripts/fig1_liver_cv.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Figure 1.
# Purpose: Ten-fold cross-validated MSE of ridge regression against log(lambda)
#          for the ten clinical variables of the liver toxicity data, with 50
#          randomly chosen genes (A) and with all 3116 genes (B).
# Inputs:  data/liver/mixOmicsData.csv.
# Outputs: figures/fig1a_liver_cv_subset.pdf, figures/fig1b_liver_cv_full.pdf,
#          results/fig1_minima.csv.
# Usage:   Rscript scripts/fig1_liver_cv.R
# Source:  fig1/CV Ridge.R.
# -----------------------------------------------------------------------------

source("R/setup.R")

SUBSET_SEED    <- 889191 + 1  # seed of panel A
SUBSET_SIZE    <- 50
SUBSET_LAMBDAS <- 10^seq(4, -2.5, by = -0.001)
# Panel B continues the random number stream of the original session: seed
# 889191 + 2, one sample of 50 genes, then eleven fold assignments.
FULL_SEED          <- 889191 + 2
FULL_SKIPPED_FOLDS <- 11
FULL_LAMBDAS   <- 10^seq(9, -5, by = -0.1)
N_FOLDS        <- 10
MSE_AXIS_MAX   <- 1.2

opts <- parse_args(list())  # no options; Figure 1 is always computed from the data

liver <- read_liver(PATHS$liver)
n_genes <- ncol(liver$X)

set.seed(SUBSET_SEED)
genes <- sample(seq_len(n_genes), SUBSET_SIZE, replace = FALSE)
subset_curves <- cv_curves(liver$X[, genes], liver$Y, SUBSET_LAMBDAS, N_FOLDS)

set.seed(FULL_SEED)
invisible(sample(seq_len(n_genes), SUBSET_SIZE, replace = FALSE))
skip_fold_assignments(FULL_SKIPPED_FOLDS, nrow(liver$X), N_FOLDS)
full_curves <- cv_curves(liver$X, liver$Y, FULL_LAMBDAS, N_FOLDS)

save_figure(cv_plot(subset_curves, MSE_AXIS_MAX), "fig1a_liver_cv_subset.pdf",
            PATHS$figures, PAPER$figure_size)
save_figure(cv_plot(full_curves, MSE_AXIS_MAX), "fig1b_liver_cv_full.pdf",
            PATHS$figures, PAPER$figure_size)
write.csv(rbind(cv_minima(subset_curves, "A"), cv_minima(full_curves, "B")),
          output_file(PATHS$results, "fig1_minima.csv"), row.names = FALSE)
