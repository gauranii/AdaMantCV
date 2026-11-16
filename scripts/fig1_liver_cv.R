# -----------------------------------------------------------------------------
# File:    scripts/fig1_liver_cv.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Figure 1.
# Purpose: Ten-fold cross-validated MSE of ridge regression against log(lambda)
#          for the ten clinical variables of the liver toxicity data, with 50
#          randomly chosen genes (A) and with all 3116 genes (B).
# Inputs:  <data>/liver/mixOmicsData.csv.
# Outputs: figures/fig1a_liver_cv_subset.pdf, figures/fig1b_liver_cv_full.pdf,
#          results/fig1_minima.csv.
# Usage:   Rscript scripts/fig1_liver_cv.R
# Source:  fig1/CV Ridge.R.
# -----------------------------------------------------------------------------

source("R/setup.R")

library(glmnet)

SUBSET_SEED    <- 889191 + 1  # seed of the published panel A
SUBSET_SIZE    <- 50
SUBSET_LAMBDAS <- 10^seq(4, -2.5, by = -0.001)
# The original set no seed before the fits of panel B, so its random folds are
# lost and the published panel cannot be reproduced fold for fold.
FULL_SEED      <- 889191
FULL_LAMBDAS   <- 10^seq(9, -5, by = -0.1)
N_FOLDS        <- 10
MSE_AXIS_MAX   <- 1.2

# Cross-validated MSE of every clinical variable over the penalties, long format.
cv_curves <- function(X, Y, lambdas) {
  do.call(rbind, lapply(seq_len(ncol(Y)), function(j) {
    fit <- cv.glmnet(X, Y[, j], alpha = 0, nfolds = N_FOLDS, lambda = lambdas)
    data.frame(variable = paste0("Variable", j), log_lambda = log(lambdas),
               cvm = fit$cvm, at_minimum = fit$cvm == min(fit$cvm),
               log_lambda_min = log(fit$lambda.min))
  }))
}

cv_plot <- function(curves) {
  curves$variable <- factor(curves$variable, levels = unique(curves$variable))
  ggplot(curves, aes(x = log_lambda, y = cvm, color = variable)) +
    geom_line(linewidth = 1) + theme_bw() +
    geom_point(data = curves[curves$at_minimum, ], color = "blue") +
    xlab("log(lambda)") + ylab("Cross-validated MSE") +
    theme(legend.position = "none") +
    ylim(NA, MSE_AXIS_MAX)
}

minima <- function(curves, panel) {
  lowest <- curves[curves$at_minimum & curves$log_lambda == curves$log_lambda_min, ]
  data.frame(panel = panel, variable = lowest$variable,
             log_lambda_min = lowest$log_lambda_min, min_cvm = lowest$cvm)
}

liver <- read_liver()

set.seed(SUBSET_SEED)
genes <- sample(seq_len(ncol(liver$X)), SUBSET_SIZE, replace = FALSE)
subset_curves <- cv_curves(liver$X[, genes], liver$Y, SUBSET_LAMBDAS)

set.seed(FULL_SEED)
full_curves <- cv_curves(liver$X, liver$Y, FULL_LAMBDAS)

save_figure(cv_plot(subset_curves), "fig1a_liver_cv_subset.pdf")
save_figure(cv_plot(full_curves), "fig1b_liver_cv_full.pdf")
write.csv(rbind(minima(subset_curves, "A"), minima(full_curves, "B")),
          results_path("fig1_minima.csv"), row.names = FALSE)
