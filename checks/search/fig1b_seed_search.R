# -----------------------------------------------------------------------------
# File:    checks/search/fig1b_seed_search.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
#          Figure 1B.
# Purpose: The original code set no seed before the cross-validation folds of
#          Figure 1B. Try the seeds and session states the original script
#          could have run from, and store the ten minima of each, so that the
#          one closest to the published panel can be chosen.
# Inputs:  <data>/liver/mixOmicsData.csv.
# Outputs: checks/search/fig1b_candidates.csv.
# Usage:   Rscript checks/search/fig1b_seed_search.R   (about 15 minutes)
# Source:  fig1/CV Ridge.R.
# -----------------------------------------------------------------------------

source("R/setup.R")
library(glmnet)

liver <- read_liver()
FULL   <- 10^seq(9, -5, by = -0.1)
SUBSET <- 10^seq(4, -2.5, by = -0.001)

minima <- function(X, lambdas) sapply(seq_len(ncol(liver$Y)), function(j) {
  fit <- cv.glmnet(X, liver$Y[, j], alpha = 0, nfolds = 10, lambda = lambdas)
  c(log(fit$lambda.min), min(fit$cvm))
})

found <- list()
# Plain seeds: those of the original scripts and common choices
for (s in c(889191, 889192, 889193, 889190, 8891, 1, 42, 123, 1234, 2021, 2022)) {
  set.seed(s)
  found[[paste0("seed", s)]] <- minima(liver$X, FULL)
}
# The single fit of variable 1 that comes before the loop in the original
for (s in c(889191, 889192, 889193, 8891)) {
  set.seed(s)
  invisible(cv.glmnet(liver$X, liver$Y[, 1], alpha = 0, lambda = FULL))
  found[[paste0("seed", s, "+fit1")]] <- minima(liver$X, FULL)
}
# The state left by the panel A loop (k = 1, then k = 1 and 2)
for (ks in list(1, 1:2)) {
  for (k in ks) {
    set.seed(889191 + k)
    genes <- sample(1:ncol(liver$X), 50)
    invisible(minima(liver$X[, genes], SUBSET))
  }
  found[[paste0("after_subset_k", paste(ks, collapse = ""))]] <- minima(liver$X, FULL)
}

table <- do.call(rbind, lapply(names(found), function(name) {
  data.frame(candidate = name, variable = 1:10, log_lambda_min = found[[name]][1, ],
             min_cvm = found[[name]][2, ])
}))
write.csv(table, file.path(REPO_ROOT, "checks", "search", "fig1b_candidates.csv"),
          row.names = FALSE)
