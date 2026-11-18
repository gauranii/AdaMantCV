# -----------------------------------------------------------------------------
# File:    checks/search/fig1b_seed_search2.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
#          Figure 1B.
# Purpose: Second search for the folds of Figure 1B: the session states that
#          parts of the original fig1/CV Ridge.R leave, with and without its
#          single fit of variable 1, under both sample.kind settings.
# Inputs:  <data>/liver/mixOmicsData.csv.
# Outputs: the csv named on the command line (same columns as the first search).
# Usage:   Rscript checks/search/fig1b_seed_search2.R <out.csv>   (about 10 minutes)
# -----------------------------------------------------------------------------
source("R/setup.R")
library(glmnet)
liver <- read_liver()
X <- liver$X; Y <- liver$Y; N <- nrow(X)
FULL <- 10^seq(9, -5, by = -0.1)
OUT <- commandArgs(TRUE)[1]

# The RNG draw one cv.glmnet call makes (its folds); checked against a real fit below.
fold_draw <- function() invisible(sample(rep(seq(10), length = N)))
set.seed(5); invisible(cv.glmnet(X[, 1:20], Y[, 1], alpha = 0, lambda = FULL)); a <- .Random.seed
set.seed(5); fold_draw(); b <- .Random.seed
stopifnot(identical(a, b))

panelB <- function(fit1) {
  if (fit1) fold_draw()            # the single fit of variable 1 before the loop
  sapply(1:10, function(j) {
    f <- cv.glmnet(X, Y[, j], alpha = 0, nfolds = 10, lambda = FULL)
    c(log(f$lambda.min), min(f$cvm))
  })
}
states <- list(
  seed8891 = function() set.seed(8891),
  seed889191 = function() set.seed(889191),
  seed889192 = function() set.seed(889192),
  seed889193 = function() set.seed(889193),
  seed1234 = function() set.seed(1234),
  after_panelA = function() { set.seed(889193); sample(1:ncol(X), 50); for (j in 1:10) fold_draw() },
  after_perm_loop = function() { set.seed(2234); sample(1:N); fold_draw() },
  after_whole_file = function() { set.seed(1234); for (b in 1:1000) sample(1:N) }
)
found <- list()
for (kind in c("Rejection", "Rounding")) for (nm in names(states)) for (fit1 in c(FALSE, TRUE)) {
  suppressWarnings(RNGkind(sample.kind = kind))
  states[[nm]]()
  key <- paste0(nm, if (fit1) "+fit1" else "", "|", kind)
  found[[key]] <- panelB(fit1)
  cat(key, "\n")
}
tab <- do.call(rbind, lapply(names(found), function(k) data.frame(candidate = k,
  variable = 1:10, log_lambda_min = found[[k]][1, ], min_cvm = found[[k]][2, ])))
write.csv(tab, OUT, row.names = FALSE)
