# -----------------------------------------------------------------------------
# File:    checks/check_equivalence.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Run the ORIGINAL 2021/2022 functions and the rebuilt ones on the same
#          inputs and seeds (small permutation counts, so it runs in minutes)
#          and compare every output.
# Inputs:  original scripts (env ADAMANT_ORIGINAL: run0204b/, coh7/RPC_fcns.R),
#          EEG and liver data (env ADAMANT_DATA).
# Outputs: results/checks/equivalence.csv; prints one line per comparison.
# Usage:   ADAMANT_ORIGINAL=... ADAMANT_DATA=... Rscript checks/check_equivalence.R
# Needs:   dplyr, purrr, tibble, magrittr, reshape, MASS (for the original code).
# -----------------------------------------------------------------------------

source("R/setup.R")
suppressMessages({
  library(dplyr); library(purrr); library(tibble); library(magrittr); library(reshape); library(MASS)
})
original <- Sys.getenv("ADAMANT_ORIGINAL")

# Source original files into their own environment, without library() calls.
load_original <- function(files) {
  env <- new.env()
  for (f in files) {
    code <- readLines(file.path(original, f), warn = FALSE)
    eval(parse(text = code[!grepl("^(library|source)", code)]), envir = env)
  }
  env$upper <- UPPER
  env$alpha <- ALPHA
  env
}
feb <- load_original(c("run0204b/R_fcns.R", "run0204b/R_adaMantCV.R",
                       "run0204b/R_simulations.R"))
dec <- load_original("coh7/RPC_fcns.R")

results <- list()
record <- function(what, rebuilt, original, tolerance = 0) {
  diff <- max(abs(as.numeric(rebuilt) - as.numeric(original)))
  ok <- length(original) %in% c(1, length(rebuilt)) && diff <= tolerance
  results[[length(results) + 1]] <<- data.frame(what, n = length(rebuilt),
                                                 max_abs_diff = diff, tolerance, ok)
  cat(sprintf("%-60s n=%-7d max|diff|=%-10.3g %s\n", what, length(rebuilt), diff,
              if (ok) "OK" else "DIFFERENT"))
}
relative_or_absolute <- function(new, old) {
  ifelse(abs(old) < 1e-12, new - old, new / old - 1)
}
quiet <- function(expr) { capture.output(value <- expr); value }

# ----- 1. adamant_cv on real EEG data (Table 4 code), 100 permutations -----
eeg <- read_eeg_triple(11)
for (q in c(0.25, 0.5)) {
  new <- adamant_cv(eeg$X, eeg$y, q, 50, 200, 100, SEED)
  old <- quiet(feb$adamant_cv(eeg$X, as.matrix(eeg$y), q, 50, 200, 1, 100, SEED))
  tag <- sprintf("adamant_cv triple 11, q = %.2f, 100 perms:", q)
  record(paste(tag, "lambda grid"), new$lambdas, old$lambdas)
  # At lambda = 0 the observed fit is exact, so its score is rounding noise
  # (about 1e-24); scores below 1e-12 are compared on the absolute scale.
  # Tolerance 1e-8: all permutations are fitted in one matrix product.
  record(paste(tag, "GCV scores (relative; absolute below 1e-12)"),
         relative_or_absolute(new$gcv, old$gcv), 0, 1e-8)
  record(paste(tag, "LOOCV scores (relative; absolute below 1e-12)"),
         relative_or_absolute(new$loo, old$loo), 0, 1e-8)
  record(paste(tag, "chosen lambdas GCV, LOO"),
         c(new$lambda_gcv, new$lambda_loo), old$summary_cv[, c(2, 4)])
  record(paste(tag, "p-values (4)"), new$p_values,
         c(old$Pval_gcv, old$Pval_loo, old$Pval_gam_gcv, old$Pval_gam_loo), 1e-12)
}

# ----- 2. Simulation replications (Tables 1 to 3 code), 30 permutations -----
for (beta_type in c("zero", "fixed", "random")) {
  design <- simulate_design(SIM_N, 500, "CSYM", SIM_RHO, beta_type, SIM_BETA_MEAN,
                            SIGMA, SEED)
  new <- sapply(1:2, function(r) {
    simulation_run(design, r, 0.25, SIGMA, 50, 200, 30, SEED, UPPER)$p_values
  })
  old <- if (beta_type == "zero") {
    quiet(feb$type1_simulate(SIM_N, 500, 0.25, "CSYM", SIM_RHO, NULL, TRUE, TRUE, 0,
                             SIGMA, 0, 50, 200, 1, 30, 2, ALPHA, SEED))
  } else {
    quiet(feb$power_simulate(SIM_N, 500, 0.25, "CSYM", SIM_RHO, NULL, TRUE,
                             beta_type == "random", SIM_BETA_MEAN, SIGMA, 0, 50, 200,
                             1, 30, 2, ALPHA, SEED))
  }
  record(sprintf("simulation CS p = 500 Q1, beta %s, runs 1-2, 30 perms: p-values",
                 beta_type), new, old$Pval_summary, 1e-12)
}

# ----- 3. Power and GCV curves (Figs 2, 5 to 8 code), 3 runs, 30 permutations -----
curve_inputs <- list(eeg191 = c(read_eeg_triple(191), EEG_GRID),
                     liver1 = c(list(X = read_liver()$X, y = read_liver()$Y[, 1]),
                                LIVER_GRID))
for (name in names(curve_inputs)) {
  d <- curve_inputs[[name]]
  new <- power_curve(d$X, d$y, d$negative, d$positive, SIGMA, 30, 3, ALPHA, SEED)
  old <- quiet(dec$power_ridge(d$X, as.matrix(d$y), 1, SIGMA, 30, 3, d$negative,
                               d$positive, ALPHA, SEED))
  record(sprintf("power curve %s, 3 runs, 30 perms: lambdas", name), new$lambdas,
         old$lambdas)
  record(sprintf("power curve %s: power", name), new$power, old$power[1, ])
  record(sprintf("power curve %s: power MOM", name), new$power_MOM, old$power[2, ])
  record(sprintf("power curve %s: permutation p-values (as mean -log10)", name),
         new$ave_logpval,
         colMeans(-log10(old$pvals + sqrt(.Machine$double.eps))), 1e-12)
  new_cv <- gcv_curve(d$X, d$y, d$negative, d$positive, SIGMA, SEED)
  old_dir <- tempfile(); dir.create(old_dir)
  old_cv <- local({
    here <- setwd(old_dir); on.exit(setwd(here))
    quiet(dec$lambda_optimal(d$X, as.matrix(d$y), 1, name, SIGMA, d$negative,
                             d$positive, SEED))
  })
  record(sprintf("GCV curve %s (relative)", name), new_cv$gcv / old_cv$gcv - 1, 0, 1e-10)
  record(sprintf("LOOCV curve %s (relative)", name), new_cv$loo / old_cv$loo - 1, 0, 1e-10)
}

dir.create(results_path("checks"), showWarnings = FALSE)
out <- do.call(rbind, results)
write.csv(out, results_path("checks/equivalence.csv"), row.names = FALSE)
cat(sprintf("equivalence: %d/%d comparisons OK\n", sum(out$ok), nrow(out)))
