# -----------------------------------------------------------------------------
# File:    R/simulation.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Tables 1 to 3, Figure 3.
# Purpose: Simulation study of AdaMantCV: simulated designs, one replication,
#          and the driver that runs replications in parallel, stores their
#          p-values and summarizes them into a rejection-rate table.
# Inputs:  None.
# Outputs: results/<name>_<settings>/scenarioSS_runRRRR.csv (one per
#          replication) and results/<name>_<settings>.csv (the table).
# Usage:   source("R/setup.R"); simulation_study("table1_type1_error", "zero", opts)
# Source:  R_initial.R, R_driver.R, R_simulations.R, R_fcns.R (Feb 2022).
# -----------------------------------------------------------------------------

SIM_QUANTILES   <- c(none = 0, Q1 = 0.25, Q2 = 0.5)
SIM_COVARIANCES <- c(CS = "CSYM", Heteroskedastic = "NCDIAG")
SIM_DIMENSIONS  <- c(500, 1000)

make_covariance <- function(type, p, rho) {
  switch(type,
         CSYM = {
           Sigma <- matrix(rho, nrow = p, ncol = p)
           diag(Sigma) <- 1
           Sigma
         },
         NCDIAG = diag(log(2:(p + 1)), nrow = p, ncol = p),
         stop("unknown covariance type ", type))
}

# Standardized design X and coefficients beta. beta_type: "zero" (null),
# "fixed" (all equal to mu) or "random" (N(mu, sigma), drawn right after X
# from the same random stream).
simulate_design <- function(n, p, covariance, rho, beta_type, mu, sigma, seed) {
  set.seed(seed)
  X <- scale(MASS::mvrnorm(n, rep(0, p), Sigma = make_covariance(covariance, p, rho)))
  beta <- switch(beta_type,
                 zero = rep(0, p),
                 fixed = rep(mu, p),
                 random = rnorm(p, mu, sigma),
                 stop("unknown beta type ", beta_type))
  list(X = X, beta = beta)
}

# The 12 settings of a table in the row order of the paper, with the scenario
# numbers of the original runs (2 to 18 random or zero beta, 20 to 36 fixed).
simulation_scenarios <- function(beta_fixed) {
  grid <- expand.grid(tau = names(SIM_QUANTILES), p = SIM_DIMENSIONS,
                      covariance = names(SIM_COVARIANCES), stringsAsFactors = FALSE)
  grid$scenario <- 18 * beta_fixed + 9 * (grid$p == 1000) +
    3 * (match(grid$tau, names(SIM_QUANTILES)) - 1) +
    match(grid$covariance, names(SIM_COVARIANCES)) + 1
  grid[c("scenario", "covariance", "p", "tau")]
}

# Replication `run`: a new response from the design, then one AdaMantCV test.
# The permutation seed is the same in every replication, as in the original.
simulation_run <- function(design, run, q, sigma, n_negative, n_positive,
                           n_perms, seed, upper = 1000) {
  set.seed(10 * run + seed)
  Y <- design$X %*% design$beta + rnorm(nrow(design$X), mean = 0, sd = sigma)
  adamant_cv(design$X, Y, q, n_negative, n_positive, n_perms, seed, upper)
}

# Run the replications asked for in `opts` (runs, scenarios, perms, negative,
# positive, cores) that are not yet stored, then summarize every stored one.
simulation_study <- function(name, beta_type, opts) {
  scenarios <- simulation_scenarios(beta_type == "fixed")
  if (opts$scenarios != "all") {
    scenarios <- scenarios[scenarios$scenario %in% parse_index_set(opts$scenarios), ]
  }
  tag <- settings_tag(opts)
  run_dir <- results_path(paste0(name, "_", tag))
  dir.create(run_dir, showWarnings = FALSE)
  run_file <- function(scenario, run) {
    file.path(run_dir, sprintf("scenario%02d_run%04d.csv", scenario, run))
  }

  design_keys <- unique(scenarios[c("covariance", "p")])
  designs <- Map(function(covariance, p) {
    simulate_design(SIM_N, p, SIM_COVARIANCES[[covariance]], SIM_RHO, beta_type,
                    SIM_BETA_MEAN, SIGMA, SEED)
  }, design_keys$covariance, design_keys$p)
  names(designs) <- paste(design_keys$covariance, design_keys$p)

  jobs <- expand.grid(run = parse_index_set(opts$runs, count = TRUE), row = seq_len(nrow(scenarios)))
  jobs <- jobs[!file.exists(run_file(scenarios$scenario[jobs$row], jobs$run)), ]
  log_message("%s: %d replications to run", name, nrow(jobs))
  status <- parallel::mclapply(seq_len(nrow(jobs)), function(i) {
    setting <- scenarios[jobs$row[i], ]
    ada <- simulation_run(designs[[paste(setting$covariance, setting$p)]], jobs$run[i],
                          SIM_QUANTILES[[setting$tau]], SIGMA, opts$negative,
                          opts$positive, opts$perms, SEED, UPPER)
    write.csv(data.frame(scenario = setting$scenario, run = jobs$run[i],
                         p_gcv = ada$p_values[["gcv"]], p_loo = ada$p_values[["loo"]],
                         p_ga_gcv = ada$p_values[["ga_gcv"]],
                         p_ga_loo = ada$p_values[["ga_loo"]]),
              run_file(setting$scenario, jobs$run[i]), row.names = FALSE)
  }, mc.cores = opts$cores, mc.preschedule = FALSE)
  failed <- vapply(status, inherits, logical(1), what = "try-error")
  if (any(failed)) stop(sum(failed), " replications failed: ", status[failed][[1]])

  pvalues <- do.call(rbind, lapply(list.files(run_dir, full.names = TRUE), read.csv))
  rejected <- aggregate(cbind(runs = 1, gcv = p_gcv <= ALPHA, loocv = p_loo <= ALPHA,
                              ga_gcv = p_ga_gcv <= ALPHA, ga_loocv = p_ga_loo <= ALPHA) ~
                          scenario, data = pvalues, FUN = sum)
  rates <- c("gcv", "loocv", "ga_gcv", "ga_loocv")
  rejected[rates] <- rejected[rates] / rejected$runs
  table <- merge(scenarios, rejected, by = "scenario", sort = FALSE)
  table <- table[order(match(table$scenario, scenarios$scenario)), ]
  write.csv(table, results_path(sprintf("%s_%s.csv", name, tag)), row.names = FALSE)
  print(table, row.names = FALSE)
  invisible(table)
}
