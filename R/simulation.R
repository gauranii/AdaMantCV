# -----------------------------------------------------------------------------
# File:    R/simulation.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Tables 1 to 3, Figure 3.
# Purpose: Simulation study of AdaMantCV: simulated designs, one replication,
#          the driver that runs replications in parallel and stores their
#          p-values, and the rejection-rate table of a set of p-values.
# Inputs:  archive/run0204b/table<N>_pvalues.csv, fig3_lambdas.csv (unless
#          rerun).
# Outputs: results/<name>.csv (the table); in a rerun also
#          results/<name>_<settings>/scenarioSS_runRRRR.csv (one per
#          replication).
# Usage:   source("R/setup.R")
#          simulation_study("table1_type1_error", "zero", opts, PAPER, PATHS)
# Source:  R_initial.R, R_driver.R, R_simulations.R, R_fcns.R (Feb 2022).
# -----------------------------------------------------------------------------

# Simulation design of Tables 1 to 3 and Figure 3, passed as `design`.
SIMULATION_DESIGN <- list(
  n           = 350,                                         # sample size
  rho         = 0.025,                                       # correlation, compound symmetric
  beta_mean   = 1,                                           # value or mean of the coefficients
  quantiles   = c(none = 0, Q1 = 0.25, Q2 = 0.5),            # truncation tau
  covariances = c(CS = "CSYM", Heteroskedastic = "NCDIAG"),
  dimensions  = c(500, 1000)                                 # numbers of predictors p
)

#' Covariance matrix of the simulated predictors.
#'
#' @param type "CSYM" (compound symmetric) or "NCDIAG" (diagonal with
#'   variances log(2), ..., log(p + 1)).
#' @param p Dimension.
#' @param rho Correlation of the compound symmetric matrix.
#' @return p x p matrix.
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

#' Simulated standardized design and coefficients.
#'
#' @param n,p Numbers of observations and predictors.
#' @param covariance,rho As for `make_covariance()`.
#' @param beta_type "zero" (null), "fixed" (all equal to `mu`) or "random"
#'   (N(mu, sigma), drawn right after X from the same random stream).
#' @param mu,sigma Mean and sd of the coefficients.
#' @param seed Seed of the design.
#' @return List with `X` (n x p, standardized) and `beta` (length p).
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

#' The 12 settings of a table in the row order of the paper.
#'
#' @param beta_fixed TRUE for fixed coefficients (scenario numbers 20 to 36
#'   of the original runs), FALSE for zero or random ones (2 to 18).
#' @param design Simulation design (`quantiles`, `covariances`, `dimensions`).
#' @return Data frame with columns `scenario`, `covariance`, `p`, `tau`.
simulation_scenarios <- function(beta_fixed, design) {
  grid <- expand.grid(tau = names(design$quantiles), p = design$dimensions,
                      covariance = names(design$covariances), stringsAsFactors = FALSE)
  grid$scenario <- 18 * beta_fixed + 9 * (grid$p == 1000) +
    3 * (match(grid$tau, names(design$quantiles)) - 1) +
    match(grid$covariance, names(design$covariances)) + 1
  grid[c("scenario", "covariance", "p", "tau")]
}

#' One replication: a new response from the design, then one AdaMantCV test.
#'
#' The permutation seed is the same in every replication, as in the original.
#'
#' @param simulated Output of `simulate_design()`.
#' @param run Replication number; the response is drawn with seed
#'   `10 * run + seed`.
#' @param q,n_negative,n_positive,n_perms,seed,upper As for `adamant_cv()`.
#' @param sigma Sd of the noise.
#' @return Output of `adamant_cv()`.
simulation_run <- function(simulated, run, q, sigma, n_negative, n_positive,
                           n_perms, seed, upper = 1000) {
  set.seed(10 * run + seed)
  Y <- simulated$X %*% simulated$beta + rnorm(nrow(simulated$X), mean = 0, sd = sigma)
  adamant_cv(simulated$X, Y, q, n_negative, n_positive, n_perms, seed, upper)
}

#' Simulated designs of a set of scenarios, one per covariance and dimension.
#'
#' @param scenarios Rows of `simulation_scenarios()`.
#' @param beta_type As for `simulate_design()`.
#' @param settings Settings of the paper (`seed`, `sigma`).
#' @param design Simulation design.
#' @return List of `simulate_design()` outputs named "<covariance> <p>".
simulate_designs <- function(scenarios, beta_type, settings, design) {
  keys <- unique(scenarios[c("covariance", "p")])
  simulated <- Map(function(covariance, p) {
    simulate_design(design$n, p, design$covariances[[covariance]], design$rho, beta_type,
                    design$beta_mean, settings$sigma, settings$seed)
  }, keys$covariance, keys$p)
  setNames(simulated, paste(keys$covariance, keys$p))
}

#' Run the replications that are not yet stored and return all stored ones.
#'
#' Each replication seeds itself, so the result does not depend on the number
#' of parallel processes or on the order of the replications.
#'
#' @param run_dir Folder with one csv file per replication.
#' @param beta_type As for `simulate_design()`.
#' @param scenarios Rows of `simulation_scenarios()` to run.
#' @param opts Options: `runs`, `perms`, `negative`, `positive`, `cores`.
#' @param settings Settings of the paper (`seed`, `sigma`, `upper`).
#' @param design Simulation design.
#' @return Data frame with columns `scenario`, `run` and the four p-values
#'   `p_gcv`, `p_loo`, `p_ga_gcv`, `p_ga_loo` of every stored replication.
run_replications <- function(run_dir, beta_type, scenarios, opts, settings, design) {
  dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
  run_file <- function(scenario, run) {
    file.path(run_dir, sprintf("scenario%02d_run%04d.csv", scenario, run))
  }
  simulated <- simulate_designs(scenarios, beta_type, settings, design)

  jobs <- expand.grid(run = parse_index_set(opts$runs, count = TRUE),
                      row = seq_len(nrow(scenarios)))
  jobs <- jobs[!file.exists(run_file(scenarios$scenario[jobs$row], jobs$run)), ]
  log_message("%s: %d replications to run", basename(run_dir), nrow(jobs))
  status <- parallel::mclapply(seq_len(nrow(jobs)), function(i) {
    setting <- scenarios[jobs$row[i], ]
    ada <- simulation_run(simulated[[paste(setting$covariance, setting$p)]], jobs$run[i],
                          design$quantiles[[setting$tau]], settings$sigma, opts$negative,
                          opts$positive, opts$perms, settings$seed, settings$upper)
    write.csv(data.frame(scenario = setting$scenario, run = jobs$run[i],
                         p_gcv = ada$p_values[["gcv"]], p_loo = ada$p_values[["loo"]],
                         p_ga_gcv = ada$p_values[["ga_gcv"]],
                         p_ga_loo = ada$p_values[["ga_loo"]]),
              run_file(setting$scenario, jobs$run[i]), row.names = FALSE)
  }, mc.cores = opts$cores, mc.preschedule = FALSE)
  failed <- vapply(status, inherits, logical(1), what = "try-error")
  if (any(failed)) stop(sum(failed), " replications failed: ", status[failed][[1]])

  do.call(rbind, lapply(list.files(run_dir, full.names = TRUE), read.csv))
}

#' Rejection rates of the four tests in every setting.
#'
#' A test rejects when its stored p-value is below `alpha`.
#'
#' @param scenarios Rows of `simulation_scenarios()`.
#' @param pvalues Data frame as returned by `run_replications()`.
#' @param alpha Significance level.
#' @return `scenarios` with the columns `runs` (number of replications) and
#'   the rates `gcv`, `loocv`, `ga_gcv`, `ga_loocv`.
rejection_rates <- function(scenarios, pvalues, alpha) {
  rejected <- aggregate(cbind(runs = 1, gcv = p_gcv < alpha, loocv = p_loo < alpha,
                              ga_gcv = p_ga_gcv < alpha, ga_loocv = p_ga_loo < alpha) ~
                          scenario, data = pvalues, FUN = sum)
  rates <- c("gcv", "loocv", "ga_gcv", "ga_loocv")
  rejected[rates] <- rejected[rates] / rejected$runs
  table <- merge(scenarios, rejected, by = "scenario", sort = FALSE)
  table[order(match(table$scenario, scenarios$scenario)), ]
}

#' Build one of Tables 1 to 3.
#'
#' Rejection rates from the p-values of the archived run
#' (archive/run0204b/table<N>_pvalues.csv) or, with `opts$rerun`, of
#' replications run here. The table is written to results/<name>.csv.
#'
#' @param name "table1_type1_error", "table2_power_fixed" or
#'   "table3_power_random".
#' @param beta_type As for `simulate_design()`.
#' @param opts Options: `rerun`, `scenarios` ("all" or an index text) and
#'   those of `run_replications()`.
#' @param settings Settings of the paper.
#' @param paths Project paths (`archive`, `results`).
#' @param design Simulation design.
#' @return The table, invisibly.
simulation_study <- function(name, beta_type, opts, settings, paths,
                             design = SIMULATION_DESIGN) {
  scenarios <- simulation_scenarios(beta_type == "fixed", design)
  if (opts$scenarios != "all") {
    scenarios <- scenarios[scenarios$scenario %in% parse_index_set(opts$scenarios), ]
  }
  pvalues <- if (opts$rerun) {
    run_dir <- file.path(paths$results, paste0(name, "_", settings_tag(opts)))
    run_replications(run_dir, beta_type, scenarios, opts, settings, design)
  } else {
    read.csv(file.path(paths$archive, "run0204b", sub("_.*", "_pvalues.csv", name)))
  }
  table <- rejection_rates(scenarios, pvalues, settings$alpha)
  write.csv(table, output_file(paths$results, paste0(name, ".csv")), row.names = FALSE)
  print(table, row.names = FALSE)
  invisible(table)
}

#' Penalties chosen by GCV and LOOCV in one replication (Figure 3).
#'
#' From the archived run (archive/run0204b/fig3_lambdas.csv) or, with
#' `opts$rerun`, simulated and chosen here (stored in results/).
#'
#' @param panels Data frame with columns `panel`, `beta_type`, `tau`.
#' @param replication Replication number.
#' @param p Number of predictors.
#' @param covariance Covariance type, as for `make_covariance()`.
#' @param opts Options: `rerun`, `perms`, `negative`, `positive`, `cores`.
#' @param settings Settings of the paper.
#' @param paths Project paths (`archive`, `results`).
#' @param design Simulation design.
#' @return Data frame with columns `panel`, `b` (1 = observed data, then the
#'   permutations), `GCV`, `LOOCV`.
chosen_penalties <- function(panels, replication, p, covariance, opts, settings, paths,
                             design = SIMULATION_DESIGN) {
  if (!opts$rerun) return(read.csv(file.path(paths$archive, "run0204b", "fig3_lambdas.csv")))
  cached_csv(output_file(paths$results, sprintf("fig3_lambdas_%s.csv", settings_tag(opts))), function() {
    simulated <- lapply(c(zero = "zero", fixed = "fixed"), function(beta_type) {
      simulate_design(design$n, p, covariance, design$rho, beta_type, design$beta_mean,
                      settings$sigma, settings$seed)
    })
    chosen <- parallel::mclapply(seq_len(nrow(panels)), function(i) {
      ada <- simulation_run(simulated[[panels$beta_type[i]]], replication,
                            design$quantiles[[panels$tau[i]]], settings$sigma, opts$negative,
                            opts$positive, opts$perms, settings$seed, settings$upper)
      data.frame(panels[i, ], b = seq_along(ada$lambda_gcv), GCV = ada$lambda_gcv,
                 LOOCV = ada$lambda_loo, row.names = NULL)
    }, mc.cores = opts$cores)
    do.call(rbind, chosen)
  })
}
