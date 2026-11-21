# -----------------------------------------------------------------------------
# File:    tests/testthat/test-simulation.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Unit tests of R/simulation.R: the covariance builders, the simulated designs, the
#          scenario numbers, one replication, the stored replications of a
#          small rerun, the rejection rates of a hand-made table, and the
#          tables read from the archived runs.
# Inputs:  archive/run0204b/ (p-values and chosen penalties of the archived runs);
#          files are written to temporary folders.
# Outputs: None.
# Usage:   Rscript tests/testthat.R      # all test files
# -----------------------------------------------------------------------------

source("../../R/setup.R")

# A design with 12 observations and 4 predictors, in place of
# SIMULATION_DESIGN.
SMALL_DESIGN <- list(n = 12, rho = 0.025, beta_mean = 1,
                     quantiles = c(none = 0, Q1 = 0.25, Q2 = 0.5),
                     covariances = c(CS = "CSYM", Heteroskedastic = "NCDIAG"),
                     dimensions = 4)

test_that("make_covariance builds the two structures", {
  expect_equal(make_covariance("CSYM", 3, 0.025),
               rbind(c(1, 0.025, 0.025), c(0.025, 1, 0.025), c(0.025, 0.025, 1)))
  expect_equal(make_covariance("NCDIAG", 3, 0.025), diag(log(c(2, 3, 4))))
  expect_equal(dim(make_covariance("NCDIAG", 1, 0)), c(1, 1))
  expect_error(make_covariance("AR1", 3, 0.5), "unknown covariance type AR1")
})

test_that("simulate_design returns a standardized design and the coefficients", {
  fixed <- simulate_design(30, 4, "CSYM", 0.025, "fixed", mu = 1, sigma = 0.5, seed = 8891)
  expect_equal(dim(fixed$X), c(30, 4))
  expect_equal(unname(colMeans(fixed$X)), rep(0, 4))
  expect_equal(unname(apply(fixed$X, 2, sd)), rep(1, 4))
  expect_equal(fixed$beta, rep(1, 4))
  expect_equal(simulate_design(30, 4, "CSYM", 0.025, "zero", 1, 0.5, 8891)$beta, rep(0, 4))
  expect_error(simulate_design(30, 4, "CSYM", 0.025, "other", 1, 0.5, 8891),
               "unknown beta type other")
})

test_that("simulate_design depends on its seed only", {
  fixed <- simulate_design(30, 4, "CSYM", 0.025, "fixed", 1, 0.5, seed = 8891)
  runif(5)
  random <- simulate_design(30, 4, "CSYM", 0.025, "random", 1, 0.5, seed = 8891)
  expect_identical(random$X, fixed$X)
  expect_false(identical(simulate_design(30, 4, "CSYM", 0.025, "fixed", 1, 0.5, 8892)$X, fixed$X))
  # the random coefficients are the draws that follow the design in the stream
  set.seed(8891)
  MASS::mvrnorm(30, rep(0, 4), Sigma = make_covariance("CSYM", 4, 0.025))
  expect_equal(random$beta, rnorm(4, 1, 0.5))
})

test_that("simulation_scenarios numbers the 12 settings as the original runs", {
  null <- simulation_scenarios(FALSE, SIMULATION_DESIGN)
  expect_named(null, c("scenario", "covariance", "p", "tau"))
  expect_equal(null$scenario, c(2, 5, 8, 11, 14, 17, 3, 6, 9, 12, 15, 18))
  expect_equal(null$covariance, rep(c("CS", "Heteroskedastic"), each = 6))
  expect_equal(null$p, rep(rep(c(500, 1000), each = 3), 2))
  expect_equal(null$tau, rep(c("none", "Q1", "Q2"), 4))
  expect_equal(simulation_scenarios(TRUE, SIMULATION_DESIGN)$scenario, null$scenario + 18)
})

test_that("simulation_run draws the response with seed 10 * run + seed", {
  simulated <- simulate_design(20, 5, "CSYM", 0.025, "fixed", 1, 0.5, seed = 8891)
  ada <- simulation_run(simulated, run = 3, q = 0.25, sigma = 0.5, n_negative = 3,
                        n_positive = 3, n_perms = 9, seed = 8891)
  set.seed(10 * 3 + 8891)
  Y <- simulated$X %*% simulated$beta + rnorm(20, mean = 0, sd = 0.5)
  expect_identical(ada, adamant_cv(simulated$X, Y, 0.25, 3, 3, 9, 8891))
  expect_named(ada$p_values, TESTS)
  other <- simulation_run(simulated, run = 4, 0.25, 0.5, 3, 3, 9, seed = 8891)
  expect_false(identical(other$statistics, ada$statistics))
})

test_that("simulate_designs makes one design per covariance and dimension", {
  design <- modifyList(SMALL_DESIGN, list(dimensions = c(4, 6)))
  scenarios <- simulation_scenarios(FALSE, design)
  simulated <- simulate_designs(scenarios, "zero", SMALL_SETTINGS, design)
  expect_named(simulated, c("CS 4", "CS 6", "Heteroskedastic 4", "Heteroskedastic 6"))
  expect_equal(dim(simulated[["Heteroskedastic 6"]]$X), c(12, 6))
  expect_identical(simulated[["CS 4"]],
                   simulate_design(12, 4, "CSYM", 0.025, "zero", 1, 0.5, 8891))
})

test_that("run_replications stores one file per replication and reuses them", {
  run_dir <- file.path(new_folder(), "runs")
  scenarios <- simulation_scenarios(FALSE, SMALL_DESIGN)[1:2, ]
  opts <- modifyList(SMALL_OPTS, list(runs = "2"))
  expect_output(pvalues <- run_replications(run_dir, "zero", scenarios, opts, SMALL_SETTINGS,
                                            SMALL_DESIGN), "4 replications to run")
  expect_named(pvalues, c("scenario", "run", "p_gcv", "p_loo", "p_ga_gcv", "p_ga_loo"))
  expect_equal(pvalues$scenario, c(2, 2, 5, 5))
  expect_equal(pvalues$run, c(1, 2, 1, 2))
  expect_equal(sort(list.files(run_dir)),
               c("scenario02_run0001.csv", "scenario02_run0002.csv",
                 "scenario05_run0001.csv", "scenario05_run0002.csv"))
  expect_true(all(unlist(pvalues[3:6]) >= 0 & unlist(pvalues[3:6]) <= 1))

  # a replication is the test of simulation_run() on the design of its setting
  simulated <- simulate_design(12, 4, "CSYM", 0.025, "zero", 1, 0.5, 8891)
  ada <- simulation_run(simulated, run = 2, q = 0.25, 0.5, 3, 3, 9, 8891)
  expect_equal(unlist(pvalues[pvalues$scenario == 5 & pvalues$run == 2, 3:6]),
               ada$p_values, check.attributes = FALSE)

  # stored replications are not run again; new ones are added
  opts$runs <- "3"
  expect_output(more <- run_replications(run_dir, "zero", scenarios, opts, SMALL_SETTINGS,
                                         SMALL_DESIGN), "2 replications to run")
  expect_equal(nrow(more), 6)
  expect_equal(more[more$run <= 2, ], pvalues, check.attributes = FALSE)
})

test_that("rejection_rates counts the p-values below alpha per setting", {
  scenarios <- data.frame(scenario = c(5, 2), covariance = "CS", p = c(4, 4),
                          tau = c("Q1", "none"))
  pvalues <- data.frame(scenario = c(2, 2, 2, 2, 5, 5),
                        run = c(1, 2, 3, 4, 1, 2),
                        p_gcv = c(0.01, 0.04, 0.05, 0.9, 0.5, 0.6),
                        p_loo = c(0.01, 0.2, 0.3, 0.9, 0.01, 0.6),
                        p_ga_gcv = c(0.01, 0.02, 0.03, 0.04, 0.5, 0.6),
                        p_ga_loo = c(0.5, 0.5, 0.5, 0.5, 0.01, 0.02))
  rates <- rejection_rates(scenarios, pvalues, alpha = 0.05)
  expect_equal(rates$scenario, c(5, 2))        # the order of `scenarios`
  expect_equal(rates$runs, c(2, 4))
  expect_equal(rates$gcv, c(0, 0.5))           # 0.05 itself is not below alpha
  expect_equal(rates$loocv, c(0.5, 0.25))
  expect_equal(rates$ga_gcv, c(0, 1))
  expect_equal(rates$ga_loocv, c(1, 0))
  expect_equal(rates$tau, c("Q1", "none"))
})

test_that("simulation_study builds Table 1 from the archived p-values", {
  paths <- project_paths(REPO_ROOT, output = new_folder())
  opts <- list(rerun = FALSE, scenarios = "all")
  table <- quietly(simulation_study("Table1_type1-error", "zero", opts, PAPER, paths))
  expect_equal(dim(table), c(12, 9))
  expect_equal(table$scenario, simulation_scenarios(FALSE, SIMULATION_DESIGN)$scenario)
  expect_true(all(table$runs == 500))
  # under the null every rejection rate is near the level 0.05
  rates <- unlist(table[c("gcv", "loocv", "ga_gcv", "ga_loocv")])
  expect_true(all(rates >= 0.01 & rates <= 0.1))
  expect_equal(read.csv(file.path(paths$results, "Table1_type1-error.csv")), table,
               check.attributes = FALSE)

  opts$scenarios <- "2,17"
  expect_equal(quietly(simulation_study("Table1_type1-error", "zero", opts, PAPER, paths))$scenario,
               c(2, 17))
})

test_that("simulation_study reruns a small design", {
  paths <- project_paths(REPO_ROOT, output = new_folder())
  opts <- modifyList(SMALL_OPTS, list(runs = "2", scenarios = "20,23"))
  table <- quietly(simulation_study("Table2_power_fixed-effects", "fixed", opts, SMALL_SETTINGS,
                                    paths, design = SMALL_DESIGN))
  expect_equal(table$scenario, c(20, 23))
  expect_equal(table$runs, c(2, 2))
  expect_true(dir.exists(file.path(paths$results,
                                   "Table2_power_fixed-effects_perms9_grid3x3")))
})

test_that("chosen_penalties reads the archive or chooses the penalties again", {
  panels <- data.frame(panel = c("a", "b"), beta_type = c("zero", "fixed"), tau = c("Q1", "Q2"))
  paths <- project_paths(REPO_ROOT, output = new_folder())
  archived <- chosen_penalties(panels, 1, 500, "CSYM", list(rerun = FALSE), PAPER, paths)
  expect_named(archived, c("panel", "b", "GCV", "LOOCV"))
  expect_equal(unique(archived$panel), c("a", "b", "c", "d"))
  expect_equal(as.vector(table(archived$panel)), rep(1001, 4))

  chosen <- chosen_penalties(panels, replication = 1, p = 5, covariance = "CSYM", SMALL_OPTS,
                             SMALL_SETTINGS, paths, design = SMALL_DESIGN)
  expect_named(chosen, c("panel", "beta_type", "tau", "b", "GCV", "LOOCV"))
  expect_equal(chosen$panel, rep(c("a", "b"), each = 10))
  expect_equal(chosen$b, rep(1:10, 2))
  simulated <- simulate_design(12, 5, "CSYM", 0.025, "fixed", 1, 0.5, 8891)
  ada <- simulation_run(simulated, 1, q = 0.5, 0.5, 3, 3, 9, 8891)
  expect_equal(chosen$GCV[chosen$panel == "b"], ada$lambda_gcv)
  expect_true(file.exists(file.path(paths$results, "fig3_lambdas_perms9_grid3x3.csv")))
})
