# -----------------------------------------------------------------------------
# File:    scripts/table1_type1_error.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Table 1.
# Purpose: Type I error of AdaMantCV (GCV, LOOCV and their gamma approximations)
#          in 12 simulated settings with no association (beta = 0).
# Inputs:  archive/run0204b/table1_pvalues.csv (p-values of the 500 replications
#          of the archived run); none with --rerun (simulated data).
# Outputs: results/table1_type1_error.csv (the table); with --rerun also
#          results/table1_type1_error_<settings>/scenarioSS_runRRRR.csv.
# Usage:   Rscript scripts/table1_type1_error.R [--rerun] [--runs=500] [--scenarios=all]
#          [--perms=1000] [--negative=50] [--positive=200] [--cores=4]
#          --rerun: run the replications here instead of reading the archive.
#          --runs: N (runs 1 to N), a:b or a comma list; finished runs are kept.
#          --scenarios: original scenario numbers, a:b or a comma list.
# Source:  run0204b: R_type1.R, R_driver.R (driver_type1), R_simulations.R.
# -----------------------------------------------------------------------------

source("R/setup.R")

opts <- parse_args(list(runs = "500", scenarios = "all", perms = N_PERMS,
                        negative = 50, positive = 200, cores = 4))

simulation_study("table1_type1_error", beta_type = "zero", opts)
