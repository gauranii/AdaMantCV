# -----------------------------------------------------------------------------
# File:    scripts/table1_type1_error.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Table 1.
# Purpose: Type I error of AdaMantCV (GCV, LOOCV and their gamma approximations)
#          in 12 simulated settings with no association (beta = 0).
# Inputs:  None (simulated data).
# Outputs: results/table1_type1_error_<settings>/scenarioSS_runRRRR.csv (p-values of each
#          replication), results/table1_type1_error_<settings>.csv (the table).
# Usage:   Rscript scripts/table1_type1_error.R [--runs=500] [--scenarios=all] [--perms=1000]
#          [--negative=50] [--positive=200] [--cores=4]
#          --runs: N (runs 1 to N), a:b or a comma list; finished runs are kept.
#          --scenarios: original scenario numbers, a:b or a comma list.
#          --method=paper: the method as the paper states it (README).
# Source:  run0204b: R_type1.R, R_driver.R (driver_type1), R_simulations.R.
# -----------------------------------------------------------------------------

source("R/setup.R")

opts <- parse_args(list(method = "original",
                        runs = "500", scenarios = "all", perms = N_PERMS,
                        negative = 50, positive = 200, cores = 4))

simulation_study("table1_type1_error", beta_type = "zero", opts)
