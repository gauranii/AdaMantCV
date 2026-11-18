# -----------------------------------------------------------------------------
# File:    scripts/table3_power_random.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Table 3.
# Purpose: Power of AdaMantCV in 12 simulated settings with random effects
#          (coefficients drawn from N(1, 0.5)).
# Inputs:  None (simulated data).
# Outputs: results/table3_power_random_<settings>/scenarioSS_runRRRR.csv (p-values of each
#          replication), results/table3_power_random_<settings>.csv (the table).
# Usage:   Rscript scripts/table3_power_random.R [--runs=500] [--scenarios=all] [--perms=1000]
#          [--negative=50] [--positive=200] [--cores=4]
#          --runs: N (runs 1 to N), a:b or a comma list; finished runs are kept.
#          --scenarios: original scenario numbers, a:b or a comma list.
#          --method=paper: the method as the paper states it (README).
# Source:  run0204b: R_power.R, R_driver.R (driver_power, scenarios 2 to 18).
# -----------------------------------------------------------------------------

source("R/setup.R")

opts <- parse_args(list(method = "original",
                        runs = "500", scenarios = "all", perms = N_PERMS,
                        negative = 50, positive = 200, cores = 4))

simulation_study("table3_power_random", beta_type = "random", opts)
