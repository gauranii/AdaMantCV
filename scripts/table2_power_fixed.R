# -----------------------------------------------------------------------------
# File:    scripts/table2_power_fixed.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Table 2.
# Purpose: Power of AdaMantCV in 12 simulated settings with fixed effects
#          (every coefficient equal to 1).
# Inputs:  None (simulated data).
# Outputs: results/table2_power_fixed_<settings>/scenarioSS_runRRRR.csv (p-values of each
#          replication), results/table2_power_fixed_<settings>.csv (the table).
# Usage:   Rscript scripts/table2_power_fixed.R [--runs=500] [--scenarios=all] [--perms=1000]
#          [--negative=50] [--positive=200] [--cores=4]
#          --runs: N (runs 1 to N), a:b or a comma list; finished runs are kept.
#          --scenarios: original scenario numbers, a:b or a comma list.
# Source:  run0204b: R_power.R, R_driver.R (driver_power, scenarios 20 to 36).
# -----------------------------------------------------------------------------

source("R/setup.R")

opts <- parse_args(list(runs = "500", scenarios = "all", perms = N_PERMS,
                        negative = 50, positive = 200, cores = 4))

simulation_study("table2_power_fixed", beta_type = "fixed", opts)
