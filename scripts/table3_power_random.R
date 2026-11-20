# -----------------------------------------------------------------------------
# File:    scripts/table3_power_random.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Table 3.
# Purpose: Power of AdaMantCV in 12 simulated settings with random effects
#          (coefficients drawn from N(1, 0.5)).
# Inputs:  archive/run0204b/table3_pvalues.csv (p-values of the 500 replications
#          of the archived run); none with --rerun (simulated data).
# Outputs: results/Table3_power_random-effects.csv (the table); with --rerun also
#          results/Table3_power_random-effects_<settings>/scenarioSS_runRRRR.csv.
# Usage:   Rscript scripts/table3_power_random.R [--rerun] [--runs=500] [--scenarios=all]
#          [--perms=1000] [--negative=50] [--positive=200] [--cores=4]
#          --rerun: run the replications here instead of reading the archive.
#          --runs: N (runs 1 to N), a:b or a comma list; finished runs are kept.
#          --scenarios: original scenario numbers, a:b or a comma list.
# Source:  run0204b: R_power.R, R_driver.R (driver_power, scenarios 2 to 18).
# -----------------------------------------------------------------------------

source("R/setup.R")

opts <- parse_args(c(list(runs = "500", scenarios = "all"), test_options(PAPER)))

simulation_study("Table3_power_random-effects", beta_type = "random", opts, PAPER, PATHS)
