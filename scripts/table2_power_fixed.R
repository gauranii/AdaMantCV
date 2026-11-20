# -----------------------------------------------------------------------------
# File:    scripts/table2_power_fixed.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Table 2.
# Purpose: Power of AdaMantCV in 12 simulated settings with fixed effects
#          (every coefficient equal to 1).
# Inputs:  archive/run0204b/table2_pvalues.csv (p-values of the 500 replications
#          of the archived run); none with --rerun (simulated data).
# Outputs: results/Table2_power_fixed-effects.csv (the table); with --rerun also
#          results/Table2_power_fixed-effects_<settings>/scenarioSS_runRRRR.csv.
# Usage:   Rscript scripts/table2_power_fixed.R [--rerun] [--runs=500] [--scenarios=all]
#          [--perms=1000] [--negative=50] [--positive=200] [--cores=4]
#          --rerun: run the replications here instead of reading the archive.
#          --runs: N (runs 1 to N), a:b or a comma list; finished runs are kept.
#          --scenarios: original scenario numbers, a:b or a comma list.
# Source:  run0204b: R_power.R, R_driver.R (driver_power, scenarios 20 to 36).
# -----------------------------------------------------------------------------

source("R/setup.R")

opts <- parse_args(c(list(runs = "500", scenarios = "all"), test_options(PAPER)))

simulation_study("Table2_power_fixed-effects", beta_type = "fixed", opts, PAPER, PATHS)
