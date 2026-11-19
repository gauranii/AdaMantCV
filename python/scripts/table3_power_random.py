# -----------------------------------------------------------------------------
# File:    python/scripts/table3_power_random.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Table 3.
# Purpose: Power of AdaMantCV in 12 simulated settings with random effects
#          (coefficients drawn from N(1, 0.5)).
# Inputs:  archive/run0204b/table3_pvalues.csv (p-values of the 500 replications
#          of the archived run); none with --rerun (simulated data).
# Outputs: results/table3_power_random.csv (the table); with --rerun also
#          results/table3_power_random_<settings>/scenarioSS_runRRRR.csv.
# Usage:   python python/scripts/table3_power_random.py [--rerun] [--runs=500] [--scenarios=all]
#          [--perms=1000] [--negative=50] [--positive=200] [--cores=4]
#          --rerun: run the replications here instead of reading the archive.
#          --runs: N (runs 1 to N), a:b or a comma list; finished runs are kept.
#          --scenarios: original scenario numbers, a:b or a comma list.
# Source:  scripts/table3_power_random.R.
# -----------------------------------------------------------------------------
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from adamant.setup import PAPER, PATHS, parse_args, test_options  # noqa: E402
from adamant.simulation import simulation_study  # noqa: E402

opts = parse_args({"runs": "500", "scenarios": "all", **test_options(PAPER)})

simulation_study("table3_power_random", "random", opts, PAPER, PATHS)
