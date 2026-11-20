# -----------------------------------------------------------------------------
# File:    python/scripts/table1_type1_error.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Table 1.
# Purpose: Type I error of AdaMantCV (GCV, LOOCV and their gamma approximations)
#          in 12 simulated settings with no association (beta = 0).
# Inputs:  archive/run0204b/table1_pvalues.csv (p-values of the 500 replications
#          of the archived run); none with --rerun (simulated data).
# Outputs: results/Table1_type1-error.csv (the table); with --rerun also
#          results/Table1_type1-error_<settings>/scenarioSS_runRRRR.csv.
# Usage:   python python/scripts/table1_type1_error.py [--rerun] [--runs=500] [--scenarios=all]
#          [--perms=1000] [--negative=50] [--positive=200] [--cores=4]
#          --rerun: run the replications here instead of reading the archive.
#          --runs: N (runs 1 to N), a:b or a comma list; finished runs are kept.
#          --scenarios: original scenario numbers, a:b or a comma list.
# Source:  scripts/table1_type1_error.R.
# -----------------------------------------------------------------------------
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from adamant.setup import PAPER, PATHS, parse_args, test_options  # noqa: E402
from adamant.simulation import simulation_study  # noqa: E402

opts = parse_args({"runs": "500", "scenarios": "all", **test_options(PAPER)})

simulation_study("Table1_type1-error", "zero", opts, PAPER, PATHS)
