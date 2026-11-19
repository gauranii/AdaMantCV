# -----------------------------------------------------------------------------
# File:    checks/write_python.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Write the comparison of the R and the Python implementation as
#          Markdown (checks/PYTHON.md) from the tables of
#          checks/check_python.py.
# Usage:   python3 checks/write_python.py <python_archive.csv>
#          <python_rerun.csv> > checks/PYTHON.md
#          A table that does not exist is left out.
# -----------------------------------------------------------------------------
import csv
import os
import sys

NAMES = {"table1": "Table 1", "table2": "Table 2", "table3": "Table 3", "table4": "Table 4",
         "fig1": "Figure 1", "fig2": "Figure 2", "fig3": "Figure 3", "fig5": "Figure 5",
         "fig6": "Figure 6", "fig7": "Figure 7", "fig8": "Figure 8"}
INTRO = {
    "archive": """## From the archived runs

Both implementations run with the defaults of their scripts, as
`nix flake check` runs them. Figure 1 is computed from the data in both.
Tolerance: none for the tables and csv files, 1e-9 (relative) for the plotted
numbers.
""",
    "rerun": """## Recomputed from the data

Both implementations run with `--rerun` (`bash checks/check_python.sh rerun
<work dir>` with `ADAMANT_DATA` set, 4 cores, one BLAS thread). Table 4 and
Figures 1, 2, 3, 7 and 8 are run in full. Tables 1 to 3 are run with 2
replications per setting and Figures 5 and 6 with 3 replications per triple,
with the seeds of the full runs. Tolerance: 1e-8 (relative). Compared are the
tables, the stored parts of the rerun (curves, chosen penalties, p-values of
every replication) and the plotted numbers.
""",
}


def section(file):
    rows = list(csv.DictReader(open(file)))
    mode = rows[0]["mode"]
    good = sum(row["verdict"] != "DIFFERENT" for row in rows)
    identical = sum(row["verdict"] == "identical" for row in rows)
    lines = [INTRO[mode], f"{good} of {len(rows)} outputs equal, {identical} of them identical "
             "in every value.", "",
             "| Output | Run | Files | Values | Not identical | Largest difference | "
             "Largest relative difference | Verdict | R (s) | Python (s) |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for row in rows:
        lines.append(f"| {NAMES[row['output']]} | {row['settings']} | {row['files']} | "
                     f"{row['values']} | {row['not_identical']} | {row['max_abs_diff']} | "
                     f"{row['max_rel_diff']} | {row['verdict']} | {row['r_seconds']} | "
                     f"{row['python_seconds']} |")
    problems = [f"- {NAMES[row['output']]}: {row['problems']}" for row in rows if row["problems"]]
    if problems:
        lines += ["", "Differences:", ""] + problems
    return "\n".join(lines) + "\n"


print("""# R and Python

Written by `checks/check_python.sh`: the outputs of the R scripts (`scripts/`)
and of the Python scripts (`python/scripts/`) compared by code
(`checks/check_python.py`). For a table, every cell of its csv file is
compared. For a figure, the numbers it plots are compared mark by mark
(tiles, bars, points, lines), written by both implementations as csv.
"Values" counts the compared cells; "identical" means every one is the same
double as read from the files, which both implementations write with 15
significant digits. Times are wall-clock seconds on an 8 core Apple M2
laptop; the Python times include compiling its numba functions on first use.

The tests in `python/tests/` (also run by `nix flake check`) compare with R,
bit for bit, the random numbers of `python/adamant/rng.py`, the gamma
probabilities of `rmath.py` and the cross-validated ridge regression of
`glmnet.py`.
""")
for file in sys.argv[1:3]:
    if os.path.exists(file):
        print(section(file))
