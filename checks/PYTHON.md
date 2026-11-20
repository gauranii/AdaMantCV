# R and Python

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

## From the archived runs

Both implementations run with the defaults of their scripts, as
`nix flake check` runs them. Figure 1 is computed from the data in both.
Tolerance: none for the tables and csv files, 1e-9 (relative) for the plotted
numbers.

11 of 11 outputs equal, 11 of them identical in every value.

| Output | Run | Files | Values | Not identical | Largest difference | Largest relative difference | Verdict | R (s) | Python (s) |
|---|---|---|---|---|---|---|---|---|---|
| Table 1 | archived runs | 1 | 108 | 0 | 0 | 0 | identical | 1.0 | 1.0 |
| Table 2 | archived runs | 1 | 108 | 0 | 0 | 0 | identical | 1.0 | 1.0 |
| Table 3 | archived runs | 1 | 108 | 0 | 0 | 0 | identical | 1.0 | 1.0 |
| Table 4 | archived runs | 1 | 50 | 0 | 0 | 0 | identical | 1.0 | 1.0 |
| Figure 1 | from the data | 3 | 598040 | 0 | 0 | 0 | identical | 18.5 | 8.2 |
| Figure 2 | archived runs | 2 | 48016 | 0 | 0 | 0 | identical | 1.0 | 2.1 |
| Figure 3 | archived runs | 4 | 2196 | 0 | 0 | 0 | identical | 2.1 | 2.1 |
| Figure 5 | archived runs | 2 | 480016 | 0 | 0 | 0 | identical | 2.1 | 2.1 |
| Figure 6 | archived runs | 2 | 480016 | 0 | 0 | 0 | identical | 2.1 | 2.1 |
| Figure 7 | archived runs | 6 | 73935 | 0 | 0 | 0 | identical | 2.1 | 2.1 |
| Figure 8 | archived runs | 6 | 71955 | 0 | 0 | 0 | identical | 2.1 | 2.1 |

## Recomputed from the data

Both implementations run with `--rerun` (`bash checks/check_python.sh rerun
<work dir>` with `ADAMANT_DATA` set, 4 cores, one BLAS thread). Table 4 and
Figures 1, 2, 3, 7 and 8 are run in full. Tables 1 to 3 are run with 2
replications per setting and Figures 5 and 6 with 3 replications per triple,
with the seeds of the full runs. Tolerance: 1e-8 (relative). Compared are the
tables, the stored parts of the rerun (curves, chosen penalties, p-values of
every replication) and the plotted numbers. A run whose processes pass 3 GB
of memory is stopped and run again with fewer cores (shown beside its time);
the number of cores does not change a result.

11 of 11 outputs equal, 11 of them identical in every value.

| Output | Run | Files | Values | Not identical | Largest difference | Largest relative difference | Verdict | R (s) | Python (s) |
|---|---|---|---|---|---|---|---|---|---|
| Table 1 | 2 replications | 25 | 252 | 0 | 0 | 0 | identical | 127.5 | 135.5 |
| Table 2 | 2 replications | 25 | 252 | 0 | 0 | 0 | identical | 228.3 | 148.0 |
| Table 3 | 2 replications | 25 | 252 | 0 | 0 | 0 | identical | 275.1 | 174.6 |
| Table 4 | full | 1 | 50 | 0 | 0 | 0 | identical | 57.2 | 40.3 |
| Figure 1 | full | 3 | 598040 | 0 | 0 | 0 | identical | 17.5 | 8.2 |
| Figure 2 | full | 22 | 84016 | 0 | 0 | 0 | identical | 160.4 | 136.7 |
| Figure 3 | full | 5 | 26472 | 0 | 0 | 0 | identical | 53.3 | 35.3 |
| Figure 5 | 3 replications | 26 | 840016 | 0 | 0 | 0 | identical | 963.4 | 925.8 |
| Figure 6 | 3 replications | 26 | 840016 | 0 | 0 | 0 | identical | 950.3 | 914.2 |
| Figure 7 | full | 12 | 163935 | 0 | 0 | 0 | identical | 423.6 | 332.8 |
| Figure 8 | full | 12 | 161955 | 0 | 0 | 0 | identical | 641.5 (2 cores) | 361.8 |

