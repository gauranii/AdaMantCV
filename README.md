# AdaMantCV: Adaptive Mantel test with cross-validation

R and Python code for

> Gauran, I. I., Xue, G., Chen, C., Ombao, H., and Yu, Z. (2022). Ridge
> Penalization in High Dimensional Testing with Applications to Imaging
> Genetics. *Frontiers in Neuroscience* 16:836100.

The test is `adamant_cv()` in `R/adamant_cv.R`: a ridge-penalized Mantel test
of association between a response and a high dimensional set of predictors,
with the penalty chosen by generalized or leave-one-out cross-validation for
the observed data and for every permutation. Each table and figure of the
paper that comes from code has one script in `scripts/`, and `nix flake check`
builds all of them and compares them with the paper. Figure 4 of the paper is
a diagram and has no code. `python/` holds the same code in Python; both
implementations give the same outputs, checked by `nix flake check`.

## Layout

| Path | Content |
|---|---|
| `R/setup.R` | Loader sourced by every script: settings of the paper, paths, argument parser |
| `R/mantel.R` | Truncated SVD, grid of ridge penalties, Mantel statistics, permutation and gamma p-values |
| `R/ridge_cv.R` | GCV and LOOCV scores of ridge regression over the penalty grid |
| `R/adamant_cv.R` | `adamant_cv()`, the test with penalties chosen by GCV or LOOCV |
| `R/simulation.R` | Simulated designs and the simulation study of Tables 1 to 3 |
| `R/power_curve.R` | Power and GCV curves over the penalty grid (Figures 2, 5 to 8) |
| `R/data.R` | Readers for the data and for the archived curves, labels of the EEG triples |
| `R/eeg_tests.R` | The tests of the EEG triples (Table 4) |
| `R/liver_cv.R` | Cross-validated ridge regression of the liver data (Figure 1) |
| `R/plots.R` | Heat maps, scatter plots and histograms |
| `scripts/` | One script per table or figure |
| `python/adamant/` | The same modules in Python (`mantel.py`, `ridge_cv.py`, `adamant_cv.py`, ...), with `rng.py` (the random numbers of R), `rbase.py` and `rmath.py` (R's arithmetic, matrix products and decompositions, `pgamma`) and `glmnet.py` (cross-validated ridge regression as glmnet computes it) |
| `python/scripts/` | One Python script per table or figure, with the names and options of `scripts/` |
| `python/tests/` | Comparisons of `rng.py`, `rmath.py` and `glmnet.py` with R |
| `archive/` | Outputs of the original runs that the scripts read, the published figure files (`figures/`), and `reduce_archive.R`, which made the former from the run folders |
| `data/liver/` | The liver toxicity data |
| `reference/` | The 250 EEG triples (`triple.csv`), the channel names, and `table4_cluster_lambda0.csv` |
| `checks/` | The numbers printed in the paper (`paper_values.csv`), the comparison scripts and their result (`RESULTS.md`), and `check_equivalence.R`, which compares the outputs of two versions of the code; `check_python.sh`, which compares the outputs of the R and the Python scripts, and its result (`PYTHON.md`) |
| `flake.nix`, `flake.lock` | Pinned environment (R with ggplot2, MASS and glmnet; Python with numpy, scipy, pandas, matplotlib and numba), one command per output, and the checks |

## Running

With [Nix](https://nixos.org) (flakes enabled), from the repository root:

```sh
nix run .#table4        # one output: table1 to table4, fig1 to fig3, fig5 to fig8
nix run .#all           # every output
nix run .#py-table4     # the same output with the Python scripts (py-table1, ..., py-fig8, py-all)
nix flake check         # build every output, compare it with the paper, compare R with Python
```

Tables are written to `results/` (csv) and figures to `figures/` (pdf).
`nix flake check` compares every cell of Tables 1 to 4 with the paper and
every figure with the published file; `checks/RESULTS.md` holds the comparison.
It also builds every output with the Python scripts and compares it with the
output of the R scripts, cell by cell for the tables and mark by mark for the
numbers the figures plot; `checks/PYTHON.md` holds that comparison.

To recompute an output from the data instead of the archived runs, add
`--rerun`. The EEG outputs then need the EEG data (see Data):

```sh
export ADAMANT_DATA=/path/to/data      # the folder that holds eeg/
nix run .#table4 -- --rerun --cores=5
```

`nix develop` opens a shell with the same R and Python, where the scripts can
be run directly: `Rscript scripts/table4_eeg_pvalues.R [--rerun]` or
`python python/scripts/table4_eeg_pvalues.py [--rerun]`. Each script lists
its options in its header. The settings are those of the paper: seed 8891, 1000
permutations, significance level 0.05.

## Outputs

| Paper | Script | Command | From archive | With `--rerun` |
|---|---|---|---|---|
| Table 1 | `table1_type1_error.R` | `nix run .#table1` | 1 s | about 9 h, 8 cores (estimate) |
| Table 2 | `table2_power_fixed.R` | `nix run .#table2` | 1 s | about 9 h, 8 cores (estimate) |
| Table 3 | `table3_power_random.R` | `nix run .#table3` | 1 s | about 9 h, 8 cores (estimate) |
| Table 4 | `table4_eeg_pvalues.R` | `nix run .#table4` | 1 s | 35 s, 5 cores |
| Figure 1 | `fig1_liver_cv.R` | `nix run .#fig1` | 20 s (always computed from the data) | same |
| Figure 2 | `fig2_liver_heat.R` | `nix run .#fig2` | 1 s | 2 min, 7 cores |
| Figure 3 | `fig3_lambda_hist.R` | `nix run .#fig3` | 1 s | 1 min, 4 cores |
| Figure 5 | `fig5_eeg_heat.R` | `nix run .#fig5` | 1 s | about 55 min, 7 cores (estimate) |
| Figure 6 | `fig6_eeg_heat.R` | `nix run .#fig6` | 1 s | about 55 min, 7 cores (estimate) |
| Figure 7 | `fig7_power_gcv_decreasing.R` | `nix run .#fig7` | 1 s | 14 min, 7 cores |
| Figure 8 | `fig8_power_gcv_increasing.R` | `nix run .#fig8` | 1 s | 14 min, 7 cores |

Times are from an 8 core Apple M2 laptop; estimates are scaled from 3
replications (Tables 1 to 3) or 3 triples (Figures 5 and 6). A rerun stores its parts in
`results/` and continues from them: Figures 5 to 8 share the curves of their
triples, and the simulation scripts write one file per replication, so a
table can be run in pieces, for example
`nix run .#table1 -- --rerun --runs=1:100 --cores=8`.

Tables 1 to 4 and Figures 2, 3, 5 to 8 are computed from the archived runs of
the original analysis (`archive/`), run on the Ibex cluster and a KAUST
laptop; `--rerun` recomputes them from the data. Table 4 uses the observed
test statistics saved by the original run
(`reference/table4_cluster_lambda0.csv`). On a different machine or R build,
recomputed results can differ.

## Data

- **Liver toxicity data** (`data/liver/mixOmicsData.csv`): 64 samples, 10
  clinical variables and 3116 gene expressions; the public `liver.toxicity`
  data of the R package [mixOmics](https://mixomics.org), exported to csv.
  Used by Figures 1 and 2.
- **EEG coherence and SNPs** of 350 subjects: one file `eeg/dat_ch<k>.csv`
  per triple k of channel pair and frequency band (a row id, the coherence,
  then 497 SNPs). These subject-level data are not in the repository and are
  needed only for `--rerun` of Table 4 and Figures 5 to 8; the scripts look
  for them in the folder named by `ADAMANT_DATA`.
- **Archived runs** (`archive/`): per-replication p-values of Tables 1 to 3,
  chosen penalties of Figure 3, per-permutation p-values of Table 4, and the
  power and cross-validation curves of Figures 2 and 5 to 8. They hold no
  subject-level data.
