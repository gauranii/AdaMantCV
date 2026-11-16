# AdaMantCV: Adaptive Mantel test with cross-validation

R code for

> Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High Dimensional
> Testing with Applications to Imaging Genetics. *Frontiers in Neuroscience*
> 16:836100.

The test is `adamant_cv()` in `R/adamant_cv.R`. Every figure and table of the
paper that comes from code has one script in `scripts/`. Figure 4 is a picture
from Pluta et al. (2021) and has no code.

This branch rebuilds the code that made the paper. The original scripts (run
folders `Coh7`, `Run0204b` and `Real0210a` of 2021 and 2022) were rewritten
into shared functions, and the outputs were checked against the archived run
outputs and the paper by code (`checks/`). The check results are in
`checks/RESULTS.md`. The sections below say what agrees, what does not, and why.

## Layout

| Path | Content |
|---|---|
| `R/setup.R` | Loader sourced by every script: settings, seeds, paths, argument parser |
| `R/mantel.R` | Truncated SVD, grid of ridge penalties, Mantel statistics, permutation and gamma p-values |
| `R/ridge_cv.R` | GCV and LOOCV scores of ridge regression over the penalty grid |
| `R/adamant_cv.R` | `adamant_cv()`, the test with penalties chosen by GCV or LOOCV |
| `R/simulation.R` | Simulated designs and the simulation study of Tables 1 to 3 |
| `R/power_curve.R` | Power and GCV curves over the penalty grid (Figures 2, 5 to 8) |
| `R/data.R` | Readers for the EEG and liver data, labels of the EEG triples |
| `R/plots.R` | Heat maps and scatter plots |
| `scripts/` | One script per figure or table |
| `reference/` | `triple.csv` (the 250 channel pair and band triples) and `chan_names.csv` |
| `checks/` | Check scripts, the numbers printed in the paper (`paper_values.csv`) and `RESULTS.md` |
| `shell.nix` | Pinned R environment (R 4.6.1 with ggplot2, MASS and glmnet) |

## Data

The data are not in the repository.

- EEG coherence and SNPs of 350 subjects: one file `eeg/dat_ch<k>.csv` per
  triple k (a row id, the coherence, then 497 SNPs), as in the run folder
  `Real0210a`.
- Liver toxicity data (64 samples, 3116 genes, 10 clinical variables) from
  the mixOmics package, saved as `liver/mixOmicsData.csv`.

Scripts look for the data in `data/` or in the folder named by `ADAMANT_DATA`.

## Running

From the repository root:

```
nix-shell
export ADAMANT_DATA=/path/to/data
Rscript scripts/table4_eeg_pvalues.R
```

Results go to `results/` (csv) and `figures/` (pdf). Each script lists its
options in its header. The defaults are the settings of the paper (seed 8891,
1000 permutations, alpha 0.05, noise sd 0.5, largest penalty 1000).

| Paper | Script | Full setting | Time measured here |
|---|---|---|---|
| Figure 1 | `fig1_liver_cv.R` | as published | 1.2 min |
| Figure 2 | `fig2_liver_heat.R` | 1000 runs, 50 + 250 penalties | 6.3 min, 3 cores |
| Figure 3 | `fig3_lambda_hist.R` | replication 500, 50 + 200 penalties | 2.0 min, 6 cores |
| Figure 5 | `fig5_eeg_heat.R` | 1000 runs, 500 + 2000 penalties | not run in full |
| Figure 6 | `fig6_eeg_heat.R` | same | not run in full |
| Figure 7 | `fig7_power_gcv_decreasing.R` | same | not run in full |
| Figure 8 | `fig8_power_gcv_increasing.R` | same | not run in full |
| Table 1 | `table1_type1_error.R` | 500 runs x 12 settings | runs 1 to 5: 16 min, 5 cores |
| Table 2 | `table2_power_fixed.R` | same | runs 1 to 5 only |
| Table 3 | `table3_power_random.R` | same | runs 1 to 5 only |
| Table 4 | `table4_eeg_pvalues.R` | 5 triples | 1.4 min, 5 cores |

The machine used for these times was heavily loaded. At these rates a full
run of one of Tables 1 to 3 (6000 tests) takes more than a day, and the full
EEG power curves (26 triples, 1000 runs each) take several days. The
simulation scripts write one file per replication and skip the ones already
written, so a full run can be split, for example
`Rscript scripts/table1_type1_error.R --runs=1:100 --cores=8`. The EEG curves
are cached in `results/` and shared by Figures 5 to 8; `--runs=10` gives a
quick, noisy version.

## Checks

`bash checks/run_checks.sh`, with `ADAMANT_ARCHIVE` set to a copy of the
archived run outputs and `ADAMANT_PUBLISHED` to the published figure files,
writes `checks/RESULTS.md`. It checks four things.

1. The archived run outputs give the numbers printed in Tables 1 to 4.
2. The rebuilt functions give the same outputs as the original functions on
   the same inputs and seeds (`checks/check_equivalence.R`, small permutation
   counts; needs the original scripts and dplyr).
3. Each published figure equals the same figure redrawn by the rebuilt scripts
   from the archived data behind it: heat map tile colours (Figures 2, 5, 6),
   histogram bar heights (Figure 3), point positions (Figures 1, 7, 8).
4. The rebuilt results equal the archived results where the same computation
   could be repeated here.

## What was reproduced, and what was not

- **Figures 2, 3, 5 to 8, drawn from the archived data**: every tile colour,
  bar and point of the published files is reproduced.
- **Figure 1A**: the ten minimum points match the published figure.
  **Figure 1B**: the original code set no seed before the cross-validation
  folds, so the published panel cannot be reproduced fold for fold; the script
  fixes a seed.
- **Figure 2, recomputed in full**: equal to the archived run to about 1e-13,
  except at lambda = 0 (note 1).
- **Figures 5 to 8, recomputed**: the GCV and LOOCV curves equal the archive
  (largest relative difference about 1e-7, at penalties where the scores reach
  about 1e10). The power curves were run with far fewer runs than the paper.
- **Table 4, recomputed in full**: the rebuilt code gives exactly what the
  original code gives on this machine, and the penalties chosen for every
  permutation equal the archive, but the p-values differ from the paper
  (note 1).
- **Tables 1 to 3 and Figure 3, recomputed**: only the heteroskedastic
  settings can be repeated exactly here (note 2), and only a few replications
  were run.

Notes, each found by the checks:

1. **Rounding decides the statistic at lambda = 0.** With n < p, U from the
   SVD of X is square, so at lambda = 0 the Mantel statistic equals the
   squared norm of Y for every permutation, and its rank among the
   permutations is set by rounding. In the EEG data the response is centred,
   its ridge fit at lambda = 0 is exact (GCV score near 0), and GCV and LOOCV
   choose lambda = 0 for the observed data in all 250 triples of `Real0210a`.
   The Table 4 p-values therefore depend on rounding: the original code run
   here gives 0.790, 0.785, 0.720, 0.719 for Fp1 - T8 (paper: 0.016, 0.015,
   0.012, 0.012), with the same chosen penalties as the archive. The same tie
   affects the lambda = 0 point of the power curves. The rebuilt code keeps
   the original behaviour.
2. **Simulated designs depend on the LAPACK build.** `MASS::mvrnorm` uses an
   eigen decomposition. The compound symmetric covariance has p - 1 equal
   eigenvalues, so its eigenvectors, and the simulated X, differ between
   builds. The diagonal (heteroskedastic) covariance has distinct eigenvalues,
   and its replications match the archive exactly.
3. **Rejections at p = 0.05.** The original runs used dplyr 0.8.5, whose
   compiled `percent_rank()` adds the steps one at a time, so a p-value of
   0.050 came out slightly above 0.05 and was not counted. `adamant_cv()`
   reproduces that arithmetic for the final p-value.
4. **Captions of Tables 2 and 3 and Figure 3.** The runs that made the
   numbers used noise sd 0.5 (`R_initial.R`) and coefficients equal to 1
   (Table 2, Figure 3C and 3D) or drawn from N(1, 0.5^2) (Table 3), as passed
   by `driver_power` in `R_driver.R`. The captions say sd 1, xi = 3 and
   N(0, sigma_b^2). The scripts follow the runs.
5. **Figure 2B** shows the power of the gamma approximation (column
   `power_MOM`), as drawn by `R_Heat3.R` under the published file name.
6. **Figure 6** was drawn from 12 triples (181, 231, 200, 144, 159, 136, 223,
   205, 25, 179, 219, 123) with colour midpoint 0.4 in panel B. The later
   `R_Heat2.R` in `Coh7` lists other triples. The list and midpoint were
   recovered by matching tile colours.
7. **Figure 3** uses the 1000 permutations of replication 500 and leaves out
   the observed data, as in the published figure.
8. **Figures 7 and 8, GCV panels**: the windows shown were set by hand and
   not saved. The scripts hold windows recovered from the published files.
9. Kept from the original code, with a comment at each place: the GCV
   residual of every permutation uses the observed response; at lambda = 0
   with n <= p the smoother is taken as I minus the inverse; the 2021 power
   code permutes the rows of U rather than the response.

## Origin of each script

Each script header names the original files it replaces. The original run
folders are on the external drive (`DUMP 010824/Active/Runs/`).
