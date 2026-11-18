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
| `reference/` | `triple.csv` (the 250 channel pair and band triples), `chan_names.csv`, and `table4_cluster_lambda0.csv` (note 1) |
| `checks/` | Check scripts, the numbers printed in the paper (`paper_values.csv`) and `RESULTS.md` |
| `checks/search/` | The searches for the best fit (Figure 1B folds, chosen penalties of the archived runs) |
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
1000 permutations, alpha 0.05, noise sd 0.5, largest penalty 1000) and the
code that best fits the published numbers. The method as the paper states it
is an option: `--method=paper` for Tables 1 to 3 and Figures 2, 3 and 5 to 8,
`--version=paper` for Table 4. Its outputs carry `paper` in their file names,
so the two versions never overwrite each other. See "Best fit and correct
version, side by side" below.

| Paper | Script | Full setting | Time measured here |
|---|---|---|---|
| Figure 1 | `fig1_liver_cv.R` | as published | 1.2 min |
| Figure 2 | `fig2_liver_heat.R` | 1000 runs, 50 + 250 penalties | 6.3 min, 3 cores |
| Figure 3 | `fig3_lambda_hist.R` | replication 500, 50 + 200 penalties | 2.0 min, 6 cores |
| Figure 5 | `fig5_eeg_heat.R` | 1000 runs, 500 + 2000 penalties | 26 min, 6 cores |
| Figure 6 | `fig6_eeg_heat.R` | same | 37 min, 6 cores |
| Figure 7 | `fig7_power_gcv_decreasing.R` | same | 10 min, 6 cores |
| Figure 8 | `fig8_power_gcv_increasing.R` | same | 6 min, 6 cores |
| Table 1 | `table1_type1_error.R` | 500 runs x 12 settings | runs 1 to 10 only (about 30 min per 5 runs) |
| Table 2 | `table2_power_fixed.R` | same | runs 1 to 10 only |
| Table 3 | `table3_power_random.R` | same | runs 1 to 10 only |
| Table 4 | `table4_eeg_pvalues.R` | 5 triples | 20 s per version, 5 cores |

The machine used for these times was shared with other jobs. The EEG power
curves take about 3 to 4 minutes per triple and are cached in `results/`, so
Figures 5 to 8 share the work (times above count only new triples). At the
rates seen here a full run of one of Tables 1 to 3 (6000 tests) takes more
than a day, so only replications 1 to 10 were run. The simulation scripts
write one file per replication and skip the ones already written, so a full
run can be split, for example
`Rscript scripts/table1_type1_error.R --runs=1:100 --cores=8`.

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
5. The side-by-side tables of this README equal the results
   (`checks/check_readme.py`; set `ADAMANT_RUNS` to the archived `Run0204b`
   folder for the third table).

## Best fit and correct version, side by side

The repository holds two versions of the method.

- **Best fit (the default).** The code that made the paper, with the settings
  and the one replayed rounding (note 1) that come closest to the published
  numbers.
- **Correct (`--method=paper`, Table 4 `--version=paper`).** The method as the
  paper writes it (Section 2, the AdaMantCV algorithm). It changes four
  things, each a bug by the paper's own definitions:
  1. Ties count as ties. P(b) = #{c : T(c) >= T(b)} / (B + 1), with values
     within a relative 1.5e-8 counted as equal, as the paper's indicator
     I(T(b) <= T(c)) says. The original code ranks tied statistics by rounding
     (note 1).
  2. Counts are exact. P_CV = #{b : P(b) <= P(0)} / (B + 1), counted, not
     summed step by step, so a p-value on the boundary is not pushed above
     0.05 by rounding (note 3).
  3. Every permutation is scored by GCV and LOOCV on its own response. The
     paper chooses lambda(b) by the CV score of y(b); the original code fits
     y(b) but measures the residual against the observed y (note 9).
  4. At lambda = 0 the smoother is X W, the true hat matrix, with the limit
     scores of Hastie et al. (2019) where 1 - tr(L) / n is zero. The original
     code uses I minus the inverse, which is not the hat matrix the paper
     defines (note 9).

  The gamma p-value of a point mass (all permuted statistics equal) is 1 in
  this version. The 2021 power code permutes the rows of U rather than the
  response; that is not a bug, since crossprod(U[s, ], y) = U' y[order(s)]
  and order(s) is again a uniform random permutation, so it is the same in
  both versions.

Table 4, rerun three ways (`Rscript scripts/table4_eeg_pvalues.R
--version=best_fit`, `original`, `paper`):

| Triple | Paper | Best fit (default) | Original code, this machine | Correct (`--version=paper`) |
|---|---|---|---|---|
| Fp1 - T8, delta | 0.016 / 0.015 / 0.012 / 0.012 | 0.016 / 0.015 / 0.012 / 0.012 | 0.790 / 0.785 / 0.720 / 0.719 | 1.000 / 1.000 / 1.000 / 1.000 |
| CPz - F8, delta | 0.009 / 0.008 / 0.007 / 0.007 | 0.009 / 0.008 / 0.007 / 0.007 | 0.971 / 0.970 / 0.855 / 0.856 | 1.000 / 1.000 / 1.000 / 1.000 |
| FC5 - O2, delta | 0.004 / 0.003 / 0.005 / 0.004 | 0.004 / 0.003 / 0.005 / 0.004 | 0.179 / 0.184 / 0.051 / 0.049 | 1.000 / 1.000 / 1.000 / 1.000 |
| F3 - FC1, alpha | 0.037 / 0.038 / 0.006 / 0.007 | 0.037 / 0.038 / 0.006 / 0.007 | 0.522 / 0.520 / 0.520 / 0.517 | 1.000 / 1.000 / 1.000 / 1.000 |
| P4 - PO8, gamma | 0.034 / 0.034 / 0.014 / 0.015 | 0.034 / 0.034 / 0.014 / 0.015 | 0.214 / 0.214 / 0.123 / 0.121 | 1.000 / 1.000 / 1.000 / 1.000 |

Each cell: GCV / LOOCV / gamma GCV / gamma LOOCV.

Every table and figure:

| Item | Best fit (default) against the paper | Correct option (`--method=paper`) |
|---|---|---|
| Table 1 | archived runs give the paper: 48/48 cells; runs 1 to 10 rerun equal to the archive: heteroskedastic 240/240, compound symmetric 4/240 | runs 1 to 2 of settings 2, 3, 11, 12 (original scenario numbers): p-values 0.134 to 0.750 |
| Table 2 | archived runs give the paper: 48/48 cells; runs 1 to 10 rerun equal to the archive: heteroskedastic 126/240, compound symmetric 142/240 | runs 1 to 2 of settings 20, 21, 29, 30 (original scenario numbers): every p-value 1 |
| Table 3 | archived runs give the paper: 48/48 cells; runs 1 to 10 rerun equal to the archive: heteroskedastic 134/240, compound symmetric 122/240 | not run yet (job N14) |
| Table 4 | 20/20 cells equal to the paper | all 20 p-values 1 |
| Figure 1A | 10 minima within 0.006 pt of the published points | same code |
| Figure 1B | closest of 17 session states tried: largest offset 14.9 pt, log penalty within 0.92, CV MSE within 0.032 (folds lost, see note 10) | same code |
| Figure 2 | drawn from the archived runs: 4/4 checks OK | power at lambda != 0 equal to the default (largest difference 0); at lambda = 0 power 0 for all 10 variables |
| Figure 3 | drawn from the archived runs: 4/4 checks OK; rerun mean chosen penalty, GCV and LOOCV (paper text in brackets): C -9.706 (-9.867), -10.023 (-10.023), D -10.620 (-10.623), -10.687 (-10.689) | not run yet (job N14) |
| Figures 5, 6 | drawn from the archived runs: Figure 5 4/4 checks OK, Figure 6 4/4 checks OK | not run |
| Figures 7, 8 | drawn from the archived runs: Figure 7 6/6 checks OK, Figure 8 6/6 checks OK | not run |

How often the observed data chose lambda = 0 in the archived runs of Tables 1
to 3 (all 500 runs of each setting, from the chosen-penalty files of
`Run0204b`; `checks/search/sim_lambda0_count.py`):

| Table, p | Observed data chose lambda = 0 | Chose lambda = +-sqrt(eps) | Runs |
|---|---|---|---|
| Table 1, p = 500 | 313 | 0 | 3000 |
| Table 1, p = 1000 | 0 | 189 | 3000 |
| Table 2, p = 500 | 2670 | 0 | 3000 |
| Table 2, p = 1000 | 0 | 2753 | 3000 |
| Table 3, p = 500 | 2888 | 0 | 3000 |
| Table 3, p = 1000 | 0 | 2936 | 3000 |

What these tables say.

- **Table 4** is matched in all 20 cells only by replaying the cluster's
  rounding (note 1). The same code on this machine gives 0.790 for Fp1 - T8.
  By the paper's method every p-value is 1: the coherence is standardized and
  X is centred, so with n < p the ridge fit at lambda = 0 interpolates every
  permuted response exactly, GCV is 0 for all of them, all 1000 permutations
  choose lambda = 0 (check in `checks/RESULTS.md`) and every statistic ties.
- **Tables 2 and 3, p = 500**: in most runs the observed data chose lambda = 0
  (table above), where the p-value is set by rounding. These rows of the paper
  measure rounding, not power, and cannot be repeated on other hardware.
  **p = 1000**: the observed data chose lambda = sqrt(eps), the statistic is
  not tied, and the rows (power 0.996 to 1.000 in the paper) are robust. The
  compound symmetric rows also depend on the LAPACK build (note 2). No other
  version of the scripts or settings came closer: the archived runs of
  `Run0204b` already give all 144 cells of Tables 1 to 3, so what differs is
  the machine, not the code.
- **Correct version of Tables 1 to 3**: only runs 1 and 2 of four settings of
  Tables 1 and 2 were run tonight (job N14 runs them in full). In those runs
  of Table 2 every p-value is 1, for the same reason as Table 4: the response
  is X beta plus small noise, and the fit at lambda near 0 is the same for
  every permutation up to the mean, which a permutation does not change.
- **Figure 1B**: the original set no seed before its folds (note 10). The
  default now uses the closest of 17 session states tried, not a recovered
  seed.
- **Figures 2 and 5 to 8, correct version**: only the lambda = 0 point
  changes (power 0, since every permuted statistic ties). Figure 2 was rerun
  that way; Figures 5 to 8 were not rerun tonight.

## What was reproduced, and what was not

- **Figures 2, 3, 5 to 8, drawn from the archived data**: every tile colour,
  bar and point of the published files is reproduced.
- **Figure 1A**: the ten minimum points match the published figure.
  **Figure 1B**: see note 10.
- **Figure 2, recomputed in full**: equal to the archived run to about 1e-13,
  except at lambda = 0 (note 1).
- **Figures 5 to 8, recomputed in full** (26 triples, 1000 runs each): the
  power curves equal the archived runs to about 1e-14, except at lambda = 0
  (note 1); the GCV and LOOCV curves equal the archive (largest relative
  difference about 1e-7, at penalties where the scores reach about 1e10).
  Because the lambda = 0 value can set the end of the colour range, the
  recomputed heat maps can differ in colour from the published ones even where
  the numbers agree.
- **Table 4, recomputed in full**: see the side-by-side table above.
- **Tables 1 to 3 and Figure 3, recomputed**: only the heteroskedastic
  settings of Table 1 repeat exactly here (note 2), and only replications 1 to
  10 of 500 were run. For Figure 3 (compound symmetric design) the mean chosen
  penalties are close to the values in the paper text but not equal.

Notes, each found by the checks:

1. **Rounding decides the statistic at lambda = 0.** With n < p, U from the
   SVD of X is square, so at lambda = 0 the Mantel statistic equals the
   squared norm of Y for every permutation, and its rank among the
   permutations is set by rounding. In the EEG data GCV and LOOCV
   choose lambda = 0 for the observed data in 250 of the 250 triples of
   `Real0210a`. The five triples of Table 4 are five of the 7 smallest
   p-values of those 250. The original code run here gives 0.790, 0.785, 0.720, 0.719 for
   Fp1 - T8 with the same chosen penalties as the archive. The per-permutation
   p-values the cluster saved (`PctRk_<k>q_0.5.csv` in `Real0210a`) differ
   from this machine only in the observed data at lambda = 0. The best-fit
   default therefore replays those two numbers per triple (how many permuted
   statistics fell below the observed one, and the gamma p-value), stored in
   `reference/table4_cluster_lambda0.csv`, and computes everything else; it
   gives all 20 published cells. They cannot be recomputed: an attempt to run
   the 2019 R 3.6.0 build under emulation was not allowed on this machine.
2. **Simulated designs depend on the LAPACK build.** `MASS::mvrnorm` uses an
   eigen decomposition. The compound symmetric covariance has p - 1 equal
   eigenvalues, so its eigenvectors, and the simulated X, differ between
   builds. For the diagonal (heteroskedastic) covariance only the signs of the
   eigenvectors can differ (they are all positive here). The heteroskedastic
   replications of Table 1 match the archive exactly. Those of Tables 2 and 3
   at p = 500 differ because their observed p-values are set by rounding at
   lambda = 0 (table above).
3. **Rejections at p = 0.05.** The original runs used dplyr 0.8.5, whose
   compiled `percent_rank()` adds the steps one at a time, so a p-value of
   0.050 came out slightly above 0.05 and was not counted. The default
   reproduces that arithmetic; the correct version counts exactly.
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
9. Kept from the original code in the default, with a comment at each place,
   and changed in the correct version: the GCV residual of every permutation
   uses the observed response (a bug: the paper scores each permuted response
   by its own fit); at lambda = 0 with n <= p the smoother is taken as I minus
   the inverse (a bug: that is not the hat matrix, and its near-zero residual
   for the observed data is why the observed data choose lambda = 0). The
   2021 power code permutes the rows of U rather than the response; that is
   kept in both versions, since it gives the same permutation distribution.
10. **Figure 1B folds.** The original code set no seed before the
    cross-validation folds of panel B, so the published panel cannot be
    reproduced fold for fold. `checks/search/fig1b_seed_search.R` tried 17
    seeds and session states; `checks/search/fig1b_score.py` reads the
    published minima on the published axes and ranks them. The closest, now
    the default, is the state the original script leaves after its panel A
    loop over k = 1 and 2.

## Origin of each script

Each script header names the original files it replaces. The original run
folders are on the external drive (`DUMP 010824/Active/Runs/`).
