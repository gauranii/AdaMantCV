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
| `archive/` | Outputs of the original runs that the scripts read, the published figure files (`figures/`, named as the figures the scripts draw), and `reduce_archive.R`, which made the former from the run folders |
| `data/liver/` | The liver toxicity data |
| `coherence/` | The scripts that computed the EEG coherence and built the per-triple files, with how each step connects to the data folders (not run by the flake) |
| `reference/` | The 250 EEG triples (`triple.csv`), the channel names, and `table4_cluster_lambda0.csv` |
| `checks/` | The numbers printed in the paper (`paper_values.csv`), the comparison scripts and their result (`RESULTS.md`), and `check_equivalence.R`, which compares the outputs of two versions of the code; `check_python.sh`, which compares the outputs of the R and the Python scripts, and its result (`PYTHON.md`); `versions.sh`, which prints the software versions |
| `flake.nix`, `flake.lock` | Pinned environment (R with ggplot2, MASS and glmnet; Python with numpy, scipy, pandas, matplotlib and numba), one command per output, and the checks |

## Running

With [Nix](https://nixos.org) (flakes enabled), from the repository root:

```sh
nix run .#table4        # one output: table1 to table4, fig1 to fig3, fig5 to fig8
nix run .#all           # every output
nix run .#py-table4     # the same output with the Python scripts (py-table1, ..., py-fig8, py-all)
nix flake check         # build every output, compare it with the paper, compare R with Python
nix run .#versions      # print the software versions
```

Tables are written to `results/` (csv) and figures to `figures/` (pdf), each
under the name of its table or figure in the paper (see Outputs).
`nix flake check` compares every cell of Tables 1 to 4 with the paper and
every figure with the published file; `checks/RESULTS.md` holds the comparison.
It also builds every output with the Python scripts and compares it with the
output of the R scripts, cell by cell for the tables and mark by mark for the
numbers the figures plot; `checks/PYTHON.md` holds that comparison.

To recompute an output from the data instead of the archived runs, add
`--rerun`. The EEG outputs then need the EEG data (see Data):

```sh
export ADAMANT_DATA=/path/to/Data      # the folder that holds 04_coherence_SNPs/
nix run .#table4 -- --rerun --cores=5
```

`nix develop` opens a shell with the same R and Python, where the scripts can
be run directly: `Rscript scripts/table4_eeg_pvalues.R [--rerun]` or
`python python/scripts/table4_eeg_pvalues.py [--rerun]`. Each script lists
its options in its header. The settings are those of the paper: seed 8891, 1000
permutations, significance level 0.05.

## Outputs

Each table and figure of the paper, the scripts that build it, what they read,
the files they write (figures under `figures/`) and how `nix flake check`
checks them. The R and the Python script of an output write the same file
names.

| Paper | R script | Python script | Input | Output | Check |
|---|---|---|---|---|---|
| Table 1 | `scripts/table1_type1_error.R` | `python/scripts/table1_type1_error.py` | `archive/run0204b/table1_pvalues.csv` | `results/Table1_type1-error.csv` | every cell against the paper (`checks/paper_values.csv`); R against Python, cell by cell |
| Table 2 | `scripts/table2_power_fixed.R` | `python/scripts/table2_power_fixed.py` | `archive/run0204b/table2_pvalues.csv` | `results/Table2_power_fixed-effects.csv` | every cell against the paper (`checks/paper_values.csv`); R against Python, cell by cell |
| Table 3 | `scripts/table3_power_random.R` | `python/scripts/table3_power_random.py` | `archive/run0204b/table3_pvalues.csv` | `results/Table3_power_random-effects.csv` | every cell against the paper (`checks/paper_values.csv`); R against Python, cell by cell |
| Table 4 | `scripts/table4_eeg_pvalues.R` | `python/scripts/table4_eeg_pvalues.py` | `archive/real0210a/table4_pvalues.csv`, `reference/triple.csv`, `reference/chan_names.csv` | `results/Table4_EEG_pvalues.csv` | every cell against the paper (`checks/paper_values.csv`); R against Python, cell by cell |
| Figure 1 | `scripts/fig1_liver_cv.R` | `python/scripts/fig1_liver_cv.py` | `data/liver/mixOmicsData.csv` | `Fig1A_liver_ridge-CV-error_p50-random-genes.pdf`<br>`Fig1B_liver_ridge-CV-error_p3116-all-genes.pdf` | minimum points of the published files in `archive/figures/` against the computed minima (`results/fig1_minima.csv`); R against Python, the plotted numbers |
| Figure 2 | `scripts/fig2_liver_heat.R` | `python/scripts/fig2_liver_heat.py` | `archive/liver/curves_liver<j>.csv.gz`, j = 1 to 10 | `Fig2A_liver_heatmap_mean-log10-pvalue.pdf`<br>`Fig2B_liver_heatmap_power.pdf` | marks of each PDF against the published file of the same name in `archive/figures/`; R against Python, the plotted numbers |
| Figure 3 | `scripts/fig3_lambda_hist.R` | `python/scripts/fig3_lambda_hist.py` | `archive/run0204b/fig3_lambdas.csv` | `Fig3A_simulation_lambda-histogram_null_threshold-Q1.pdf`<br>`Fig3B_simulation_lambda-histogram_null_threshold-median.pdf`<br>`Fig3C_simulation_lambda-histogram_alternative_threshold-Q1.pdf`<br>`Fig3D_simulation_lambda-histogram_alternative_threshold-median.pdf` | marks of each PDF against the published file of the same name in `archive/figures/`; R against Python, the plotted numbers |
| Figure 5 | `scripts/fig5_eeg_heat.R` | `python/scripts/fig5_eeg_heat.py` | `archive/coh7/curves_eeg<k>.csv.gz` (12 triples), `reference/triple.csv`, `reference/chan_names.csv` | `Fig5A_EEG_heatmap_mean-log10-pvalue_top-triples.pdf`<br>`Fig5B_EEG_heatmap_power_top-triples.pdf` | marks of each PDF against the published file of the same name in `archive/figures/`; R against Python, the plotted numbers |
| Figure 6 | `scripts/fig6_eeg_heat.R` | `python/scripts/fig6_eeg_heat.py` | `archive/coh7/curves_eeg<k>.csv.gz` (12 triples), `reference/triple.csv`, `reference/chan_names.csv` | `Fig6A_EEG_heatmap_mean-log10-pvalue_other-triples.pdf`<br>`Fig6B_EEG_heatmap_power_other-triples.pdf` | marks of each PDF against the published file of the same name in `archive/figures/`; R against Python, the plotted numbers |
| Figure 7 | `scripts/fig7_power_gcv_decreasing.R` | `python/scripts/fig7_power_gcv_decreasing.py` | `archive/coh7/curves_eeg<k>.csv.gz`, k = 223, 22, 25, `reference/triple.csv`, `reference/chan_names.csv` | `Fig7A_EEG_power-decreasing_F3-FT8-beta.pdf`<br>`Fig7B_EEG_power-decreasing_Iz-Oz-delta.pdf`<br>`Fig7C_EEG_power-decreasing_CPz-AF4-delta.pdf`<br>`Fig7D_EEG_GCV_F3-FT8-beta.pdf`<br>`Fig7E_EEG_GCV_Iz-Oz-delta.pdf`<br>`Fig7F_EEG_GCV_CPz-AF4-delta.pdf` | marks of each PDF against the published file of the same name in `archive/figures/`; R against Python, the plotted numbers |
| Figure 8 | `scripts/fig8_power_gcv_increasing.R` | `python/scripts/fig8_power_gcv_increasing.py` | `archive/coh7/curves_eeg<k>.csv.gz`, k = 191, 184, 160, `reference/triple.csv`, `reference/chan_names.csv` | `Fig8A_EEG_power-increasing_FC3-P7-beta.pdf`<br>`Fig8B_EEG_power-increasing_F3-O2-alpha.pdf`<br>`Fig8C_EEG_power-increasing_Pz-PO4-alpha.pdf`<br>`Fig8D_EEG_GCV_FC3-P7-beta.pdf`<br>`Fig8E_EEG_GCV_F3-O2-alpha.pdf`<br>`Fig8F_EEG_GCV_Pz-PO4-alpha.pdf` | marks of each PDF against the published file of the same name in `archive/figures/`; R against Python, the plotted numbers |

With `--rerun`, Table 4 and Figures 5 to 8 read
`$ADAMANT_DATA/04_coherence_SNPs/dat_ch<k>.csv` instead of the archive (Table 4
also `reference/table4_cluster_lambda0.csv`), Figure 2 reads
`data/liver/mixOmicsData.csv`, and Tables 1 to 3 and Figure 3 simulate their
data. The comparisons with the paper are in `checks/RESULTS.md` and those of R
with Python in `checks/PYTHON.md`.

Commands and run times:

| Paper | Command | From archive | With `--rerun` |
|---|---|---|---|
| Table 1 | `nix run .#table1` | 1 s | about 9 h, 8 cores (estimate) |
| Table 2 | `nix run .#table2` | 1 s | about 9 h, 8 cores (estimate) |
| Table 3 | `nix run .#table3` | 1 s | about 9 h, 8 cores (estimate) |
| Table 4 | `nix run .#table4` | 1 s | 35 s, 5 cores |
| Figure 1 | `nix run .#fig1` | 20 s (always computed from the data) | same |
| Figure 2 | `nix run .#fig2` | 1 s | 2 min, 7 cores |
| Figure 3 | `nix run .#fig3` | 1 s | 1 min, 4 cores |
| Figure 5 | `nix run .#fig5` | 1 s | about 55 min, 7 cores (estimate) |
| Figure 6 | `nix run .#fig6` | 1 s | about 55 min, 7 cores (estimate) |
| Figure 7 | `nix run .#fig7` | 1 s | 14 min, 7 cores |
| Figure 8 | `nix run .#fig8` | 1 s | 14 min, 7 cores |

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

- **Liver toxicity data** (`data/liver/mixOmicsData.csv`): 64 rats, 10
  clinical chemistry measurements and 3116 gene expressions; the public
  `liver.toxicity` data of the Bioconductor package
  [mixOmics](https://mixomics.org), exported to csv on 2021-12-05
  (`library(mixOmics); data(liver.toxicity)`). The data are from Bushel,
  Wolfinger and Gibson (2007), BMC Systems Biology 1:15,
  https://doi.org/10.1186/1752-0509-1-15, and Heinloth et al. (2004),
  Toxicological Sciences 80:193-202, https://doi.org/10.1093/toxsci/kfh145.
  Used by Figures 1 and 2.
- **EEG coherence and SNPs** of 350 subjects: one file
  `04_coherence_SNPs/dat_ch<k>.csv` per triple k of channel pair and frequency
  band (a row id, the coherence, then 497 SNPs), 250 files. These
  subject-level data are not in the repository and are needed only for
  `--rerun` of Table 4 and Figures 5 to 8; the scripts look for the folder
  `04_coherence_SNPs/` in the folder named by `ADAMANT_DATA`. That data
  folder is laid out in the order of the pipeline:

  | Folder | Content |
  |---|---|
  | `01_raw/` | Genotypes of the 350 subjects (`eeg_dat_snp_all.RData`) and the 500 candidate SNPs (`cups_top_500_SNPs.txt`) |
  | `02_coherence/` | Band coherence of the 64 channels per subject (`vwm_coh_clean.RData`, 350 x 64 x 64 x 5), computed by the scripts in `coherence/`, and the earlier test results that came with it (`coh_results.RData`) |
  | `03_triples_channels/` | `triple.csv` and `chan_names.csv` (in the repository as `reference/`) |
  | `04_coherence_SNPs/` | `dat_ch<k>.csv`, k = 1 to 250, read by the scripts |
  | `liver_toxicity/` | The liver toxicity data (in the repository as `data/liver/`) |
- **Archived runs** (`archive/`): per-replication p-values of Tables 1 to 3,
  chosen penalties of Figure 3, per-permutation p-values of Table 4, and the
  power and cross-validation curves of Figures 2 and 5 to 8. They hold no
  subject-level data.

## Software versions

`flake.lock` pins the nixpkgs revision, and with it the version of R, Python
and every package (there is no `renv.lock`). `nix run .#versions` (or
`bash checks/versions.sh` inside `nix develop`) prints them; `nix flake check`
compares this section with its output.

```text
nixpkgs  NixOS/nixpkgs f3aa7d3d72c2bc66e0f30f4262d7ec2e33c9eab2

Python   3.14.7
  numpy        2.5.2
  scipy        1.18.1
  pandas       3.0.4
  matplotlib   3.11.1
  numba        0.67.0

R        4.6.1
  boot         1.3-32
  class        7.3-23
  cli          3.6.6
  cluster      2.1.8.2
  codetools    0.2-20
  cpp11        0.5.5
  farver       2.1.2
  foreach      1.5.2
  foreign      0.8-91
  ggplot2      4.0.3
  glmnet       5.0
  glue         1.8.1
  gtable       0.3.6
  isoband      0.3.0
  iterators    1.0.14
  KernSmooth   2.23-26
  labeling     0.4.3
  lattice      0.22-9
  lifecycle    1.0.5
  MASS         7.3-65
  Matrix       1.7-5
  mgcv         1.9-4
  nlme         3.1-169
  nnet         7.3-20
  R6           2.6.1
  RColorBrewer 1.1-3
  Rcpp         1.1.2
  RcppEigen    0.3.4.0.2
  rlang        1.3.0
  rpart        4.1.27
  S7           0.2.2
  scales       1.4.0
  shape        1.4.6.1
  spatial      7.3-18
  survival     3.8-6
  vctrs        0.7.3
  viridisLite  0.4.3
  withr        3.0.3

OpenBLAS 0.3.33
```
