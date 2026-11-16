# Check results (2026-11-16)

Written by `checks/run_checks.sh`. "differs" is not always a failure: the
README explains each expected difference.

## 1. Archived run outputs against the paper

```
archive vs paper: 164/164 cells equal (3 decimals)
Tables 1-3 from the 500 per-run p-values: 144/144 with p < 0.05, 106/144 with p <= 0.05
```

## 2. Rebuilt functions against the original functions (same inputs and seeds)

25 of 25 checks OK.

| what | n | max_abs_diff | tolerance | ok |
|---|---|---|---|---|
| adamant_cv triple 11, q = 0.25, 100 perms: lambda grid | 250 | 0 | 0 | OK |
| adamant_cv triple 11, q = 0.25, 100 perms: GCV scores (relative; absolute below 1e-12) | 25250 | 1.20581111762164e-09 | 1e-08 | OK |
| adamant_cv triple 11, q = 0.25, 100 perms: LOOCV scores (relative; absolute below 1e-12) | 25250 | 1.20579235485252e-09 | 1e-08 | OK |
| adamant_cv triple 11, q = 0.25, 100 perms: chosen lambdas GCV, LOO | 202 | 0 | 0 | OK |
| adamant_cv triple 11, q = 0.25, 100 perms: p-values (4) | 4 | 5.55111512312578e-16 | 1e-12 | OK |
| adamant_cv triple 11, q = 0.50, 100 perms: lambda grid | 250 | 0 | 0 | OK |
| adamant_cv triple 11, q = 0.50, 100 perms: GCV scores (relative; absolute below 1e-12) | 25250 | 1.20581111762164e-09 | 1e-08 | OK |
| adamant_cv triple 11, q = 0.50, 100 perms: LOOCV scores (relative; absolute below 1e-12) | 25250 | 1.20579235485252e-09 | 1e-08 | OK |
| adamant_cv triple 11, q = 0.50, 100 perms: chosen lambdas GCV, LOO | 202 | 0 | 0 | OK |
| adamant_cv triple 11, q = 0.50, 100 perms: p-values (4) | 4 | 4.44089209850063e-16 | 1e-12 | OK |
| simulation CS p = 500 Q1, beta zero, runs 1-2, 30 perms: p-values | 8 | 2.77555756156289e-17 | 1e-12 | OK |
| simulation CS p = 500 Q1, beta fixed, runs 1-2, 30 perms: p-values | 8 | 0 | 1e-12 | OK |
| simulation CS p = 500 Q1, beta random, runs 1-2, 30 perms: p-values | 8 | 0 | 1e-12 | OK |
| power curve eeg191, 3 runs, 30 perms: lambdas | 2500 | 0 | 0 | OK |
| power curve eeg191: power | 2500 | 0 | 0 | OK |
| power curve eeg191: power MOM | 2500 | 0 | 0 | OK |
| power curve eeg191: permutation p-values (as mean -log10) | 2500 | 8.88178419700125e-16 | 1e-12 | OK |
| GCV curve eeg191 (relative) | 2500 | 0 | 1e-10 | OK |
| LOOCV curve eeg191 (relative) | 2500 | 0 | 1e-10 | OK |
| power curve liver1, 3 runs, 30 perms: lambdas | 300 | 0 | 0 | OK |
| power curve liver1: power | 300 | 0 | 0 | OK |
| power curve liver1: power MOM | 300 | 0 | 0 | OK |
| power curve liver1: permutation p-values (as mean -log10) | 300 | 8.88178419700125e-16 | 1e-12 | OK |
| GCV curve liver1 (relative) | 300 | 0 | 1e-10 | OK |
| LOOCV curve liver1 (relative) | 300 | 0 | 1e-10 | OK |

## 3. Published figures against the same figures redrawn from the archived data

47 of 48 checks OK.

| figure | check | value | ok |
|---|---|---|---|
| Fig 2A | tiles: count, largest colour difference (0-1) | 3000/3000, 0.0 | OK |
| Fig 2A | row labels equal (ignoring case) | True | OK |
| Fig 2B | tiles: count, largest colour difference (0-1) | 3000/3000, 0.0 | OK |
| Fig 2B | row labels equal (ignoring case) | True | OK |
| Fig 5A | tiles: count, largest colour difference (0-1) | 30000/30000, 0.0 | OK |
| Fig 5A | row labels equal (ignoring case) | True | OK |
| Fig 5B | tiles: count, largest colour difference (0-1) | 30000/30000, 0.0 | OK |
| Fig 5B | row labels equal (ignoring case) | True | OK |
| Fig 6A | tiles: count, largest colour difference (0-1) | 30000/30000, 0.0 | OK |
| Fig 6A | row labels equal (ignoring case) | True | OK |
| Fig 6B | tiles: count, largest colour difference (0-1) | 30000/30000, 0.0 | OK |
| Fig 6B | row labels equal (ignoring case) | True | OK |
| Fig 3A | bars: count, spread of height ratios, same empty bins | 56/56, 0.0009, True | OK |
| Fig 3B | bars: count, spread of height ratios, same empty bins | 62/62, 0.0006, True | OK |
| Fig 3C | bars: count, spread of height ratios, same empty bins | 56/56, 0.0010, True | OK |
| Fig 3D | bars: count, spread of height ratios, same empty bins | 70/70, 0.0012, True | OK |
| Fig 7 | Power223: points, largest offset after scaling (pt) | 2500, 0.009 | OK |
| Fig 7 | Power22: points, largest offset after scaling (pt) | 2500, 0.009 | OK |
| Fig 7 | Power25: points, largest offset after scaling (pt) | 2500, 0.009 | OK |
| Fig 7 | GCV223: points, largest offset after scaling (pt) | 275, 0.008 | OK |
| Fig 7 | GCV22: points, largest offset after scaling (pt) | 166, 0.009 | OK |
| Fig 7 | GCV25: points, largest offset after scaling (pt) | 274, 0.009 | OK |
| Fig 8 | Power191s0.5: points, largest offset after scaling (pt) | 2500, 0.008 | OK |
| Fig 8 | Power184s0.5: points, largest offset after scaling (pt) | 2500, 0.009 | OK |
| Fig 8 | Power160s0.5: points, largest offset after scaling (pt) | 2500, 0.031 | OK |
| Fig 8 | GCV191: points, largest offset after scaling (pt) | 162, 0.009 | OK |
| Fig 8 | GCV184: points, largest offset after scaling (pt) | 158, 0.008 | OK |
| Fig 8 | GCV160: points, largest offset after scaling (pt) | 175, 0.008 | OK |
| Fig 1A | blue minimum points vs rebuilt minima, largest offset (pt) | 10, 0.006 | OK |
| Fig 1B | blue minimum points vs rebuilt minima, largest offset (pt) | 10, 19.929 | differs |
| fig7_gcv_22 | GCV panel drawn from rebuilt table: points | 166 | OK |
| fig7_gcv_223 | GCV panel drawn from rebuilt table: points | 275 | OK |
| fig7_gcv_25 | GCV panel drawn from rebuilt table: points | 274 | OK |
| fig7_power_22 | points vs rebuilt table, largest offset (pt) | 2500, 0.004 | OK |
| fig7_power_223 | points vs rebuilt table, largest offset (pt) | 2500, 0.004 | OK |
| fig7_power_25 | points vs rebuilt table, largest offset (pt) | 2500, 0.004 | OK |
| fig8_gcv_160 | GCV panel drawn from rebuilt table: points | 175 | OK |
| fig8_gcv_184 | GCV panel drawn from rebuilt table: points | 158 | OK |
| fig8_gcv_191 | GCV panel drawn from rebuilt table: points | 162 | OK |
| fig8_power_160 | points vs rebuilt table, largest offset (pt) | 2500, 0.006 | OK |
| fig8_power_184 | points vs rebuilt table, largest offset (pt) | 2500, 0.004 | OK |
| fig8_power_191 | points vs rebuilt table, largest offset (pt) | 2500, 0.004 | OK |
| Fig 2A | rebuilt run vs published: tiles that differ (report only) | 2703 of 3000 | OK |
| Fig 2B | rebuilt run vs published: tiles that differ (report only) | 10 of 3000 | OK |
| Fig 5A | rebuilt run vs published: tiles that differ (report only) | 29517 of 30000 | OK |
| Fig 5B | rebuilt run vs published: tiles that differ (report only) | 11 of 30000 | OK |
| Fig 6A | rebuilt run vs published: tiles that differ (report only) | 23926 of 30000 | OK |
| Fig 6B | rebuilt run vs published: tiles that differ (report only) | 12 of 30000 | OK |

## 4. Rebuilt results against the archive and the paper

148 of 149 checks OK.

| item | check | value | ok |
|---|---|---|---|
| Fig 2 v1 | lambda grid, largest difference | 300, 1.0e-12 | OK |
| Fig 2 v1 | GCV and LOOCV curves, largest relative difference | 5.6e-12 | OK |
| Fig 2 v1 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 2.4e-14 | OK |
| Fig 2 v1 | power at lambda = 0 (tied statistic), difference | 0.003 | OK |
| Fig 2 v2 | lambda grid, largest difference | 300, 1.0e-12 | OK |
| Fig 2 v2 | GCV and LOOCV curves, largest relative difference | 5.5e-12 | OK |
| Fig 2 v2 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 1.0e-13 | OK |
| Fig 2 v2 | power at lambda = 0 (tied statistic), difference | 0.001 | OK |
| Fig 2 v3 | lambda grid, largest difference | 300, 1.0e-12 | OK |
| Fig 2 v3 | GCV and LOOCV curves, largest relative difference | 4.8e-12 | OK |
| Fig 2 v3 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 9.0e-14 | OK |
| Fig 2 v3 | power at lambda = 0 (tied statistic), difference | 0.157 | OK |
| Fig 2 v4 | lambda grid, largest difference | 300, 1.0e-12 | OK |
| Fig 2 v4 | GCV and LOOCV curves, largest relative difference | 5.6e-12 | OK |
| Fig 2 v4 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 7.3e-14 | OK |
| Fig 2 v4 | power at lambda = 0 (tied statistic), difference | 0.498 | OK |
| Fig 2 v5 | lambda grid, largest difference | 300, 1.0e-12 | OK |
| Fig 2 v5 | GCV and LOOCV curves, largest relative difference | 5.5e-12 | OK |
| Fig 2 v5 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 2.2e-14 | OK |
| Fig 2 v5 | power at lambda = 0 (tied statistic), difference | 0.109 | OK |
| Fig 2 v6 | lambda grid, largest difference | 300, 1.0e-12 | OK |
| Fig 2 v6 | GCV and LOOCV curves, largest relative difference | 5.5e-12 | OK |
| Fig 2 v6 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 7.4e-14 | OK |
| Fig 2 v6 | power at lambda = 0 (tied statistic), difference | 0.001 | OK |
| Fig 2 v7 | lambda grid, largest difference | 300, 1.0e-12 | OK |
| Fig 2 v7 | GCV and LOOCV curves, largest relative difference | 5.5e-12 | OK |
| Fig 2 v7 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 2.1e-14 | OK |
| Fig 2 v7 | power at lambda = 0 (tied statistic), difference | 0.100 | OK |
| Fig 2 v8 | lambda grid, largest difference | 300, 1.0e-12 | OK |
| Fig 2 v8 | GCV and LOOCV curves, largest relative difference | 5.5e-12 | OK |
| Fig 2 v8 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 1.7e-14 | OK |
| Fig 2 v8 | power at lambda = 0 (tied statistic), difference | 0.001 | OK |
| Fig 2 v9 | lambda grid, largest difference | 300, 1.0e-12 | OK |
| Fig 2 v9 | GCV and LOOCV curves, largest relative difference | 5.5e-12 | OK |
| Fig 2 v9 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 2.5e-14 | OK |
| Fig 2 v9 | power at lambda = 0 (tied statistic), difference | 0.042 | OK |
| Fig 2 v10 | lambda grid, largest difference | 300, 1.0e-12 | OK |
| Fig 2 v10 | GCV and LOOCV curves, largest relative difference | 5.5e-12 | OK |
| Fig 2 v10 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 2.0e-14 | OK |
| Fig 2 v10 | power at lambda = 0 (tied statistic), difference | 0.197 | OK |
| EEG 123 | lambda grid, largest difference | 2500, 1.0e-12 | OK |
| EEG 123 | GCV and LOOCV curves, largest relative difference | 9.1e-08 | OK |
| EEG 123 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 1.0e-14 | OK |
| EEG 123 | power at lambda = 0 (tied statistic), difference | 0.037 | OK |
| EEG 136 | lambda grid, largest difference | 2500, 1.0e-12 | OK |
| EEG 136 | GCV and LOOCV curves, largest relative difference | 9.1e-08 | OK |
| EEG 136 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 1.5e-14 | OK |
| EEG 136 | power at lambda = 0 (tied statistic), difference | 0.234 | OK |
| EEG 144 | lambda grid, largest difference | 2500, 1.0e-12 | OK |
| EEG 144 | GCV and LOOCV curves, largest relative difference | 9.1e-08 | OK |
| EEG 144 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 6.1e-15 | OK |
| EEG 144 | power at lambda = 0 (tied statistic), difference | 0.012 | OK |
| EEG 148 | lambda grid, largest difference | 2500, 1.0e-12 | OK |
| EEG 148 | GCV and LOOCV curves, largest relative difference | 9.1e-08 | OK |
| EEG 148 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 2.0e-14 | OK |
| EEG 148 | power at lambda = 0 (tied statistic), difference | 0.005 | OK |
| EEG 151 | lambda grid, largest difference | 2500, 1.0e-12 | OK |
| EEG 151 | GCV and LOOCV curves, largest relative difference | 9.1e-08 | OK |
| EEG 151 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 1.0e-14 | OK |
| EEG 151 | power at lambda = 0 (tied statistic), difference | 0.018 | OK |
| EEG 157 | lambda grid, largest difference | 2500, 1.0e-12 | OK |
| EEG 157 | GCV and LOOCV curves, largest relative difference | 9.1e-08 | OK |
| EEG 157 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 1.0e-14 | OK |
| EEG 157 | power at lambda = 0 (tied statistic), difference | 0.007 | OK |
| EEG 159 | lambda grid, largest difference | 2500, 1.0e-12 | OK |
| EEG 159 | GCV and LOOCV curves, largest relative difference | 9.1e-08 | OK |
| EEG 159 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 1.0e-14 | OK |
| EEG 159 | power at lambda = 0 (tied statistic), difference | 0.001 | OK |
| EEG 160 | lambda grid, largest difference | 2500, 1.0e-12 | OK |
| EEG 160 | GCV and LOOCV curves, largest relative difference | 9.1e-08 | OK |
| EEG 160 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 2.0e-14 | OK |
| EEG 160 | power at lambda = 0 (tied statistic), difference | 0.052 | OK |
| EEG 177 | lambda grid, largest difference | 2500, 1.0e-12 | OK |
| EEG 177 | GCV and LOOCV curves, largest relative difference | 9.1e-08 | OK |
| EEG 177 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 1.0e-14 | OK |
| EEG 177 | power at lambda = 0 (tied statistic), difference | 0.000 | OK |
| EEG 179 | lambda grid, largest difference | 2500, 1.0e-12 | OK |
| EEG 179 | GCV and LOOCV curves, largest relative difference | 9.1e-08 | OK |
| EEG 179 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 1.0e-14 | OK |
| EEG 179 | power at lambda = 0 (tied statistic), difference | 0.042 | OK |
| EEG 181 | lambda grid, largest difference | 2500, 1.0e-12 | OK |
| EEG 181 | GCV and LOOCV curves, largest relative difference | 9.1e-08 | OK |
| EEG 181 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 1.5e-14 | OK |
| EEG 181 | power at lambda = 0 (tied statistic), difference | 0.001 | OK |
| EEG 183 | lambda grid, largest difference | 2500, 1.0e-12 | OK |
| EEG 183 | GCV and LOOCV curves, largest relative difference | 9.1e-08 | OK |
| EEG 183 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 1.0e-14 | OK |
| EEG 183 | power at lambda = 0 (tied statistic), difference | 0.023 | OK |
| EEG 184 | lambda grid, largest difference | 2500, 1.0e-12 | OK |
| EEG 184 | GCV and LOOCV curves, largest relative difference | 9.1e-08 | OK |
| EEG 184 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 2.0e-14 | OK |
| EEG 184 | power at lambda = 0 (tied statistic), difference | 0.013 | OK |
| EEG 191 | lambda grid, largest difference | 2500, 1.0e-12 | OK |
| EEG 191 | GCV and LOOCV curves, largest relative difference | 9.1e-08 | OK |
| EEG 191 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 1.0e-14 | OK |
| EEG 191 | power at lambda = 0 (tied statistic), difference | 0.131 | OK |
| EEG 200 | lambda grid, largest difference | 2500, 1.0e-12 | OK |
| EEG 200 | GCV and LOOCV curves, largest relative difference | 9.1e-08 | OK |
| EEG 200 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 4.1e-15 | OK |
| EEG 200 | power at lambda = 0 (tied statistic), difference | 0.002 | OK |
| EEG 205 | lambda grid, largest difference | 2500, 1.0e-12 | OK |
| EEG 205 | GCV and LOOCV curves, largest relative difference | 9.1e-08 | OK |
| EEG 205 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 1.0e-14 | OK |
| EEG 205 | power at lambda = 0 (tied statistic), difference | 0.157 | OK |
| EEG 219 | lambda grid, largest difference | 2500, 1.0e-12 | OK |
| EEG 219 | GCV and LOOCV curves, largest relative difference | 9.1e-08 | OK |
| EEG 219 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 3.0e-14 | OK |
| EEG 219 | power at lambda = 0 (tied statistic), difference | 0.017 | OK |
| EEG 220 | lambda grid, largest difference | 2500, 1.0e-12 | OK |
| EEG 220 | GCV and LOOCV curves, largest relative difference | 9.1e-08 | OK |
| EEG 220 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 1.0e-14 | OK |
| EEG 220 | power at lambda = 0 (tied statistic), difference | 0.105 | OK |
| EEG 223 | lambda grid, largest difference | 2500, 1.0e-12 | OK |
| EEG 223 | GCV and LOOCV curves, largest relative difference | 9.1e-08 | OK |
| EEG 223 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 2.0e-14 | OK |
| EEG 223 | power at lambda = 0 (tied statistic), difference | 0.041 | OK |
| EEG 22 | lambda grid, largest difference | 2500, 1.0e-12 | OK |
| EEG 22 | GCV and LOOCV curves, largest relative difference | 9.1e-08 | OK |
| EEG 22 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 2.0e-14 | OK |
| EEG 22 | power at lambda = 0 (tied statistic), difference | 0.021 | OK |
| EEG 231 | lambda grid, largest difference | 2500, 1.0e-12 | OK |
| EEG 231 | GCV and LOOCV curves, largest relative difference | 9.1e-08 | OK |
| EEG 231 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 1.6e-14 | OK |
| EEG 231 | power at lambda = 0 (tied statistic), difference | 0.009 | OK |
| EEG 248 | lambda grid, largest difference | 2500, 1.0e-12 | OK |
| EEG 248 | GCV and LOOCV curves, largest relative difference | 9.1e-08 | OK |
| EEG 248 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 3.0e-14 | OK |
| EEG 248 | power at lambda = 0 (tied statistic), difference | 0.105 | OK |
| EEG 25 | lambda grid, largest difference | 2500, 1.0e-12 | OK |
| EEG 25 | GCV and LOOCV curves, largest relative difference | 9.1e-08 | OK |
| EEG 25 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 2.0e-14 | OK |
| EEG 25 | power at lambda = 0 (tied statistic), difference | 0.165 | OK |
| EEG 31 | lambda grid, largest difference | 2500, 1.0e-12 | OK |
| EEG 31 | GCV and LOOCV curves, largest relative difference | 9.1e-08 | OK |
| EEG 31 | power table (1000 runs) vs archive (1000 runs), lambda != 0 | 1.0e-14 | OK |
| EEG 31 | power at lambda = 0 (tied statistic), difference | 0.095 | OK |
| Table 4 | archived run outputs equal to paper (cells) | 20/20 | OK |
| Table 4 | rebuilt equal to paper (cells) | 0/20 | differs |
| Table 4 | rebuilt equal to ORIGINAL code rerun on this machine | 20/20 | OK |
| Table 4 | original rerun: lambdas chosen per permutation = archive | True | OK |
| Table 4 | archive: lambda chosen by GCV for the observed data | 0 0 0 0 0 | OK |
| Fig 3C | mean chosen lambda GCV, LOOCV (paper text in brackets) | -9.706 (-9.867), -10.023 (-10.023) | OK |
| Fig 3D | mean chosen lambda GCV, LOOCV (paper text in brackets) | -10.620 (-10.623), -10.687 (-10.689) | OK |
| Table1 | CS rows: rebuilt per-run p-values = archive (runs 1-10) | 4/240 | OK |
| Table1 | Het rows: rebuilt per-run p-values = archive (runs 1-10) | 240/240 | OK |
| Table2 | CS rows: rebuilt per-run p-values = archive (runs 1-10) | 142/240 | OK |
| Table2 | Het rows: rebuilt per-run p-values = archive (runs 1-10) | 126/240 | OK |
| Table3 | CS rows: rebuilt per-run p-values = archive (runs 1-10) | 122/240 | OK |
| Table3 | Het rows: rebuilt per-run p-values = archive (runs 1-10) | 134/240 | OK |

