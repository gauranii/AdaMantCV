# Check results

Written by `checks/run_checks.sh`: the outputs of the scripts in `scripts/`
(run with their defaults) compared with the paper by code.

## Tables 1 to 4

Every cell is compared with the number printed in the paper
(`checks/paper_values.csv`). Each cell below reads paper / computed, in the
order GCV, LOOCV, gamma GCV, gamma LOOCV.

### Table 1: 48 of 48 cells equal

| Row | GCV | LOOCV | gamma GCV | gamma LOOCV | Check |
|---|---|---|---|---|---|
| CS, p = 500, tau = none | 0.042 / 0.042 | 0.042 / 0.042 | 0.036 / 0.036 | 0.036 / 0.036 | OK |
| CS, p = 500, tau = Q1 | 0.030 / 0.030 | 0.028 / 0.028 | 0.032 / 0.032 | 0.028 / 0.028 | OK |
| CS, p = 500, tau = Q2 | 0.048 / 0.048 | 0.050 / 0.050 | 0.052 / 0.052 | 0.052 / 0.052 | OK |
| CS, p = 1000, tau = none | 0.048 / 0.048 | 0.048 / 0.048 | 0.038 / 0.038 | 0.038 / 0.038 | OK |
| CS, p = 1000, tau = Q1 | 0.044 / 0.044 | 0.044 / 0.044 | 0.042 / 0.042 | 0.042 / 0.042 | OK |
| CS, p = 1000, tau = Q2 | 0.036 / 0.036 | 0.036 / 0.036 | 0.036 / 0.036 | 0.036 / 0.036 | OK |
| Heteroskedastic, p = 500, tau = none | 0.038 / 0.038 | 0.038 / 0.038 | 0.038 / 0.038 | 0.038 / 0.038 | OK |
| Heteroskedastic, p = 500, tau = Q1 | 0.026 / 0.026 | 0.026 / 0.026 | 0.028 / 0.028 | 0.028 / 0.028 | OK |
| Heteroskedastic, p = 500, tau = Q2 | 0.046 / 0.046 | 0.046 / 0.046 | 0.050 / 0.050 | 0.050 / 0.050 | OK |
| Heteroskedastic, p = 1000, tau = none | 0.050 / 0.050 | 0.050 / 0.050 | 0.040 / 0.040 | 0.040 / 0.040 | OK |
| Heteroskedastic, p = 1000, tau = Q1 | 0.046 / 0.046 | 0.046 / 0.046 | 0.046 / 0.046 | 0.046 / 0.046 | OK |
| Heteroskedastic, p = 1000, tau = Q2 | 0.040 / 0.040 | 0.040 / 0.040 | 0.038 / 0.038 | 0.038 / 0.038 | OK |

### Table 2: 48 of 48 cells equal

| Row | GCV | LOOCV | gamma GCV | gamma LOOCV | Check |
|---|---|---|---|---|---|
| CS, p = 500, tau = none | 0.062 / 0.062 | 0.062 / 0.062 | 0.104 / 0.104 | 0.104 / 0.104 | OK |
| CS, p = 500, tau = Q1 | 0.342 / 0.342 | 0.290 / 0.290 | 0.378 / 0.378 | 0.326 / 0.326 | OK |
| CS, p = 500, tau = Q2 | 0.348 / 0.348 | 0.288 / 0.288 | 0.382 / 0.382 | 0.324 / 0.324 | OK |
| CS, p = 1000, tau = none | 0.998 / 0.998 | 0.998 / 0.998 | 0.998 / 0.998 | 0.998 / 0.998 | OK |
| CS, p = 1000, tau = Q1 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | OK |
| CS, p = 1000, tau = Q2 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | OK |
| Heteroskedastic, p = 500, tau = none | 0.064 / 0.064 | 0.064 / 0.064 | 0.178 / 0.178 | 0.178 / 0.178 | OK |
| Heteroskedastic, p = 500, tau = Q1 | 0.074 / 0.074 | 0.072 / 0.072 | 0.186 / 0.186 | 0.184 / 0.184 | OK |
| Heteroskedastic, p = 500, tau = Q2 | 0.118 / 0.118 | 0.068 / 0.068 | 0.218 / 0.218 | 0.180 / 0.180 | OK |
| Heteroskedastic, p = 1000, tau = none | 0.996 / 0.996 | 0.996 / 0.996 | 0.998 / 0.998 | 0.998 / 0.998 | OK |
| Heteroskedastic, p = 1000, tau = Q1 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | OK |
| Heteroskedastic, p = 1000, tau = Q2 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | OK |

### Table 3: 48 of 48 cells equal

| Row | GCV | LOOCV | gamma GCV | gamma LOOCV | Check |
|---|---|---|---|---|---|
| CS, p = 500, tau = none | 0.070 / 0.070 | 0.070 / 0.070 | 0.112 / 0.112 | 0.112 / 0.112 | OK |
| CS, p = 500, tau = Q1 | 0.196 / 0.196 | 0.152 / 0.152 | 0.230 / 0.230 | 0.192 / 0.192 | OK |
| CS, p = 500, tau = Q2 | 0.150 / 0.150 | 0.136 / 0.136 | 0.188 / 0.188 | 0.176 / 0.176 | OK |
| CS, p = 1000, tau = none | 0.998 / 0.998 | 0.998 / 0.998 | 0.998 / 0.998 | 0.998 / 0.998 | OK |
| CS, p = 1000, tau = Q1 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | OK |
| CS, p = 1000, tau = Q2 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | OK |
| Heteroskedastic, p = 500, tau = none | 0.092 / 0.092 | 0.092 / 0.092 | 0.268 / 0.268 | 0.268 / 0.268 | OK |
| Heteroskedastic, p = 500, tau = Q1 | 0.084 / 0.084 | 0.092 / 0.092 | 0.266 / 0.266 | 0.268 / 0.268 | OK |
| Heteroskedastic, p = 500, tau = Q2 | 0.088 / 0.088 | 0.092 / 0.092 | 0.268 / 0.268 | 0.268 / 0.268 | OK |
| Heteroskedastic, p = 1000, tau = none | 0.996 / 0.996 | 0.996 / 0.996 | 0.996 / 0.996 | 0.996 / 0.996 | OK |
| Heteroskedastic, p = 1000, tau = Q1 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | OK |
| Heteroskedastic, p = 1000, tau = Q2 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | 1.000 / 1.000 | OK |

### Table 4: 20 of 20 cells equal

| Row | GCV | LOOCV | gamma GCV | gamma LOOCV | Check |
|---|---|---|---|---|---|
| Fp1 - T8, delta | 0.016 / 0.016 | 0.015 / 0.015 | 0.012 / 0.012 | 0.012 / 0.012 | OK |
| CPz - F8, delta | 0.009 / 0.009 | 0.008 / 0.008 | 0.007 / 0.007 | 0.007 / 0.007 | OK |
| FC5 - O2, delta | 0.004 / 0.004 | 0.003 / 0.003 | 0.005 / 0.005 | 0.004 / 0.004 | OK |
| F3 - FC1, alpha | 0.037 / 0.037 | 0.038 / 0.038 | 0.006 / 0.006 | 0.007 / 0.007 | OK |
| P4 - PO8, gamma | 0.034 / 0.034 | 0.034 / 0.034 | 0.014 / 0.014 | 0.015 / 0.015 | OK |

## Figures 1 to 3 and 5 to 8

Every figure is compared with the published file (`archive/figures/`) by the
marks in the two PDF files. 30 of 30 checks OK.

| Figure | Check | Value | Result |
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
| Fig 1A | minimum points vs computed minima, largest offset (pt) | 10, 0.006 | OK |
| Fig 1B | minimum points vs computed minima, largest offset (pt) | 10, 0.005 | OK |
