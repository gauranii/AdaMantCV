# -----------------------------------------------------------------------------
# File:    checks/archive_to_cache.py
# Purpose: Write the archived 2021/2022 run outputs in the cache format of the
#          rebuilt scripts, so the figure scripts can redraw the published
#          figures from the very data that made them (no recomputation).
# Usage:   python3 checks/archive_to_cache.py <archive dir> <results dir>
#          then run the figure scripts in a copy of the repo whose results/
#          is <results dir>.
# -----------------------------------------------------------------------------
import csv
import os
import sys

archive, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
POWER = ['lambdas', 'ave_logpval', 'ave_logpval_MOM', 'sd_logpval', 'sd_logpval_MOM',
         'power', 'power_MOM', 'power_check', 'sign']


def write(name, header, rows):
    with open(os.path.join(out, name), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def curve(src, tag, grid):
    rows = list(csv.DictReader(open(src)))
    write(f'power_{tag}_runs1000_perms1000_grid{grid}.csv', POWER,
          [[r[c] for c in POWER] for r in rows])
    write(f'gcv_{tag}_grid{grid}.csv', ['lambdas', 'gcv', 'loo'],
          [[r['lambdas'], r['l_opt$gcv'], r['l_opt$loo']] for r in rows])


for f in os.listdir(os.path.join(archive, 'coh7')):
    if f.startswith('CV_ch'):
        curve(os.path.join(archive, 'coh7', f), 'eeg' + f[5:-4], '500x2000')
for j in range(1, 11):
    curve(os.path.join(archive, 'liver', f'liver_CV_ch{j}.csv'), f'liver{j}', '50x250')

# Fig 3: per-permutation lambdas of replication 500 (full precision), from
# T1_summary (beta = 0) and P_summary (fixed beta) of the CS, p = 1000 runs.
rows = []
for panel, beta, tau, name in [('a', 'zero', 'Q1', 'T1_summary_14'), ('b', 'zero', 'Q2', 'T1_summary_17'),
                               ('c', 'fixed', 'Q1', 'P_summary_32'), ('d', 'fixed', 'Q2', 'P_summary_35')]:
    data = list(csv.reader(open(os.path.join(archive, 'run0204b', name + '.csv'))))[1:]
    rows += [[panel, beta, tau, r[1], r[-2], r[-1]] for r in data]
write('fig3_lambdas_perms1000_grid50x200.csv', ['panel', 'beta_type', 'tau', 'b', 'GCV', 'LOOCV'], rows)
