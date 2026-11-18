# -----------------------------------------------------------------------------
# File:    checks/search/sim_lambda0_count.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
#          Tables 1 to 3.
# Purpose: Count, over the 500 archived runs of each setting, how often GCV
#          chose lambda = 0 (where every permutation ties and the p-value is
#          set by rounding) or |lambda| = sqrt(machine epsilon) for the
#          observed data.
# Usage:   python3 checks/search/sim_lambda0_count.py <Run0204b folder>
#          (the folder with Scenario2 ... Scenario36, from the external drive)
# Outputs: prints one line per table and dimension.
# -----------------------------------------------------------------------------
import csv
import os
import sys


def count(runs):
    """{(table, p): [runs choosing 0, runs choosing +-sqrt(eps), runs]}"""
    counts = {}
    for s in range(2, 37):
        if s % 3 == 1:          # identity covariance: not in the paper
            continue
        for prefix in ['P', 'T1']:
            if prefix == 'T1' and s > 18:   # Table 1 uses scenarios 2 to 18
                continue
            folder = os.path.join(runs, f'Scenario{s}')
            lambdas = [float(r[1]) for r in
                       list(csv.reader(open(os.path.join(folder, f'{prefix}_lambdas_{s}.csv'))))[1:]]
            chosen = [lambdas[int(float(r[1])) - 1] for r in
                      list(csv.reader(open(os.path.join(folder, f'{prefix}_lopt_{s}.csv'))))[1:]]
            table = 'Table 1' if prefix == 'T1' else ('Table 3' if s <= 18 else 'Table 2')
            p = 500 if (s - 1) % 18 < 9 else 1000
            c = counts.setdefault((table, p), [0, 0, 0])
            c[0] += sum(v == 0 for v in chosen)
            c[1] += sum(0 < abs(v) < 1e-7 for v in chosen)
            c[2] += len(chosen)
    return counts


if __name__ == '__main__':
    for (table, p), (zero, eps, total) in sorted(count(sys.argv[1]).items()):
        print(f'{table}, p = {p}: observed data chose lambda = 0 in {zero} of {total} runs, '
              f'|lambda| = sqrt(eps) in {eps}')
