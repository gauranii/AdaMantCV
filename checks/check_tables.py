# -----------------------------------------------------------------------------
# File:    checks/check_tables.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Compare every cell of Tables 1 to 4 written by the table scripts
#          with the number printed in the paper (three decimals).
# Usage:   python3 checks/check_tables.py <results dir> <checks output dir>
# Inputs:  checks/paper_values.csv, <results dir>/Table<N>_*.csv.
# Outputs: <checks output dir>/tables.csv; prints one line per table; exit
#          status 1 unless every cell is equal.
# -----------------------------------------------------------------------------
import csv
import os
import sys

results, out = sys.argv[1:3]
here = os.path.dirname(os.path.abspath(__file__))
FILES = {'1': 'Table1_type1-error', '2': 'Table2_power_fixed-effects',
         '3': 'Table3_power_random-effects', '4': 'Table4_EEG_pvalues'}
COLUMN = {'gcv': 'gcv', 'loo': 'loocv', 'ga_gcv': 'ga_gcv', 'ga_loo': 'ga_loocv'}
COVARIANCE = {'CS': 'CS', 'Het': 'Heteroskedastic'}


def read(path):
    return list(csv.DictReader(open(path)))


tables = {n: read(os.path.join(results, name + '.csv')) for n, name in FILES.items()}
rows = []
for cell in read(os.path.join(here, 'paper_values.csv')):
    row = tables[cell['table']][int(cell['row']) - 1]
    if cell['table'] == '4':
        computed = row[cell['method']]
        label = f"{row['channel1']} - {row['channel2']}, {row['band']}"
    else:
        assert (row['covariance'], row['p'], row['tau']) == \
            (COVARIANCE[cell['covariance']], cell['p'], cell['tau'])
        computed = row[COLUMN[cell['method']]]
        label = f"{row['covariance']}, p = {row['p']}, tau = {row['tau']}"
    computed = f'{float(computed):.3f}'
    rows.append({'table': cell['table'], 'row': label, 'method': cell['method'],
                 'paper': cell['value'], 'computed': computed, 'ok': computed == cell['value']})

os.makedirs(out, exist_ok=True)
with open(os.path.join(out, 'tables.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
for n in FILES:
    cells = [r for r in rows if r['table'] == n]
    print(f"Table {n}: {sum(r['ok'] for r in cells)}/{len(cells)} cells equal to the paper")
sys.exit(0 if all(r['ok'] for r in rows) else 1)
