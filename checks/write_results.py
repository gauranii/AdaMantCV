# -----------------------------------------------------------------------------
# File:    checks/write_results.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Turn the outputs of check_tables.py and check_figures.py into
#          checks/RESULTS.md.
# Usage:   python3 checks/write_results.py <checks output dir> > checks/RESULTS.md
# -----------------------------------------------------------------------------
import csv
import os
import sys

out = sys.argv[1]
METHODS = ['gcv', 'loo', 'ga_gcv', 'ga_loo']
tables = list(csv.DictReader(open(os.path.join(out, 'tables.csv'))))
figures = list(csv.DictReader(open(os.path.join(out, 'figures.csv'))))


def mark(ok):
    return 'OK' if ok == 'True' else 'FAIL'


md = ['# Check results', '',
      'Written by `checks/run_checks.sh`: the outputs of the scripts in `scripts/`',
      '(run with their defaults) compared with the paper by code.', '',
      '## Tables 1 to 4', '',
      'Every cell is compared with the number printed in the paper',
      '(`checks/paper_values.csv`). Each cell below reads paper / computed, in the',
      'order GCV, LOOCV, gamma GCV, gamma LOOCV.', '']
for n in '1234':
    cells = [r for r in tables if r['table'] == n]
    md += [f"### Table {n}: {sum(r['ok'] == 'True' for r in cells)} of {len(cells)} cells equal", '',
           '| Row | GCV | LOOCV | gamma GCV | gamma LOOCV | Check |', '|---|---|---|---|---|---|']
    for i in range(0, len(cells), 4):
        row = cells[i:i + 4]
        assert [r['method'] for r in row] == METHODS
        md.append(f"| {row[0]['row']} | " + ' | '.join(f"{r['paper']} / {r['computed']}" for r in row) +
                  f" | {mark(str(all(r['ok'] == 'True' for r in row)))} |")
    md.append('')
md += ['## Figures 1 to 3 and 5 to 8', '',
       'Every figure is compared with the published file (`archive/figures/`) by the',
       f"marks in the two PDF files. {sum(r['ok'] == 'True' for r in figures)} of {len(figures)} checks OK.", '',
       '| Figure | Check | Value | Result |', '|---|---|---|---|']
md += [f"| {r['figure']} | {r['check']} | {r['value']} | {mark(r['ok'])} |" for r in figures]
print('\n'.join(md))
