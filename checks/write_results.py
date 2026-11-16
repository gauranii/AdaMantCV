# -----------------------------------------------------------------------------
# File:    checks/write_results.py
# Purpose: Turn the check outputs in results/checks/ into checks/RESULTS.md.
# Usage:   python3 checks/write_results.py results/checks > checks/RESULTS.md
# -----------------------------------------------------------------------------
import csv
import datetime
import os
import sys

out = sys.argv[1]


def table(path, columns):
    if not os.path.exists(path):
        return ['(not run)', '']
    rows = list(csv.DictReader(open(path)))
    lines = ['| ' + ' | '.join(columns) + ' |', '|' + '---|' * len(columns)]
    for r in rows:
        cells = [str(r[c]).replace('|', '/') for c in columns]
        cells[-1] = 'OK' if cells[-1] in ('True', 'TRUE') else 'differs'
        lines.append('| ' + ' | '.join(cells) + ' |')
    ok = sum(r[columns[-1]] in ('True', 'TRUE') for r in rows)
    return [f'{ok} of {len(rows)} checks OK.', ''] + lines + ['']


def text(path):
    return open(path).read().strip().split('\n') if os.path.exists(path) else ['(not run)']


md = [f'# Check results ({datetime.date.today()})', '',
      'Written by `checks/run_checks.sh`. "differs" is not always a failure: the',
      'README explains each expected difference.', '',
      '## 1. Archived run outputs against the paper', '']
md += ['```'] + text(os.path.join(out, 'archive_vs_paper.txt')) + ['```', '']
md += ['## 2. Rebuilt functions against the original functions (same inputs and seeds)', '']
md += table(os.path.join(out, 'equivalence.csv'), ['what', 'n', 'max_abs_diff', 'tolerance', 'ok'])
md += ['## 3. Published figures against the same figures redrawn from the archived data', '']
md += table(os.path.join(out, 'figures.csv'), ['figure', 'check', 'value', 'ok'])
md += ['## 4. Rebuilt results against the archive and the paper', '']
md += table(os.path.join(out, 'results.csv'), ['item', 'check', 'value', 'ok'])
print('\n'.join(md))
