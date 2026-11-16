# -----------------------------------------------------------------------------
# File:    checks/extract_paper_values.py
# Purpose: Pull the numbers of Tables 1 to 4 out of the accepted manuscript
#          (frontiers.tex) into checks/paper_values.csv.
# Usage:   python3 checks/extract_paper_values.py path/to/frontiers.tex > checks/paper_values.csv
# -----------------------------------------------------------------------------
import re
import sys

tex = open(sys.argv[1]).read()
methods = ['gcv', 'loo', 'ga_gcv', 'ga_loo']
print('table,row,covariance,p,tau,band,channels,method,value')
for table, label in [(1, 'tab:TypeIFE'), (2, 'tab:PowerFE'), (3, 'tab:PowerRE'), (4, 'tab:SigEEG')]:
    body = tex[:tex.index('\\label{' + label + '}')]
    body = body[body.rindex('\\begin{tabular}'):]
    rows = [l for l in body.split('\n') if len(re.findall(r'\d\.\d{3}', l)) == 4]
    cov = p = ''
    for i, l in enumerate(rows, 1):
        cells = [c.strip() for c in l.replace('\\\\', '').split('&')]
        vals = re.findall(r'\d\.\d{3}', l)
        if table < 4:
            cov = {'Compound Symmetric': 'CS', 'Heteroskedastic': 'Het'}.get(cells[0], cov)
            p = cells[1] or p
            tau = {'None': 'none', '$Q_1$': 'Q1', '$Q_2$': 'Q2'}[cells[2]]
            key = [cov, p, tau, '', '']
        else:
            key = ['', '', '', cells[0], cells[1]]
        for m, v in zip(methods, vals):
            print(','.join(map(str, [table, i] + key + [m, v])))
