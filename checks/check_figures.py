# -----------------------------------------------------------------------------
# File:    checks/check_figures.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Compare every figure drawn by the scripts with the published file
#          by the marks in the two PDF files: tile colours and row labels of
#          the heat maps (Figures 2, 5, 6), bar heights of the histograms
#          (Figure 3), point positions of the scatter plots (Figures 7, 8) and
#          the minimum points of Figure 1.
# Usage:   python3 checks/check_figures.py <published dir> <figures dir>
#          <results dir> <checks output dir>
# Outputs: <checks output dir>/figures.csv; prints one line per check; exit
#          status 1 unless every check is OK.
# -----------------------------------------------------------------------------
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figure_compare import _fit, bars_vs_bars, labels, marks, tiles_vs_tiles  # noqa: E402

published, drawn, results, out = sys.argv[1:5]
rows = []


def report(figure, check, value, ok):
    rows.append({'figure': figure, 'check': check, 'value': value, 'ok': ok})
    print(f"{figure:10} {check:58} {value:>28}  {'OK' if ok else 'FAIL'}")


def row_labels(pdf):
    return [t.lower() for t in labels(pdf) if any(c.isalpha() for c in t) and t != 'lambda']


# A drawn figure has the name of its published file.
names = sorted(os.listdir(published))


def panels(prefix, word=''):
    """(label, file name) of the published files of a figure, by panel."""
    return [(f'Fig {name[3:5]}', name) for name in names
            if name.startswith(prefix) and word in name]


def both(name):
    return os.path.join(published, name), os.path.join(drawn, name)


for fig, name in panels('Fig2') + panels('Fig5') + panels('Fig6'):
    a, b = both(name)
    na, nb, diff = tiles_vs_tiles(a, b)
    report(fig, 'tiles: count, largest colour difference (0-1)', f'{na}/{nb}, {diff}', na == nb and diff == 0)
    same = row_labels(a) == row_labels(b)
    report(fig, 'row labels equal (ignoring case)', str(same), same)

for fig, name in panels('Fig3'):
    na, nb, spread, zeros = bars_vs_bars(*both(name))
    report(fig, 'bars: count, spread of height ratios, same empty bins',
           f'{na}/{nb}, {spread:.4f}, {zeros}', na == nb and zeros and spread < 0.01)

for fig, name in panels('Fig7') + panels('Fig8'):
    a, b = both(name)
    pa, pb = marks(a)['points'], marks(b)['points']
    panel = name.split('_', 2)[2][:-4]  # for example power-decreasing_F3-FT8-beta
    if len(pa) != len(pb):
        report(fig, f'{panel}: points', f'{len(pa)}/{len(pb)}', False)
        continue
    rx = _fit([p[0] for p in pb], [p[0] for p in pa])[2]
    ry = _fit([p[1] for p in pb], [p[1] for p in pa])[2]
    colours = [p[2] for p in pa] == [p[2] for p in pb]
    report(fig, f'{panel}: points, largest offset after scaling (pt)', f'{len(pa)}, {max(rx, ry):.3f}',
           max(rx, ry) < 0.05 and colours)

# Fig 1: the ten blue minimum points against the computed minima
minima = list(csv.DictReader(open(os.path.join(results, 'fig1_minima.csv'))))
for fig, name in panels('Fig1'):
    pts = [p for p in marks(os.path.join(published, name))['points'] if p[2] == '0.000 0.000 1.000']
    m = [r for r in minima if r['panel'] == fig[-1]]
    m.sort(key=lambda r: int(''.join(c for c in r['variable'] if c.isdigit())))
    rx = _fit([float(r['log_lambda_min']) for r in m], [p[0] for p in pts])[2]
    ry = _fit([float(r['min_cvm']) for r in m], [p[1] for p in pts])[2]
    report(fig, 'minimum points vs computed minima, largest offset (pt)', f'{len(pts)}, {max(rx, ry):.3f}',
           max(rx, ry) < 0.05)

os.makedirs(out, exist_ok=True)
with open(os.path.join(out, 'figures.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['figure', 'check', 'value', 'ok'])
    w.writeheader()
    w.writerows(rows)
print(f'figures: {sum(r["ok"] for r in rows)}/{len(rows)} checks OK')
sys.exit(0 if all(r['ok'] for r in rows) else 1)
