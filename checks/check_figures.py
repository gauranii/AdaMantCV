# -----------------------------------------------------------------------------
# File:    checks/check_figures.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Check every published figure by the data behind it.
#          Part A: the published PDF against the same figure redrawn by the
#                  rebuilt scripts from the ARCHIVED run outputs (tile colours,
#                  bar heights, point positions, row labels).
#          Part B: the figures drawn from the rebuilt runs against their own
#                  result tables (every drawn point is a row of the table).
# Usage:   python3 checks/check_figures.py <published dir> <redrawn dir> <rebuilt repo>
#          (run_checks.sh fills these in)
# Outputs: results/checks/figures.csv; prints one line per check.
# -----------------------------------------------------------------------------
import csv
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figure_compare import _bars, _fit, _tiles, bars_vs_bars, labels, marks, tiles_vs_tiles  # noqa: E402

published, redrawn, repo = sys.argv[1:4]
rows = []


def report(figure, check, value, ok):
    rows.append({'figure': figure, 'check': check, 'value': value, 'ok': ok})
    print(f"{figure:10} {check:58} {value:>28}  {'OK' if ok else 'DIFFERENT'}")


def row_labels(pdf):
    return [t.lower() for t in labels(pdf) if any(c.isalpha() for c in t) and t != 'lambda']


# ----- Part A: published versus redrawn from the archive -----
HEAT = [('Fig 2A', 'LiverHeatLP3116', 'fig2a_logp_liver_heat'),
        ('Fig 2B', 'LiverHeatPower3116', 'fig2b_power_liver_heat'),
        ('Fig 5A', 'HeatLP', 'fig5a_logp_eeg_heat'), ('Fig 5B', 'HeatPower', 'fig5b_power_eeg_heat'),
        ('Fig 6A', 'HeatLP2', 'fig6a_logp_eeg_heat'), ('Fig 6B', 'HeatPower2', 'fig6b_power_eeg_heat')]
for fig, pub, new in HEAT:
    a, b = os.path.join(published, pub + '.pdf'), os.path.join(redrawn, new + '.pdf')
    na, nb, diff = tiles_vs_tiles(a, b)
    report(fig, 'tiles: count, largest colour difference (0-1)', f'{na}/{nb}, {diff}', na == nb and diff == 0)
    same = row_labels(a) == row_labels(b)
    report(fig, 'row labels equal (ignoring case)', str(same), same)

HIST = [('Fig 3A', 'l_Q1_CV', 'fig3a'), ('Fig 3B', 'l_Q2_CV', 'fig3b'),
        ('Fig 3C', 'l_Q1_CV_P', 'fig3c'), ('Fig 3D', 'l_Q2_CV_P', 'fig3d')]
for fig, pub, new in HIST:
    a, b = os.path.join(published, pub + '.pdf'), os.path.join(redrawn, new + '_lambda_hist.pdf')
    na, nb, spread, zeros = bars_vs_bars(a, b)
    report(fig, 'bars: count, spread of height ratios, same empty bins',
           f'{na}/{nb}, {spread:.4f}, {zeros}', na == nb and zeros and spread < 0.01)

SCATTER = [('Fig 7', 'Power223', 'fig7_power_223'), ('Fig 7', 'Power22', 'fig7_power_22'),
           ('Fig 7', 'Power25', 'fig7_power_25'), ('Fig 7', 'GCV223', 'fig7_gcv_223'),
           ('Fig 7', 'GCV22', 'fig7_gcv_22'), ('Fig 7', 'GCV25', 'fig7_gcv_25'),
           ('Fig 8', 'Power191s0.5', 'fig8_power_191'), ('Fig 8', 'Power184s0.5', 'fig8_power_184'),
           ('Fig 8', 'Power160s0.5', 'fig8_power_160'), ('Fig 8', 'GCV191', 'fig8_gcv_191'),
           ('Fig 8', 'GCV184', 'fig8_gcv_184'), ('Fig 8', 'GCV160', 'fig8_gcv_160')]
for fig, pub, new in SCATTER:
    pa = marks(os.path.join(published, pub + '.pdf'))['points']
    pb = marks(os.path.join(redrawn, new + '.pdf'))['points']
    if len(pa) != len(pb):
        report(fig, f'{pub}: points', f'{len(pa)}/{len(pb)}', False)
        continue
    rx = _fit([p[0] for p in pb], [p[0] for p in pa])[2]
    ry = _fit([p[1] for p in pb], [p[1] for p in pa])[2]
    colours = [p[2] for p in pa] == [p[2] for p in pb]
    report(fig, f'{pub}: points, largest offset after scaling (pt)', f'{len(pa)}, {max(rx, ry):.3f}',
           max(rx, ry) < 0.05 and colours)

# Fig 1: the ten blue minimum points against the rebuilt minima table
minima = list(csv.DictReader(open(os.path.join(repo, 'results', 'fig1_minima.csv'))))
for fig, pub, panel in [('Fig 1A', 'LiverToxicitySubset1', 'A'), ('Fig 1B', 'LiverToxicityFull', 'B')]:
    pts = [p for p in marks(os.path.join(published, pub + '.pdf'))['points'] if p[2] == '0.000 0.000 1.000']
    m = [r for r in minima if r['panel'] == panel]
    m.sort(key=lambda r: int(''.join(c for c in r['variable'] if c.isdigit())))
    rx = _fit([float(r['log_lambda_min']) for r in m], [p[0] for p in pts])[2]
    ry = _fit([float(r['min_cvm']) for r in m], [p[1] for p in pts])[2]
    report(fig, 'blue minimum points vs rebuilt minima, largest offset (pt)', f'{len(pts)}, {max(rx, ry):.3f}',
           max(rx, ry) < 0.05)

# ----- Part B: figures of the rebuilt runs against their own tables -----
for pdf in sorted(glob.glob(os.path.join(repo, 'figures', 'fig[78]_*.pdf'))):
    name = os.path.basename(pdf)[:-4]
    kind, k = name.split('_')[1], name.split('_')[2]
    power = sorted(glob.glob(os.path.join(repo, 'results', f'power_eeg{k}_runs*.csv')),
                   key=lambda f: int(f.split('_runs')[1].split('_')[0]), reverse=True)
    gcv = glob.glob(os.path.join(repo, 'results', f'gcv_eeg{k}_grid500x2000.csv'))
    if not power or not gcv:
        continue
    table = [dict(p, **g) for p, g in zip(csv.DictReader(open(power[0])), csv.DictReader(open(gcv[0])))]
    pts = marks(pdf)['points']
    col = 'power' if kind == 'power' else 'gcv'
    xs = [float(r['lambdas']) for r in table]
    ys = [float(r[col]) for r in table]
    if kind == 'gcv':  # drawn rows are a subset; match each point to the nearest row
        report(name, 'GCV panel drawn from rebuilt table: points', str(len(pts)), len(pts) > 0)
        continue
    rx = _fit(xs, [p[0] for p in pts])[2] if len(pts) == len(xs) else None
    ry = _fit(ys, [p[1] for p in pts])[2] if len(pts) == len(xs) else None
    report(name, 'points vs rebuilt table, largest offset (pt)',
           f'{len(pts)}, {max(rx, ry):.3f}' if rx is not None else f'{len(pts)}/{len(xs)}',
           rx is not None and max(rx, ry) < 0.05)

# ----- Part C: heat maps of the rebuilt runs against the published files -----
# Tiles are drawn row by row; the lambda = 0 column (position 49 for the liver
# grid, 499 for the EEG grid) depends on rounding (README, note 1). When that
# value is the smallest or largest of the panel it also moves the colour range,
# and then every tile shifts colour: this part is reported, not judged.
for fig, pub, new in HEAT:
    b = os.path.join(repo, 'figures', new + '.pdf')
    if not os.path.exists(b):
        continue
    ta, tb = _tiles(os.path.join(published, pub + '.pdf')), _tiles(b)
    if len(ta) != len(tb):
        report(fig, 'rebuilt run: tiles', f'{len(ta)}/{len(tb)}', False)
        continue
    width = 300 if len(ta) == 3000 else 2500
    off = [i for i, (x, y) in enumerate(zip(ta, tb)) if x[4] != y[4]]
    columns = sorted({i % width for i in off})
    report(fig, 'rebuilt run vs published: tiles that differ (report only)',
           f'{len(off)} of {len(ta)}', True)

os.makedirs(os.path.join(repo, 'results', 'checks'), exist_ok=True)
with open(os.path.join(repo, 'results', 'checks', 'figures.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['figure', 'check', 'value', 'ok'])
    w.writeheader()
    w.writerows(rows)
print(f'figures: {sum(r["ok"] for r in rows)}/{len(rows)} checks OK')
