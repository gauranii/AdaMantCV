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


HEAT = [('Fig 2A', 'LiverHeatLP3116', 'fig2a_logp_liver_heat'),
        ('Fig 2B', 'LiverHeatPower3116', 'fig2b_power_liver_heat'),
        ('Fig 5A', 'HeatLP', 'fig5a_logp_eeg_heat'), ('Fig 5B', 'HeatPower', 'fig5b_power_eeg_heat'),
        ('Fig 6A', 'HeatLP2', 'fig6a_logp_eeg_heat'), ('Fig 6B', 'HeatPower2', 'fig6b_power_eeg_heat')]
for fig, pub, new in HEAT:
    a, b = os.path.join(published, pub + '.pdf'), os.path.join(drawn, new + '.pdf')
    na, nb, diff = tiles_vs_tiles(a, b)
    report(fig, 'tiles: count, largest colour difference (0-1)', f'{na}/{nb}, {diff}', na == nb and diff == 0)
    same = row_labels(a) == row_labels(b)
    report(fig, 'row labels equal (ignoring case)', str(same), same)

HIST = [('Fig 3A', 'l_Q1_CV', 'fig3a'), ('Fig 3B', 'l_Q2_CV', 'fig3b'),
        ('Fig 3C', 'l_Q1_CV_P', 'fig3c'), ('Fig 3D', 'l_Q2_CV_P', 'fig3d')]
for fig, pub, new in HIST:
    a, b = os.path.join(published, pub + '.pdf'), os.path.join(drawn, new + '_lambda_hist.pdf')
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
    pb = marks(os.path.join(drawn, new + '.pdf'))['points']
    if len(pa) != len(pb):
        report(fig, f'{pub}: points', f'{len(pa)}/{len(pb)}', False)
        continue
    rx = _fit([p[0] for p in pb], [p[0] for p in pa])[2]
    ry = _fit([p[1] for p in pb], [p[1] for p in pa])[2]
    colours = [p[2] for p in pa] == [p[2] for p in pb]
    report(fig, f'{pub}: points, largest offset after scaling (pt)', f'{len(pa)}, {max(rx, ry):.3f}',
           max(rx, ry) < 0.05 and colours)

# Fig 1: the ten blue minimum points against the computed minima
minima = list(csv.DictReader(open(os.path.join(results, 'fig1_minima.csv'))))
for fig, pub, panel in [('Fig 1A', 'LiverToxicitySubset1', 'A'), ('Fig 1B', 'LiverToxicityFull', 'B')]:
    pts = [p for p in marks(os.path.join(published, pub + '.pdf'))['points'] if p[2] == '0.000 0.000 1.000']
    m = [r for r in minima if r['panel'] == panel]
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
