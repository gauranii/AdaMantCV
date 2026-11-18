# -----------------------------------------------------------------------------
# File:    checks/search/fig1b_score.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
#          Figure 1B.
# Purpose: Read the ten blue minimum points of the published Figure 1B on its
#          own axes (tick labels give the scales; panel A, reproduced exactly,
#          gives the label offset), and rank the rebuilt minima and the
#          candidates of fig1b_seed_search.R by their largest offset in pt.
# Usage:   python3 checks/search/fig1b_score.py <LiverToxicitySubset1.pdf>
#          <LiverToxicityFull.pdf> results/fig1_minima.csv
#          [checks/search/fig1b_candidates.csv]
# Outputs: prints one line per candidate.
# -----------------------------------------------------------------------------
import csv
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from figure_compare import _fit, marks  # noqa: E402

STEP = 0.1 * math.log(10)        # log-penalty step of 10^seq(9, -5, by = -0.1)
LOWEST = -5 * math.log(10)       # smallest log penalty
BLUE = '0.000 0.000 1.000'


def _marks(pdf):
    m = marks(pdf)
    return [p for p in m['points'] if p[2] == BLUE], m['texts']


def label_offset(pdf_1a, minima_1a):
    """Vertical offset (pt) between a y tick label and its value, from panel A,
    whose ten minima are reproduced exactly (same theme and fonts as B)."""
    pts, texts = _marks(pdf_1a)
    a, b, _ = _fit(minima_1a, [p[1] for p in pts])
    labels = [(float(t[3]), t[1]) for t in texts if t[0] < 30]
    return sum(y - (a + b * v) for v, y in labels) / len(labels)


def published_minima(pdf_1b, offset):
    """The ten published minima of panel B as (log penalty, CV MSE), and the
    pt per unit of each axis, read on the axes fixed by the tick labels."""
    pts, texts = _marks(pdf_1b)
    xt = {t[3]: t[0] for t in texts if t[0] >= 30}
    x_per = (xt['20'] - xt['10']) / 10
    x0 = min(p[0] for p in pts)  # the smallest penalty (five variables sit there)
    yl = [(float(t[3]), t[1] - offset) for t in texts if t[0] < 30]
    a, y_per, _ = _fit([v for v, _ in yl], [y for _, y in yl])
    return [(LOWEST + (p[0] - x0) / x_per, (p[1] - a) / y_per) for p in pts], x_per, y_per


def closeness(candidate, published, x_per, y_per):
    """Largest offset in pt on the published axes, largest log-penalty and
    CV MSE differences, and minima on the published grid penalty."""
    dx = [abs(c[0] - p[0]) for c, p in zip(candidate, published)]
    dy = [abs(c[1] - p[1]) for c, p in zip(candidate, published)]
    same = sum(round(d / STEP) == 0 for d in dx)
    return max(max(dx) * x_per, max(dy) * y_per), max(dx), max(dy), same


def main(pdf_1a, pdf_1b, minima_csv, candidates_csv):
    rows = list(csv.DictReader(open(minima_csv)))
    by_panel = {k: sorted([r for r in rows if r['panel'] == k], key=lambda r: int(r['variable'][8:]))
                for k in 'AB'}
    offset = label_offset(pdf_1a, [float(r['min_cvm']) for r in by_panel['A']])
    published, x_per, y_per = published_minima(pdf_1b, offset)
    groups = {'default (results/fig1_minima.csv)': by_panel['B']}
    if candidates_csv:
        for r in csv.DictReader(open(candidates_csv)):
            groups.setdefault(r['candidate'], []).append(r)
    scored = sorted((closeness([(float(r['log_lambda_min']), float(r['min_cvm'])) for r in g],
                               published, x_per, y_per), name) for name, g in groups.items())
    for (pt, dl, dm, same), name in scored:
        print(f'{name:34} largest offset {pt:5.1f} pt, log penalty within {dl:.2f}, '
              f'CV MSE within {dm:.3f}, {same}/10 on the published penalty')


if __name__ == '__main__':
    main(*sys.argv[1:4], sys.argv[4] if len(sys.argv) > 4 else None)
