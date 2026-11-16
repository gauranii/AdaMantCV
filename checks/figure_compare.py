# -----------------------------------------------------------------------------
# File:    checks/figure_compare.py
# Purpose: Compare a figure with the data behind it, or two figures with each
#          other, using the marks read by pdf_data.py.
#          - scatter_vs_data: drawn points, in drawing order, against data rows;
#            fits one affine map per axis and reports the largest leftover (pt).
#          - tiles_vs_tiles: heat map tiles, fill colour tile by tile.
#          - bars_vs_bars: histogram bars, heights must be proportional.
# -----------------------------------------------------------------------------
import csv
from pdf_data import marks


def read_csv(path):
    return list(csv.DictReader(open(path)))


def _fit(u, v):
    n = len(u)
    mu, mv = sum(u) / n, sum(v) / n
    sxx = sum((a - mu) ** 2 for a in u)
    b = sum((a - mu) * (c - mv) for a, c in zip(u, v)) / sxx if sxx else 0.0
    a = mv - b * mu
    return a, b, max(abs(a + b * x - y) for x, y in zip(u, v))


def scatter_vs_data(pdf, xs, ys, colour=None):
    """Points drawn in `pdf` (optionally only those with fill `colour`) against
    data (xs, ys) in the same order. Returns (n_drawn, n_data, max_dx_pt, max_dy_pt)."""
    pts = [p for p in marks(pdf)['points'] if colour is None or p[2] == colour]
    if len(pts) != len(xs):
        return len(pts), len(xs), None, None
    _, _, rx = _fit(xs, [p[0] for p in pts])
    _, _, ry = _fit(ys, [p[1] for p in pts])
    return len(pts), len(xs), rx, ry


def _tiles(pdf):
    rects = marks(pdf)['rects']
    w = sorted(r[2] for r in rects)[len(rects) // 2]  # typical tile width
    return [r for r in rects if abs(r[2] - w) < 1e-6]


def _rgb(s):
    return [float(v) for v in s.split()]


def tiles_vs_tiles(pdf_a, pdf_b):
    """Tile count in each file and largest colour channel difference (0..1)."""
    a, b = _tiles(pdf_a), _tiles(pdf_b)
    if len(a) != len(b):
        return len(a), len(b), None
    diff = max(max(abs(x - y) for x, y in zip(_rgb(ra[4]), _rgb(rb[4]))) for ra, rb in zip(a, b))
    return len(a), len(b), diff


def labels(pdf):
    return [t[3] for t in marks(pdf)['texts']]


def _bars(pdf):
    # histogram bars are drawn as "re B"; drop the panel border and background
    rects = marks(pdf)['rects']
    w = sorted(r[2] for r in rects)[len(rects) // 2]
    return [r for r in rects if abs(r[2] - w) < 0.05]


def bars_vs_bars(pdf_a, pdf_b):
    """Bar count in each file, and the spread of height ratios a/b over non-zero
    bars (0 means proportional); also checks zero bars agree."""
    a, b = _bars(pdf_a), _bars(pdf_b)
    if len(a) != len(b):
        return len(a), len(b), None, None
    zero_ok = all((ra[3] == 0) == (rb[3] == 0) for ra, rb in zip(a, b))
    ratios = [ra[3] / rb[3] for ra, rb in zip(a, b) if rb[3] > 0]
    spread = (max(ratios) - min(ratios)) / (sum(ratios) / len(ratios)) if ratios else 0.0
    return len(a), len(b), spread, zero_ok
