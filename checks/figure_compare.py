# -----------------------------------------------------------------------------
# File:    checks/figure_compare.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Compare two figures by the marks read by pdf_data.py.
#          - _fit: one affine map between two coordinate lists and the largest
#            leftover (pt), for drawn points.
#          - tiles_vs_tiles: heat map tiles, fill colour tile by tile.
#          - bars_vs_bars: histogram bars, heights must be proportional.
# -----------------------------------------------------------------------------
from pdf_data import marks


def _fit(u, v):
    n = len(u)
    mu, mv = sum(u) / n, sum(v) / n
    sxx = sum((a - mu) ** 2 for a in u)
    b = sum((a - mu) * (c - mv) for a, c in zip(u, v)) / sxx if sxx else 0.0
    a = mv - b * mu
    return a, b, max(abs(a + b * x - y) for x, y in zip(u, v))


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
