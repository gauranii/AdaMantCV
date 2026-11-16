# -----------------------------------------------------------------------------
# File:    checks/pdf_data.py
# Purpose: Read the marks drawn in an R pdf() file (points, filled rectangles,
#          text) so a figure can be compared with the data behind it.
#          Standard library only.
# Usage:   python3 checks/pdf_data.py figure.pdf   (prints a summary)
# -----------------------------------------------------------------------------
import re
import sys
import zlib
from collections import Counter

NUM = r'(-?\d+(?:\.\d+)?)'


def content(path):
    """Concatenated, decompressed page content streams of an R pdf file."""
    raw = open(path, 'rb').read()
    parts = []
    for m in re.finditer(rb'/Length (\d+)( /Filter /FlateDecode)?\s*>>\s*stream\r?\n', raw):
        body = raw[m.end():m.end() + int(m.group(1))]
        if m.group(2):
            body = zlib.decompress(body)
        parts.append(body.decode('latin1'))
    return '\n'.join(parts)


def marks(path):
    """Return dict with 'points' [(x, y, fill)], 'rects' [(x, y, w, h, fill)],
    'texts' [(x, y, size, string)] in page coordinates (pt)."""
    s = content(path)
    lines = [l.strip() for l in s.split('\n')]
    fill = None
    points, rects, texts = [], [], []
    i = 0
    while i < len(lines):
        l = lines[i]
        m = re.fullmatch(NUM + ' ' + NUM + ' ' + NUM + ' (?:scn|rg)', l)
        if m:
            fill = ' '.join(m.groups())
        # circle: "x y m" followed by four Bezier "c" segments and a fill/stroke
        if l.endswith(' m') and i + 5 < len(lines) and all(lines[i + k].endswith(' c') for k in range(1, 5)):
            x0, y0 = map(float, l.split()[:2])
            cs = [list(map(float, lines[i + k].split()[:6])) for k in range(1, 5)]
            xs = [x0] + [c[4] for c in cs]
            ys = [y0] + [c[5] for c in cs]
            points.append(((max(xs) + min(xs)) / 2, (max(ys) + min(ys)) / 2, fill))
            i += 6
            continue
        m = re.fullmatch(NUM + ' ' + NUM + ' ' + NUM + ' ' + NUM + r' re(?: (f|B))?', l)
        if m and not m.group(5) and i + 1 < len(lines) and lines[i + 1] in ('f', 'B'):
            i += 1  # newer R writes the paint operator on its own line
        elif m and not m.group(5):
            m = None
        if m:
            rects.append(tuple(float(v) for v in m.groups()[:4]) + (fill,))
        # text: "... a 0 0 a x y Tm (string) Tj"; Dingbats "l" is a filled point (old R)
        m = re.search(NUM + r' 0(?:\.00)? -?0(?:\.00)? ' + NUM + ' ' + NUM + ' ' + NUM + r' Tm \((.*)\) Tj', l)
        if m:
            size, x, y, txt = float(m.group(1)), float(m.group(3)), float(m.group(4)), m.group(5)
            if txt == 'l' and '/F1 ' in l:
                points.append((x + 0.4 * size, y + 0.4 * size, fill))
            else:
                texts.append((x, y, size, txt))
        i += 1
    return {'points': points, 'rects': rects, 'texts': texts}


if __name__ == '__main__':
    for p in sys.argv[1:]:
        mk = marks(p)
        print(p, 'points', len(mk['points']), Counter(q[2] for q in mk['points']).most_common(3),
              'rects', len(mk['rects']), 'texts', len(mk['texts']))
