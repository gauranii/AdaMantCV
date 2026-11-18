# -----------------------------------------------------------------------------
# File:    checks/check_readme.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Build the two side-by-side tables of the README (paper, best-fit
#          default, correct option) from the result files and the check
#          outputs, and check that the README holds them exactly.
# Usage:   python3 checks/check_readme.py <repo> <LiverToxicityFull.pdf> [--print]
#          ADAMANT_RUNS: the archived Run0204b folder (Scenario2 ... 36), for
#          the third table; ADAMANT_REAL: the archived Real0210a folder
#          (all 250 triples), for two numbers in the prose
#          (after the scripts and checks/run_checks.sh have run)
#          --print writes the two tables instead of checking them.
# Outputs: prints OK, or each table that differs; exit status 1 on a miss.
# -----------------------------------------------------------------------------
import csv
import glob
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'search'))
from fig1b_score import closeness, label_offset, published_minima  # noqa: E402
from sim_lambda0_count import count  # noqa: E402

repo, published_1b = sys.argv[1], sys.argv[2]
res = os.path.join(repo, 'results')
checks = os.path.join(res, 'checks')
METHODS = ['gcv', 'loo', 'ga_gcv', 'ga_loo']
TRIPLES = [11, 27, 70, 178, 248]


def read(path):
    return list(csv.DictReader(open(path)))


def value(rows, item, check_start):
    hits = [r['value'] for r in rows if r.get('item', r.get('figure')) == item
            and r['check'].startswith(check_start)]
    return hits[0]


def cells(values):
    return ' / '.join(f'{v:.3f}' for v in values)


# ----- Table 4, side by side -----
paper = read(os.path.join(repo, 'checks', 'paper_values.csv'))
printed = [float(r['value']) for r in paper if r['table'] == '4']
versions = {v: read(os.path.join(res, f'table4_eeg_pvalues_{v}_perms1000_grid50x200.csv'))
            for v in ['best_fit', 'original', 'paper']}
t4 = ['| Triple | Paper | Best fit (default) | Original code, this machine | Correct (`--version=paper`) |',
      '|---|---|---|---|---|']
for i, k in enumerate(TRIPLES):
    row = versions['best_fit'][i]
    assert int(row['triple']) == k
    label = f"{row['channel1']} - {row['channel2']}, {row['band']}"
    t4.append(f'| {label} | {cells(printed[4 * i:4 * i + 4])} | ' +
              ' | '.join(cells([float(versions[v][i][m]) for m in METHODS])
                         for v in ['best_fit', 'original', 'paper']) + ' |')
t4.append('')
t4.append('Each cell: GCV / LOOCV / gamma GCV / gamma LOOCV.')

# ----- Per table and figure -----
avp = read(os.path.join(checks, 'archive_vs_paper.csv'))
results = read(os.path.join(checks, 'results.csv'))
figures = read(os.path.join(checks, 'figures.csv'))


def archive_cells(table):
    rows = [r for r in avp if r['table'] == str(table)]
    return f"{sum(r['match'] == 'TRUE' for r in rows)}/{len(rows)}"


def per_run(table):
    het = value(results, f'Table{table}', 'Het rows')
    cs = value(results, f'Table{table}', 'CS rows')
    return f'heteroskedastic {het}, compound symmetric {cs}'


def paper_runs(name):
    files = glob.glob(os.path.join(res, f'{name}_paper_perms1000_grid50x200', '*.csv'))
    if not files:
        return 'not run yet (job N14)'
    p = [float(r[c]) for f in files for r in read(f) for c in ['p_gcv', 'p_loo', 'p_ga_gcv', 'p_ga_loo']]
    settings = sorted({int(r['scenario']) for f in files for r in read(f)})
    runs = sorted({int(r['run']) for f in files for r in read(f)})
    spread = 'every p-value 1' if min(p) == 1 else f'p-values {min(p):.3f} to {max(p):.3f}'
    return (f'runs {runs[0]} to {runs[-1]} of settings {", ".join(map(str, settings))} '
            f'(original scenario numbers): {spread}')


def fig1b():
    m = read(os.path.join(res, 'fig1_minima.csv'))
    a = sorted([r for r in m if r['panel'] == 'A'], key=lambda r: int(r['variable'][8:]))
    b = sorted([r for r in m if r['panel'] == 'B'], key=lambda r: int(r['variable'][8:]))
    offset = label_offset(os.path.join(os.path.dirname(published_1b), 'LiverToxicitySubset1.pdf'),
                          [float(r['min_cvm']) for r in a])
    published, x_per, y_per = published_minima(published_1b, offset)
    pt, dl, dm, _ = closeness([(float(r['log_lambda_min']), float(r['min_cvm'])) for r in b],
                              published, x_per, y_per)
    return (f'closest of 17 session states tried: largest offset {pt:.1f} pt, log penalty '
            f'within {dl:.2f}, CV MSE within {dm:.3f} (folds lost, see note 10)')


def fig2_paper():
    pairs = []
    for j in range(1, 11):
        a = glob.glob(os.path.join(res, f'power_liver{j}_runs1000_*.csv'))
        b = glob.glob(os.path.join(res, f'power_liver{j}_paper_runs1000_*.csv'))
        if a and b:
            pairs.append((read(a[0]), read(b[0])))
    if len(pairs) < 10:
        return 'rerun not finished'
    other = max(abs(float(x[c]) - float(y[c])) for a, b in pairs for x, y in zip(a, b)
                if float(x['lambdas']) != 0 for c in ['power', 'power_MOM'])
    zero = [float(y[c]) for a, b in pairs for x, y in zip(a, b)
            if float(x['lambdas']) == 0 for c in ['power', 'power_MOM']]
    return (f'power at lambda != 0 equal to the default (largest difference {other:.1g}); '
            f'at lambda = 0 power {max(zero):g} for all 10 variables')


def fig_ok(prefix):
    rows = [r for r in figures if r['figure'].startswith(prefix) and 'report only' not in r['check']]
    return f"{sum(r['ok'] == 'True' for r in rows)}/{len(rows)} checks"


t4_best = value(results, 'Table 4', 'best_fit: cells equal')
t4_paper = [float(versions['paper'][i][m]) for i in range(5) for m in METHODS]
summary = [
    '| Item | Best fit (default) against the paper | Correct option (`--method=paper`) |',
    '|---|---|---|',
    f'| Table 1 | archived runs give the paper: {archive_cells(1)} cells; runs 1 to 10 rerun equal to the archive: {per_run(1)} | {paper_runs("table1_type1_error")} |',
    f'| Table 2 | archived runs give the paper: {archive_cells(2)} cells; runs 1 to 10 rerun equal to the archive: {per_run(2)} | {paper_runs("table2_power_fixed")} |',
    f'| Table 3 | archived runs give the paper: {archive_cells(3)} cells; runs 1 to 10 rerun equal to the archive: {per_run(3)} | {paper_runs("table3_power_random")} |',
    f'| Table 4 | {t4_best} cells equal to the paper | {"all 20 p-values 1" if all(v == 1 for v in t4_paper) else "see table above"} |',
    f'| Figure 1A | {value(figures, "Fig 1A", "blue minimum").split(", ")[0]} minima within {value(figures, "Fig 1A", "blue minimum").split(", ")[1]} pt of the published points | same code |',
    f'| Figure 1B | {fig1b()} | same code |',
    f'| Figure 2 | drawn from the archived runs: {fig_ok("Fig 2")} OK | {fig2_paper()} |',
    f'| Figure 3 | drawn from the archived runs: {fig_ok("Fig 3")} OK; rerun mean chosen penalty, GCV and LOOCV (paper text in brackets): C {value(results, "Fig 3C", "mean")}, D {value(results, "Fig 3D", "mean")} | not run yet (job N14) |',
    f'| Figures 5, 6 | drawn from the archived runs: Figure 5 {fig_ok("Fig 5")} OK, Figure 6 {fig_ok("Fig 6")} OK | not run |',
    f'| Figures 7, 8 | drawn from the archived runs: Figure 7 {fig_ok("Fig 7")} OK, Figure 8 {fig_ok("Fig 8")} OK | not run |',
]

tables = ['\n'.join(t4), '\n'.join(summary)]
if os.environ.get('ADAMANT_RUNS'):
    runs = ['| Table, p | Observed data chose lambda = 0 | Chose lambda = +-sqrt(eps) | Runs |',
            '|---|---|---|---|']
    for (table, p), (zero, eps, total) in sorted(count(os.environ['ADAMANT_RUNS']).items()):
        runs.append(f'| {table}, p = {p} | {zero} | {eps} | {total} |')
    tables.append('\n'.join(runs))
# Numbers in the prose, each as the sentence that holds it
p1000 = [float(r['value']) for r in paper if r['table'] in ('2', '3') and r['p'] == '1000']
n_states = len({r['candidate'] for r in read(os.path.join(repo, 'checks', 'search', 'fig1b_candidates.csv'))})
sentences = [f'the rows (power {min(p1000):.3f} to {max(p1000):.3f} in the paper) are robust',
             f'`checks/search/fig1b_seed_search.R` tried {n_states} seeds and session states',
             f'The default now uses the closest of {n_states} session states tried']
if os.environ.get('ADAMANT_REAL'):
    real = os.environ['ADAMANT_REAL']
    p_obs = {k: float(list(csv.reader(open(os.path.join(real, f'AdaPval_{k}q_0.5.csv'))))[1][1])
             for k in range(1, 251)}
    zero = sum(read(os.path.join(real, f'SummaryCV_{k}q_0.5.csv'))[0]['lhat_GCV'] == '0'
               for k in range(1, 251))
    worst = max(sorted(p_obs.values()).index(p_obs[k]) + 1 for k in TRIPLES)
    sentences += [f'choose lambda = 0 for the observed data in {zero} of the 250 triples',
                  f'five triples of Table 4 are five of the {worst} smallest p-values of those 250']

if '--print' in sys.argv:
    print('\n'.join(sentences) + '\n')
    print('\n\n'.join(tables))
    sys.exit(0)
readme = open(os.path.join(repo, 'README.md')).read()
missing = [t for t in tables if t not in readme]
flat = ' '.join(readme.split())
missing += [t for t in sentences if t not in flat]
for t in missing:
    print('README differs from the results; expected:\n' + t + '\n')
print(f'README side-by-side tables and numbers in prose ({len(tables) + len(sentences)}): ' +
      ('OK' if not missing else f'{len(missing)} differ'))
sys.exit(1 if missing else 0)
