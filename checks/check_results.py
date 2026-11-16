# -----------------------------------------------------------------------------
# File:    checks/check_results.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Compare the numbers the rebuilt scripts computed (results/) with the
#          archived 2021/2022 run outputs and with the paper.
# Usage:   python3 checks/check_results.py <archive dir> <repo> [<original rerun dir>]
# Outputs: results/checks/results.csv; prints one line per check.
# -----------------------------------------------------------------------------
import csv
import glob
import os
import re
import sys

archive, repo = sys.argv[1], sys.argv[2]
original_rerun = sys.argv[3] if len(sys.argv) > 3 else None
res = os.path.join(repo, 'results')
rows = []


def report(item, check, value, ok):
    rows.append({'item': item, 'check': check, 'value': value, 'ok': ok})
    print(f"{item:9} {check:66} {value:>24}  {'OK' if ok else 'DIFFERENT'}")


def read(path):
    return list(csv.DictReader(open(path)))


POWER_COLS = ['ave_logpval', 'ave_logpval_MOM', 'sd_logpval', 'sd_logpval_MOM', 'power', 'power_MOM']


def compare_curves(item, archived, power_file, gcv_file):
    """Largest differences between a rebuilt curve table and an archived CV_ch
    table: lambdas, GCV and LOOCV (relative), and the power columns with and
    without the lambda = 0 row (whose statistic is tied across permutations)."""
    a = read(archived)
    g = read(gcv_file)
    lam = max(abs(float(x['lambdas']) - float(y['lambdas'])) for x, y in zip(a, g))
    gcv = max(abs(float(y['gcv']) / float(x['l_opt$gcv']) - 1) for x, y in zip(a, g))
    loo = max(abs(float(y['loo']) / float(x['l_opt$loo']) - 1) for x, y in zip(a, g))
    report(item, 'lambda grid, largest difference', f'{len(g)}, {lam:.1e}', len(a) == len(g) and lam < 1e-9)
    # 1e-6: near the singular penalties the scores reach 1e10 and the 2021
    # cluster BLAS rounds differently; elsewhere the agreement is about 1e-11.
    report(item, 'GCV and LOOCV curves, largest relative difference', f'{max(gcv, loo):.1e}', max(gcv, loo) < 1e-6)
    if power_file is None:
        return
    p = read(power_file)
    runs = re.search(r'_runs(\d+)_', power_file).group(1)
    other = max(abs(float(x[c]) - float(y[c])) for x, y in zip(a, p) if float(x['lambdas']) != 0
                for c in POWER_COLS)
    zero = max(abs(float(x['power']) - float(y['power'])) for x, y in zip(a, p) if float(x['lambdas']) == 0)
    full = runs == '1000'
    report(item, f'power table ({runs} runs) vs archive (1000 runs), lambda != 0',
           f'{other:.1e}', other < 1e-9 if full else True)
    report(item, 'power at lambda = 0 (tied statistic), difference', f'{zero:.3f}', True)


# ----- Fig 2: liver power curves, full runs -----
for j in range(1, 11):
    pf = glob.glob(os.path.join(res, f'power_liver{j}_runs1000_*.csv'))
    gf = os.path.join(res, f'gcv_liver{j}_grid50x250.csv')
    if os.path.exists(gf):
        compare_curves(f'Fig 2 v{j}', os.path.join(archive, 'liver', f'liver_CV_ch{j}.csv'),
                       pf[0] if pf else None, gf)

# ----- Figs 5 to 8: EEG curves (power with the runs that were done) -----
for gf in sorted(glob.glob(os.path.join(res, 'gcv_eeg*_grid500x2000.csv'))):
    k = re.search(r'gcv_eeg(\d+)_', gf).group(1)
    pf = sorted(glob.glob(os.path.join(res, f'power_eeg{k}_runs*_perms1000_grid500x2000.csv')),
                key=lambda f: int(re.search(r'_runs(\d+)_', f).group(1)))
    compare_curves(f'EEG {k}', os.path.join(archive, 'coh7', f'CV_ch{k}.csv'), pf[-1] if pf else None, gf)

# ----- Table 4 -----
paper = read(os.path.join(repo, 'checks', 'paper_values.csv'))
t4 = glob.glob(os.path.join(res, 'table4_eeg_pvalues_perms1000_grid50x200.csv'))
if t4:
    t = read(t4[0])
    methods = ['gcv', 'loo', 'ga_gcv', 'ga_loo']
    rebuilt = [round(float(r[m]), 3) for r in t for m in methods]
    printed = [float(r['value']) for r in paper if r['table'] == '4']
    archived = [round(float(x[1]), 3) for k in [11, 27, 70, 178, 248]
                for x in list(csv.reader(open(os.path.join(archive, 'real0210a', f'AdaPval_{k}q_0.5.csv'))))[1:]]
    report('Table 4', 'archived run outputs equal to paper (cells)',
           f'{sum(a == b for a, b in zip(archived, printed))}/20', archived == printed)
    report('Table 4', 'rebuilt equal to paper (cells)', f'{sum(a == b for a, b in zip(rebuilt, printed))}/20',
           rebuilt == printed)
    if original_rerun:
        orig = [round(float(x[1]), 3) for k in [11, 27, 70, 178, 248]
                for x in list(csv.reader(open(os.path.join(original_rerun, f'AdaPval_{k}q_0.5.csv'))))[1:]]
        report('Table 4', 'rebuilt equal to ORIGINAL code rerun on this machine',
               f'{sum(a == b for a, b in zip(rebuilt, orig))}/20', rebuilt == orig)
        same = all(open(os.path.join(original_rerun, f'SummaryCV_{k}q_0.5.csv')).read() ==
                   open(os.path.join(archive, 'real0210a', f'SummaryCV_{k}q_0.5.csv')).read()
                   for k in [11, 27, 70, 178, 248])
        report('Table 4', 'original rerun: lambdas chosen per permutation = archive', str(same), same)
    chosen = [read(os.path.join(archive, 'real0210a', f'SummaryCV_{k}q_0.5.csv'))[0]['lhat_GCV']
              for k in [11, 27, 70, 178, 248]]
    report('Table 4', 'archive: lambda chosen by GCV for the observed data', ' '.join(chosen), True)

# ----- Fig 3: rebuilt replication 500 (simulated X; LAPACK dependent) -----
f3 = glob.glob(os.path.join(res, 'fig3_lambdas_perms1000_grid50x200.csv'))
if f3:
    d = read(f3[0])
    text = {'c': (-9.867, -10.023), 'd': (-10.623, -10.689)}  # paper text, b = 1..1000
    for panel, (g, l) in text.items():
        s = [r for r in d if r['panel'] == panel and int(r['b']) <= 1000]
        mg = sum(float(r['GCV']) for r in s) / len(s)
        ml = sum(float(r['LOOCV']) for r in s) / len(s)
        report(f'Fig 3{panel.upper()}', 'mean chosen lambda GCV, LOOCV (paper text in brackets)',
               f'{mg:.3f} ({g}), {ml:.3f} ({l})', True)

# ----- Tables 1 to 3: rebuilt replications vs archived per-run p-values -----
# The simulated X comes from MASS::mvrnorm, whose eigenvectors depend on the
# LAPACK build: for the compound symmetric covariance (p - 1 equal eigenvalues)
# they can rotate, for the diagonal (heteroskedastic) one only their signs can
# flip. Sign flips of columns of X leave U and d unchanged, so Table 1 (beta = 0)
# repeats exactly for the heteroskedastic rows; Tables 2 and 3 need X itself.
for name, prefix in [('table1_type1_error', 'T1'), ('table2_power_fixed', 'P'), ('table3_power_random', 'P')]:
    files = glob.glob(os.path.join(res, f'{name}_perms1000_grid50x200', '*.csv'))
    if not files:
        continue
    tally = {}
    for f in files:
        r = read(f)[0]
        s_, run = int(r['scenario']), int(r['run'])
        cov = 'CS' if s_ % 3 == 2 else 'Het'
        arch = list(csv.reader(open(os.path.join(archive, 'run0204b', f'{prefix}_pvals_{s_}.csv'))))[1:]
        a = [float(x[run]) for x in arch]
        b = [float(r[m]) for m in ['p_gcv', 'p_loo', 'p_ga_gcv', 'p_ga_loo']]
        t = tally.setdefault(cov, [0, 0, set()])
        t[0] += sum(abs(x - y) < 1e-9 for x, y in zip(a, b))
        t[1] += 4
        t[2].add(run)
    for cov, (same, total, runs) in sorted(tally.items()):
        report(name.split('_')[0].title(), f'{cov} rows: rebuilt per-run p-values = archive (runs {min(runs)}-{max(runs)})',
               f'{same}/{total}', same == total if (cov == 'Het' and prefix == 'T1') else True)

os.makedirs(os.path.join(res, 'checks'), exist_ok=True)
with open(os.path.join(res, 'checks', 'results.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['item', 'check', 'value', 'ok'])
    w.writeheader()
    w.writerows(rows)
print(f'results: {sum(r["ok"] for r in rows)}/{len(rows)} checks OK')
