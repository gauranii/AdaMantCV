# -----------------------------------------------------------------------------
# File:    checks/check_python.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Compare the outputs of the R scripts with those of the Python
#          scripts, output by output: every cell of every csv under results/
#          (tables, curves, per-replication p-values) and the numbers every
#          figure plots (plot_data/, one row per plotted mark).
# Usage:   python3 checks/check_python.py archive|rerun <R output dir>
#          <Python output dir> <comparison csv>
# Inputs:  <dir>/results/, <dir>/plot_data/, <dir>/logs/<script>.done (time
#          peak memory and cores) of both implementations, written by
#          checks/check_python.sh.
# Outputs: the comparison csv, one row per output of the paper: files and
#          values compared, how many values are not identical, the largest
#          absolute and relative difference, the tolerance, the verdict and
#          the run times. Prints the same; exit status 1 unless every output
#          is identical or within its tolerance.
#          Tolerance (relative, |a - b| <= tol * max(|a|, |b|, 1)):
#          archive mode 0 for results/ and 1e-9 for plotted numbers;
#          rerun mode 1e-8 for both. Text and the layout of the files must
#          be the same in every mode.
# -----------------------------------------------------------------------------
import csv
import glob
import gzip
import math
import os
import re
import sys

mode, r_dir, py_dir, out_file = sys.argv[1:5]
TOLERANCE = {"archive": {"results": 0.0, "plot_data": 1e-9},
             "rerun": {"results": 1e-8, "plot_data": 1e-8}}[mode]
SETTINGS = {"archive": {"fig1": "from the data"}, "rerun": {"table1": "2 replications", "table2": "2 replications",
                                     "table3": "2 replications", "fig5": "3 replications",
                                     "fig6": "3 replications"}}[mode]

# Script and file patterns (under results/ and plot_data/) of every output.
EEG = {"fig5": [160, 248, 191, 183, 184, 220, 148, 177, 22, 157, 31, 151],
       "fig6": [181, 231, 200, 144, 159, 136, 223, 205, 25, 179, 219, 123],
       "fig7": [223, 22, 25], "fig8": [191, 184, 160]}
OUTPUTS = {
    "table1": ("table1_type1_error", ["results/table1_type1_error*"]),
    "table2": ("table2_power_fixed", ["results/table2_power_fixed*"]),
    "table3": ("table3_power_random", ["results/table3_power_random*"]),
    "table4": ("table4_eeg_pvalues", ["results/table4_eeg_pvalues*"]),
    "fig1": ("fig1_liver_cv", ["results/fig1_*", "plot_data/fig1*"]),
    "fig2": ("fig2_liver_heat", ["results/*_liver*", "plot_data/fig2*"]),
    "fig3": ("fig3_lambda_hist", ["results/fig3_*", "plot_data/fig3*"]),
    "fig5": ("fig5_eeg_heat", ["plot_data/fig5*"]),
    "fig6": ("fig6_eeg_heat", ["plot_data/fig6*"]),
    "fig7": ("fig7_power_gcv_decreasing", ["plot_data/fig7*"]),
    "fig8": ("fig8_power_gcv_increasing", ["plot_data/fig8*"]),
}
for figure, triples in EEG.items():
    runs = "runs3_" if figure in ("fig5", "fig6") else "runs1000_"
    for k in triples:
        OUTPUTS[figure][1].extend([f"results/power_eeg{k}_{runs}*", f"results/gcv_eeg{k}_*"])


def files_of(folder, patterns):
    """Relative paths of the csv files under a folder that match patterns."""
    found = set()
    for pattern in patterns:
        for path in glob.glob(os.path.join(folder, pattern)):
            if os.path.isdir(path):
                found.update(glob.glob(os.path.join(path, "*.csv")))
            elif path.endswith((".csv", ".csv.gz")):
                found.add(path)
    return sorted(os.path.relpath(path, folder) for path in found)


def read_table(path):
    """A csv file as header and rows of text cells."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", newline="") as f:
        rows = list(csv.reader(f))
    return rows[0], rows[1:]


def number(cell):
    """The number of a cell, NaN for a missing value, None for text."""
    if cell in ("NA", "", "NaN", "nan"):
        return math.nan
    try:
        return float(cell)
    except ValueError:
        return None


class Comparison:
    """Running totals of the comparison of one output."""

    def __init__(self):
        self.files = self.values = self.unequal = 0
        self.max_abs = self.max_rel = 0.0
        self.within = True
        self.problems = []

    def numbers(self, a, b, tolerance):
        self.values += 1
        if a == b or (math.isnan(a) and math.isnan(b)):
            return
        self.unequal += 1
        if math.isnan(a) or math.isnan(b) or math.isinf(a) or math.isinf(b):
            self.max_abs = self.max_rel = math.inf
            self.within = False
            return
        difference = abs(a - b)
        self.max_abs = max(self.max_abs, difference)
        self.max_rel = max(self.max_rel, difference / max(abs(a), abs(b)))
        if difference > tolerance * max(abs(a), abs(b), 1.0):
            self.within = False

    def problem(self, text):
        self.problems.append(text)
        self.within = False


def compare_columns(total, name, a, b, tolerance, file):
    """Compare one column of two tables, cell by cell."""
    for x, y in zip(a, b):
        nx, ny = number(x), number(y)
        if nx is None or ny is None:
            total.values += 1
            if x.lower() != y.lower():
                total.unequal += 1
                total.problem(f"{file}: text differs in column {name}")
        else:
            total.numbers(nx, ny, tolerance)


def compare_results(total, r_file, py_file, file):
    """Compare a csv file of results/: same header, same rows, every cell."""
    (r_head, r_rows), (py_head, py_rows) = read_table(r_file), read_table(py_file)
    if r_head != py_head or len(r_rows) != len(py_rows):
        total.problem(f"{file}: different columns or number of rows")
        return
    for j, name in enumerate(r_head):
        compare_columns(total, name, [row[j] for row in r_rows], [row[j] for row in py_rows],
                        TOLERANCE["results"], file)


def sort_key(row):
    """Order of the marks of a figure: by layer, then by position."""
    return tuple((1, 0.0) if math.isnan(value) else (0, value) for value in row)


def compare_plot_data(total, r_file, py_file, file):
    """Compare the plotted numbers of a figure, mark by mark.

    The marks of every layer are put in the same order (by position); the
    numeric columns must be equal; the colours of the two figures must
    correspond one to one.
    """
    NUMERIC = ["layer", "x", "xintercept", "y", "xmin", "xmax", "count", "value"]
    tables = []
    for path in (r_file, py_file):
        head, rows = read_table(path)
        columns = {name: [row[j] for row in rows] for j, name in enumerate(head)}
        tables.append(columns)
    r, py = tables
    if len(r["layer"]) != len(py["layer"]):
        total.problem(f"{file}: {len(r['layer'])} and {len(py['layer'])} marks")
        return
    # The group of a mark is compared where both figures give one.
    with_group = all(number(g) is not None and not math.isnan(number(g))
                     for g in r["group"] + py["group"]) and \
        sorted(set(r["group"])) == sorted(set(map(str, (int(float(g)) for g in py["group"]))))
    keys = ["layer"] + (["group"] if with_group else []) + NUMERIC[1:]

    def marks(table):
        rows = [[number(table[name][i]) for name in keys] + [table["colour"][i]]
                for i in range(len(table["layer"]))]
        return sorted(rows, key=lambda row: sort_key(row[:-1]))

    r_marks, py_marks = marks(r), marks(py)
    for j, name in enumerate(keys):
        for a, b in zip(r_marks, py_marks):
            total.numbers(a[j], b[j], TOLERANCE["plot_data"])
    pairs = {(a[-1], b[-1]) for a, b in zip(r_marks, py_marks) if a[-1] != "NA" and b[-1] != "NA"}
    if len({a for a, _ in pairs}) != len(pairs) or len({b for _, b in pairs}) != len(pairs):
        total.problem(f"{file}: the colours of the marks do not correspond")


def run_time(folder, script):
    """Seconds, peak memory (MB) and cores of a script, from its .done file."""
    try:
        with open(os.path.join(folder, "logs", script + ".done")) as f:
            fields = f.read().strip().split(",")
        return fields[0], fields[1], fields[2] if len(fields) > 2 else "4"
    except OSError:
        return "", "", ""


rows = []
for output, (script, patterns) in OUTPUTS.items():
    total = Comparison()
    r_files, py_files = files_of(r_dir, patterns), files_of(py_dir, patterns)
    if r_files != py_files:
        total.problem("the two implementations wrote different files")
    for file in sorted(set(r_files) & set(py_files)):
        compare = compare_plot_data if file.startswith("plot_data") else compare_results
        compare(total, os.path.join(r_dir, file), os.path.join(py_dir, file), file)
        total.files += 1
    if total.files == 0:
        total.problem("no files to compare")
    verdict = "DIFFERENT" if not total.within else "identical" if total.unequal == 0 else \
        "within tolerance"
    r_time, py_time = run_time(r_dir, script), run_time(py_dir, script)
    rows.append({"output": output, "mode": mode, "settings": SETTINGS.get(output, "archived runs" if mode == "archive" else "full"),
                 "files": total.files, "values": total.values, "not_identical": total.unequal,
                 "max_abs_diff": f"{total.max_abs:.3g}", "max_rel_diff": f"{total.max_rel:.3g}",
                 "tolerance": f"{max(TOLERANCE.values()):g}", "verdict": verdict,
                 "r_seconds": r_time[0], "python_seconds": py_time[0],
                 "r_peak_mb": r_time[1], "python_peak_mb": py_time[1],
                 "r_cores": r_time[2], "python_cores": py_time[2],
                 "problems": "; ".join(sorted(set(total.problems))[:5])})
    print(f"{output:7} {mode:8} files {total.files:3d} values {total.values:8d} "
          f"not identical {total.unequal:7d} max abs {total.max_abs:9.3g} "
          f"max rel {total.max_rel:9.3g}  {verdict}")
    for problem in sorted(set(total.problems))[:5]:
        print("        ", problem)

os.makedirs(os.path.dirname(os.path.abspath(out_file)), exist_ok=True)
with open(out_file, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
sys.exit(0 if all(row["verdict"] != "DIFFERENT" for row in rows) else 1)
