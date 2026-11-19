# -----------------------------------------------------------------------------
# File:    python/adamant/setup.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. All figures and tables.
# Purpose: Loaded first by every script: the settings of the paper, the
#          project paths, the command line parser, and the helpers the R
#          code takes from base R (csv files, parallel runs, messages).
# Inputs:  None.
# Outputs: None.
# Usage:   from adamant.setup import *   # in a script under python/scripts/
# Source:  R/setup.R.
# -----------------------------------------------------------------------------
import math
import multiprocessing
import os
import re
import sys
import time

import pandas as pd

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

# ----- Settings of the paper, passed to the functions as `settings` -----
PAPER = {
    "seed": 8891,        # seed of designs, permutations and replications
    "alpha": 0.05,       # significance level
    "sigma": 0.5,        # sd of the noise added to the response
    "upper": 1000,       # largest ridge penalty of every grid
    "n_perms": 1000,     # permutations per test
    "figure_size": 7,    # width and height of every figure, inches
}


# ----- Paths -----

def project_paths(root, data=None, output=None):
    """Folders and files of the project.

    Parameters: root (repository root), data (folder of the raw data; its
    subfolder eeg/ holds the EEG files; ADAMANT_DATA by default), output
    (folder in which results/ and figures/ are written; ADAMANT_OUTPUT by
    default).
    Returns: dict of paths: eeg, liver, reference, archive, results, figures.
    """
    data = data or os.environ.get("ADAMANT_DATA", os.path.join(root, "data"))
    output = output or os.environ.get("ADAMANT_OUTPUT", root)
    return {"eeg": os.path.join(data, "eeg"),
            "liver": os.path.join(root, "data", "liver", "mixOmicsData.csv"),
            "reference": os.path.join(root, "reference"),
            "archive": os.path.join(root, "archive"),
            "results": os.path.join(output, "results"),
            "figures": os.path.join(output, "figures")}


PATHS = project_paths(REPO_ROOT)


def output_file(folder, name):
    """Path of an output file, creating its folder when needed.

    Parameters: folder (output folder), name (file name).
    Returns: the path folder/name.
    """
    os.makedirs(folder, exist_ok=True)
    return os.path.join(folder, name)


# ----- Command line -----

def parse_args(defaults, args=None):
    """Read `--name=value` arguments and bare `--name` flags.

    Every script accepts `--rerun`: recompute from the data instead of
    reading the archived runs.

    Parameters: defaults (dict of options and their default values; a value
    given on the command line is converted to the type of its default, a
    flag sets a boolean default to True), args (command line arguments).
    Returns: `defaults` with the values given on the command line, and rerun.
    """
    args = sys.argv[1:] if args is None else args
    opts = {"rerun": False, **defaults}
    for arg in args:
        pair = re.fullmatch(r"--([a-z_]+)(=(.*))?", arg)
        if pair is None:
            raise SystemExit("arguments must look like --name=value or --name")
        name, is_flag = pair.group(1), pair.group(2) is None
        if name not in opts:
            raise SystemExit(f"unknown argument --{name}")
        if is_flag != isinstance(opts[name], bool):
            raise SystemExit(f"--{name} " + ("needs a value" if is_flag else "takes no value"))
        opts[name] = True if is_flag else type(opts[name])(pair.group(3))
    return opts


def curve_options(grid, settings):
    """Default options of a script that computes curves over a penalty grid.

    Parameters: grid (dict with the grid sizes negative and positive),
    settings (settings of the paper, n_perms).
    Returns: dict of defaults for `parse_args()`.
    """
    return {"runs": 1000, "perms": settings["n_perms"], "negative": grid["negative"],
            "positive": grid["positive"], "cores": 4}


def test_options(settings, cores=4):
    """Default options of a script that runs AdaMantCV tests.

    Parameters: settings (settings of the paper, n_perms), cores (default
    number of parallel processes).
    Returns: dict of defaults for `parse_args()`.
    """
    return {"perms": settings["n_perms"], "negative": 50, "positive": 200, "cores": cores}


def parse_index_set(text, count=False):
    """Turn an index text into indices.

    Parameters: text ("3:7" for 3 to 7, "2,5,9", or a single number N),
    count (if True a single number N means 1 to N).
    Returns: list of integers.
    """
    if re.fullmatch(r"[0-9]+:[0-9]+", text):
        first, last = (int(end) for end in text.split(":"))
        return list(range(first, last + 1))
    index = [int(item) for item in text.split(",")]
    return list(range(1, index[0] + 1)) if count and len(index) == 1 else index


def settings_tag(opts, n_runs=None):
    """Label of the settings a rerun was computed with.

    Used in the names of the stored parts of a rerun, so that runs with
    different settings are kept apart.

    Parameters: opts (options with perms, negative, positive), n_runs
    (number of replications, or None to leave it out).
    Returns: the label.
    """
    runs = "" if n_runs is None else f"runs{n_runs}_"
    return f"{runs}perms{opts['perms']}_grid{opts['negative']}x{opts['positive']}"


# ----- csv files, parallel runs, messages -----

def read_csv(file, **options):
    """`read.csv(file)`: a csv file as a data frame.

    Parameters: file (path; .gz is read directly), options (passed to
    pandas). Returns: pandas DataFrame.
    """
    return pd.read_csv(file, float_precision="round_trip", **options)


def format_number(x, digits=15):
    """A number as R writes it with `digits` significant digits.

    R finds the significant digits a number needs by scaling it to an
    integer of `digits` digits and dropping the zeros at its end, then
    writes it in fixed or scientific notation, whichever is shorter.

    Parameters: x (float), digits (significant digits, at most 15).
    Returns: the text.
    """
    if math.isnan(x):
        return "NA"
    if math.isinf(x):
        return "Inf" if x > 0 else "-Inf"
    if x == 0:
        return "0"
    negative, r = int(x < 0), abs(x)
    kp = math.floor(math.log10(r)) - digits + 1
    if abs(kp) <= 22:  # powers of ten that are exact
        scaled = r / 10.0 ** kp if kp >= 0 else r * 10.0 ** -kp
    elif kp <= -308:
        scaled = (r * 1e+303) / math.pow(10.0, kp + 303)
    else:
        scaled = r / math.pow(10.0, kp)
    if scaled < 10.0 ** (digits - 1):
        scaled *= 10.0
        kp -= 1
    alpha = float(round(scaled))  # to the nearest integer, ties to even
    significant = digits
    for _ in range(digits):
        alpha /= 10.0
        if alpha != math.floor(alpha):
            break
        significant -= 1
    if significant == 0:
        significant = 1
        kp += 1
    power = kp + digits - 1
    left = power + 1
    cut = min(max(digits - power, 0), 22)  # the decimal place rounding cuts off
    if 0 < power <= 22 and r < 10.0 ** power - 0.5 / 10.0 ** cut:
        left -= 1
    decimals = max(0, significant - power - 1)
    fixed_width = negative + max(left, 1) + decimals + (decimals > 0)
    mantissa = significant - 1
    scientific_width = negative + (mantissa > 0) + mantissa + 4 + (2 if abs(power) >= 100 else 1)
    return f"{x:.{decimals}f}" if fixed_width <= scientific_width else f"{x:.{mantissa}e}"


def write_csv(table, file):
    """`write.csv(table, file, row.names = FALSE)`: numbers with 15
    significant digits as R writes them, text in quotes.

    Parameters: table (DataFrame), file (path). Returns: None.
    """
    def quoted(text):
        return '"' + str(text).replace('"', '""') + '"'

    columns = []
    for name in table.columns:
        values = table[name]
        if pd.api.types.is_bool_dtype(values):
            columns.append(["TRUE" if value else "FALSE" for value in values])
        elif pd.api.types.is_integer_dtype(values):
            columns.append([str(value) for value in values])
        elif pd.api.types.is_float_dtype(values):
            columns.append([format_number(float(value)) for value in values])
        else:
            columns.append([quoted(value) for value in values])
    with open(file, "w") as f:
        f.write(",".join(quoted(name) for name in table.columns) + "\n")
        for row in zip(*columns):
            f.write(",".join(row) + "\n")


def cached_csv(file, compute):
    """Read a result from a csv file if it exists; otherwise compute and
    store it.

    Parameters: file (path of the csv file), compute (function without
    arguments that returns a DataFrame).
    Returns: the DataFrame.
    """
    if os.path.exists(file):
        log_message("reusing %s", file)
        return read_csv(file)
    result = compute()
    write_csv(result, file)
    return result


_job = None


def _run_job(item):
    return _job(item)


def parallel_map(function, items, cores=1):
    """`parallel::mclapply(items, function, mc.cores = cores)`.

    The processes are forked, so each starts with the memory and the random
    number stream of this one.

    Parameters: function (of one item), items (list), cores (number of
    parallel processes).
    Returns: list of the results, in the order of `items`.
    """
    global _job
    items = list(items)
    cores = int(cores)
    if cores <= 1 or len(items) <= 1:
        return [function(item) for item in items]
    _job = function
    with multiprocessing.get_context("fork").Pool(min(cores, len(items))) as pool:
        return pool.map(_run_job, items, chunksize=1)


def log_message(msg, *values):
    """Print a time-stamped message.

    Parameters: msg (format of the message, with % fields), values (values
    for the format). Returns: None.
    """
    print(time.strftime("%Y-%m-%d %H:%M:%S"), "-", msg % values, flush=True)
