# -----------------------------------------------------------------------------
# File:    python/tests/test_rmath.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Methods (all figures and tables).
# Purpose: Check that adamant.rmath.pgamma_upper() returns the same doubles,
#          bit for bit, as R's pgamma(q, shape, rate = rate,
#          lower.tail = FALSE), over every branch of the algorithm and over
#          the degenerate moment-matched fits (huge shape and rate).
# Inputs:  None (Rscript must be on the PATH; the cases are generated here
#          and go to R and back as binary doubles).
#          Also checks values with a closed form, without running R.
# Outputs: None (with `-s` pytest shows the table of cases, mismatches and
#          branch counts).
# Usage:   nix develop -c python -m pytest python/tests/test_rmath.py
# -----------------------------------------------------------------------------
import os
import subprocess
import sys
import tempfile
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import numpy as np
import pytest

from adamant.rmath import BRANCHES, DBL_EPSILON, DBL_MIN, pgamma_branch, pgamma_upper

R_SCRIPT = """
a <- commandArgs(trailingOnly = TRUE)
con <- file(a[1], "rb")
n <- readBin(con, "integer", 1L)
x <- readBin(con, "double", n)
shape <- readBin(con, "double", n)
rate <- readBin(con, "double", n)
close(con)
p <- suppressWarnings(pgamma(x, shape, rate = rate, lower.tail = FALSE))
stopifnot(length(p) == n)
writeBin(p, a[2])
"""


def r_pgamma_upper(x, shape, rate):
    """pgamma(x, shape, rate = rate, lower.tail = FALSE) computed by R.

    Parameters: x, shape, rate (1-d float64 arrays of one length).
    Returns: float64 array, the doubles R returned.
    """
    with tempfile.TemporaryDirectory() as tmp:
        fin = os.path.join(tmp, "in.bin")
        fout = os.path.join(tmp, "out.bin")
        script = os.path.join(tmp, "pgamma.R")
        with open(script, "w") as handle:
            handle.write(R_SCRIPT)
        with open(fin, "wb") as handle:
            np.array([x.size], dtype=np.int32).tofile(handle)
            for a in (x, shape, rate):
                np.ascontiguousarray(a, dtype=np.float64).tofile(handle)
        subprocess.run(["Rscript", "--vanilla", script, fin, fout], check=True)
        return np.fromfile(fout, dtype=np.float64)


def log_uniform(rng, low, high, n):
    """n numbers whose logarithm is uniform between log(low) and log(high).

    Parameters: rng (numpy Generator), low, high (floats), n (int).
    Returns: float64 array.
    """
    return np.exp(rng.uniform(np.log(low), np.log(high), n))


def near_one(rng, n, low=0.2, high=5.0):
    """n ratios in [low, high], half of them concentrated near 1.

    Parameters: rng (numpy Generator), n (int), low, high (floats).
    Returns: float64 array.
    """
    wide = log_uniform(rng, low, high, n)
    close = np.exp(rng.normal(0.0, 1.0, n) * log_uniform(rng, 1e-6, 0.3, n))
    ratio = np.where(rng.random(n) < 0.5, wide, close)
    return np.clip(ratio, low, high)


def case_branches(rng, n=200_000):
    """Shape in [0.01, 1e4], x * rate / shape in [0.2, 5] around 1."""
    shape = log_uniform(rng, 0.01, 1e4, n)
    rate = log_uniform(rng, 1e-3, 1e3, n)
    x = near_one(rng, n) * shape / rate
    return x, shape, rate


def case_degenerate(rng, n=20_000):
    """Moment-matched fits whose variance is rounding noise."""
    m = 349.0 * (1.0 + rng.uniform(-0.05, 0.05, n))
    v = log_uniform(rng, 1e-27, 1e-20, n)
    k = rng.integers(-24, 25, n)
    shape = m * m / v
    rate = m / v
    x = m * (1.0 + k * 2.0**-52)
    return x, shape, rate


def case_small_x(rng, n=50_000):
    """x * rate < 1 with shapes from tiny to huge."""
    shape = log_uniform(rng, 1e-300, 1e300, n)
    shape[: n // 2] = log_uniform(rng, 1e-6, 1e3, n // 2)
    rate = log_uniform(rng, 1e-3, 1e3, n)
    x = log_uniform(rng, 1e-300, 1.0, n) / rate
    x[: n // 4] = rng.uniform(0.0, 1.0, n // 4) / rate[: n // 4]
    return x, shape, rate


def case_wide(rng, n=100_000):
    """Tiny and huge shapes, x from far below to far above the mean."""
    shape = log_uniform(rng, 1e-12, 1e30, n)
    rate = log_uniform(rng, 1e-6, 1e6, n)
    x = near_one(rng, n, 1e-3, 1e3) * shape / rate
    part = n // 5
    # x * rate >= 1 with shape below 1, and integer and half-integer shapes
    x[:part] = log_uniform(rng, 1.0, 1e3, part) / rate[:part]
    shape[:part] = log_uniform(rng, 1e-17, 1.0, part)
    shape[part : 2 * part] = rng.integers(1, 80, part) / 2.0
    x[part : 2 * part] = near_one(rng, part, 0.05, 20.0) * shape[part : 2 * part] / rate[part : 2 * part]
    return x, shape, rate


def case_special(rng):
    """Zeros, infinities, NaN and invalid parameters, in every combination."""
    values = np.array(
        [0.0, -0.0, -1.0, 1.0, 0.5, 2.5, 1e-320, 5e-324, DBL_MIN, 1e-300, 1e300, 1.7976931348623157e308, np.inf, -np.inf, np.nan, 349.0, 1e25]
    )
    grid = np.array(np.meshgrid(values, values, values, indexing="ij")).reshape(3, -1)
    return grid[0], grid[1], grid[2]


CASES = (
    ("1 branches of pgamma_raw", case_branches),
    ("2 degenerate fits", case_degenerate),
    ("3a small x", case_small_x),
    ("3b tiny and huge shapes", case_wide),
    ("3c zero, Inf, NaN, invalid", case_special),
)


@pytest.fixture(scope="module")
def cases():
    """The cases of every class, drawn in the order of CASES from one seed.

    Returns: dict from the name of a class to its (x, shape, rate).
    """
    rng = np.random.default_rng(20261119)
    return {name: make(rng) for name, make in CASES}


@pytest.mark.against_r
@pytest.mark.parametrize("name", [name for name, _ in CASES])
def test_every_double_equals_r(cases, name):
    x, shape, rate = cases[name]
    pgamma_upper(1.0, 1.0, 1.0)  # compile
    start = time.perf_counter()
    ours = pgamma_upper(x, shape, rate)
    seconds = time.perf_counter() - start
    theirs = r_pgamma_upper(x, shape, rate)
    both_nan = np.isnan(ours) & np.isnan(theirs)
    differ = (ours.view(np.uint64) != theirs.view(np.uint64)) & ~both_nan
    print(f"{name}: {x.size} cases, {int(differ.sum())} differ from R, {x.size / seconds:,.0f} calls per second")
    branch = pgamma_branch(x, shape, rate)
    redone = (ours < DBL_MIN / DBL_EPSILON) & (branch > 0)
    for index, label in enumerate(BRANCHES):
        chosen = branch == index
        if chosen.any():
            print(f"    {int(chosen.sum()):7d} ({int((differ & chosen).sum())} differ)  {label}")
    print(f"    {int(redone.sum()):7d} ({int((differ & redone).sum())} differ)  of these, redone through the logarithm")
    examples = [f"x={x[i].hex()} shape={shape[i].hex()} rate={rate[i].hex()} ours={ours[i]!r} R={theirs[i]!r}"
                for i in np.flatnonzero(differ)[:5]]
    assert not differ.any(), f"{int(differ.sum())} doubles differ from R: {examples}"


def test_the_cases_reach_every_branch(cases):
    reached = set()
    for x, shape, rate in cases.values():
        reached |= set(np.unique(pgamma_branch(x, shape, rate)).tolist())
    assert reached == set(range(len(BRANCHES)))


def test_pgamma_upper_has_the_closed_forms_of_small_integer_shapes():
    # shape 1 is the exponential; integer shape k: P(Poisson(x * rate) <= k - 1)
    x = np.array([0.0, 0.5, 1.0, 2.0, 30.0])
    np.testing.assert_allclose(pgamma_upper(x, 1.0, 1.0), np.exp(-x), rtol=1e-14)
    np.testing.assert_allclose(pgamma_upper(x, 1.0, 3.0), np.exp(-3 * x), rtol=1e-14)
    np.testing.assert_allclose(pgamma_upper(3.0, 4.0, 2.0), 61 * np.exp(-6.0), rtol=1e-14)
    np.testing.assert_allclose(pgamma_upper(2.0, 2.0, 1.0), 3 * np.exp(-2.0), rtol=1e-14)


def test_pgamma_upper_broadcasts_and_handles_the_edges():
    grid = pgamma_upper(np.array([[1.0], [2.0]]), np.array([1.0, 2.0, 3.0]), 1.0)
    assert grid.shape == (2, 3)
    assert (np.diff(grid, axis=1) > 0).all() and (np.diff(grid, axis=0) < 0).all()
    assert pgamma_upper(0.0, 2.0, 1.0) == 1 and pgamma_upper(np.inf, 2.0, 1.0) == 0
    assert pgamma_upper(-1.0, 2.0, 1.0) == 1
    assert np.isnan(pgamma_upper(np.nan, 2.0, 1.0)) and np.isnan(pgamma_upper(1.0, -1.0, 1.0))


def test_pgamma_branch_names_the_branch_of_a_call():
    branch = pgamma_branch(np.array([-1.0, 0.5, 2.0, 2.0, 50.0]),
                           np.array([2.0, 2.0, 10.0, 0.5, 10.0]), 1.0)
    assert [BRANCHES[i].split(" ")[0] for i in branch] == \
        ["outside", "pgamma_smallx", "pd_upper_series", "pd_lower_cf", "pd_lower_series"]


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q", "-s"]))
