# -----------------------------------------------------------------------------
# File:    python/tests/test_rng.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Check that python/adamant/rng.py draws the same numbers as R, bit
#          for bit, for several seeds and sizes and for the call sequences of
#          the R code: the permutations of adamant_cv() and power_run(), the
#          folds of cv.glmnet(), the gene sample of Figure 1, the noise and
#          the coefficients of the simulations, and a simulated design
#          (MASS::mvrnorm and scale).
#          Also checks the first draws of a seed against the numbers R
#          prints, without running R.
# Inputs:  None (R is run through Rscript; numbers are passed in binary).
# Outputs: None (pytest reports the sequences that differ).
# Usage:   python -m pytest python/tests/test_rng.py     # inside `nix develop`
# -----------------------------------------------------------------------------
import os
import subprocess
import sys
import tempfile

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from adamant import rbase, rng  # noqa: E402
from adamant.simulation import make_covariance  # noqa: E402

SEEDS = [8891, 1, 42, 889192, 889193, 18891, 2147483647]

R_CODE = r"""
out <- file(commandArgs(TRUE)[1], "wb")
put <- function(x) writeBin(as.double(x), out)
for (seed in c(%s)) {
  set.seed(seed)
  put(sample.int(10)); put(runif(3)); put(rnorm(3))
  put(sample.int(350)); put(rnorm(1000, 1, 0.5)); put(sample.int(65536, 20)); put(sample.int(65537, 20))
  put(sample.int(1)); put(sample.int(2)); put(rnorm(350, mean = 0, sd = 0.5))
  # adamant_cv(): permuted responses
  Y <- matrix(rnorm(64), 64, 1)
  put(vapply(seq_len(25), function(b) Y[sample.int(64), ], numeric(64)))
  # cv.glmnet(): folds; Figure 1: gene sample, then folds
  put(sample(rep(seq_len(10), length = 64)))
  put(sample(seq_len(3116), 50, replace = FALSE))
  for (i in 1:11) sample(rep(seq_len(10), length = 64))
  put(sample(rep(seq_len(10), length = 64)))
  # power_run(): noise, then row permutations
  put(rnorm(350, mean = 0, sd = 0.5)); for (b in 1:5) put(sample.int(350))
}
# simulate_design(): mvrnorm, scale, coefficients
for (type in c("CSYM", "NCDIAG")) {
  set.seed(8891)
  Sigma <- if (type == "CSYM") { S <- matrix(0.025, 40, 40); diag(S) <- 1; S } else diag(log(2:41), 40, 40)
  put(scale(MASS::mvrnorm(30, rep(0, 40), Sigma = Sigma))); put(rnorm(40, 1, 0.5))
}
close(out)
"""


def python_draws():
    """The draws of R_CODE from rng.py, as (label, values) pairs."""
    draws = []

    def put(label, values):
        draws.append((label, np.asarray(values, dtype=np.float64).ravel(order="F")))

    for seed in SEEDS:
        rng.set_seed(seed)
        put(f"seed {seed}: sample.int(10)", rng.sample_int(10))
        put(f"seed {seed}: runif(3)", rng.runif(3))
        put(f"seed {seed}: rnorm(3)", rng.rnorm(3))
        put(f"seed {seed}: sample.int(350)", rng.sample_int(350))
        put(f"seed {seed}: rnorm(1000, 1, 0.5)", rng.rnorm(1000, 1, 0.5))
        put(f"seed {seed}: sample.int(65536, 20)", rng.sample_int(65536, 20))
        put(f"seed {seed}: sample.int(65537, 20)", rng.sample_int(65537, 20))
        put(f"seed {seed}: sample.int(1)", rng.sample_int(1))
        put(f"seed {seed}: sample.int(2)", rng.sample_int(2))
        put(f"seed {seed}: rnorm(350, 0, 0.5)", rng.rnorm(350, 0, 0.5))
        Y = rng.rnorm(64)
        put(f"seed {seed}: 25 permuted responses", Y[rng.permutation_rows(64, 25)].T)
        folds = np.resize(np.arange(1, 11), 64)
        put(f"seed {seed}: folds", rng.sample(folds))
        put(f"seed {seed}: sample of 50 genes", rng.sample_int(3116, 50))
        for _ in range(11):
            rng.sample(folds)
        put(f"seed {seed}: folds after 11 skipped", rng.sample(folds))
        put(f"seed {seed}: noise", rng.rnorm(350, 0, 0.5))
        for b in range(5):
            put(f"seed {seed}: row permutation {b + 1}", rng.sample_int(350))
    for type in ("CSYM", "NCDIAG"):
        rng.set_seed(8891)
        put(f"design {type}: X", rbase.scale(rbase.mvrnorm(30, make_covariance(type, 40, 0.025))))
        put(f"design {type}: coefficients", rng.rnorm(40, 1, 0.5))
    return draws


@pytest.fixture(scope="module")
def r_draws():
    """The numbers R draws with R_CODE, in order."""
    with tempfile.TemporaryDirectory() as folder:
        file = os.path.join(folder, "draws.bin")
        subprocess.run(["Rscript", "-e", R_CODE % ", ".join(map(str, SEEDS)), file], check=True)
        return np.fromfile(file)


@pytest.mark.against_r
def test_every_sequence_is_identical_to_r(r_draws):
    draws = python_draws()
    assert sum(len(values) for _, values in draws) == len(r_draws), "different number of draws"
    position, different = 0, []
    for label, values in draws:
        expected = r_draws[position:position + len(values)]
        position += len(values)
        if not np.array_equal(values, expected):
            different.append(f"{label}: {int((values != expected).sum())} of {len(values)}")
    print(f"rng: {len(draws) - len(different)} of {len(draws)} sequences identical to R "
          f"({len(r_draws)} numbers, {len(SEEDS)} seeds)")
    assert not different, different


def test_the_first_draws_of_seed_1_are_those_r_prints():
    # set.seed(1); runif(2) / rnorm(2) / sample.int(5); sample.int(10, 3) in R 4.6.1
    rng.set_seed(1)
    assert rng.runif(2).tolist() == [0.26550866314209998, 0.37212389963679016]
    rng.set_seed(1)
    assert rng.rnorm(2).tolist() == [-0.62645381074233242, 0.18364332422208227]
    rng.set_seed(1)
    assert rng.sample_int(5).tolist() == [1, 4, 3, 5, 2]
    assert rng.sample_int(10, 3).tolist() == [7, 2, 3]
    rng.set_seed(42)
    assert rng.sample(np.array([10, 20, 30, 40])).tolist() == [10, 40, 30, 20]


def test_set_seed_restarts_the_stream():
    rng.set_seed(8891)
    first = (rng.runif(3), rng.rnorm(4, 1, 0.5), rng.sample_int(6))
    rng.set_seed(8891)
    again = (rng.runif(3), rng.rnorm(4, 1, 0.5), rng.sample_int(6))
    assert all(np.array_equal(a, b) for a, b in zip(first, again))
    assert not np.array_equal(rng.runif(3), first[0])


def test_the_draws_have_the_range_and_the_moments_they_should():
    rng.set_seed(2)
    uniform = rng.runif(20000)
    assert ((uniform > 0) & (uniform < 1)).all() and abs(uniform.mean() - 0.5) < 0.01
    normal = rng.rnorm(20000, 3, 2)
    assert abs(normal.mean() - 3) < 0.05 and abs(normal.std() - 2) < 0.05
    assert sorted(rng.sample_int(50)) == list(range(1, 51))
    part = rng.sample_int(50, 7)
    assert len(set(part)) == 7 and part.min() >= 1 and part.max() <= 50


def test_permutation_rows_is_consecutive_calls_of_sample_int():
    rng.set_seed(3)
    rows = rng.permutation_rows(6, 4)
    assert rows.shape == (4, 6)
    rng.set_seed(3)
    for row in rows:
        assert np.array_equal(row, rng.sample_int(6) - 1)


def test_muladd_is_a_times_b_plus_c():
    from numba import njit

    @njit
    def combine(a, b, c):
        return rng.muladd(a, b, c)

    assert combine(2.0, 3.0, 1.0) == 7.0
    assert isinstance(rng.FUSED, bool)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q", "-s"]))
