# -----------------------------------------------------------------------------
# File:    python/tests/test_glmnet.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Figure 1.
# Purpose: Compares python/adamant/glmnet.py with glmnet::cv.glmnet() run
#          through Rscript on the same data and folds: cross-validated MSE,
#          lambda.min, the coefficients of the full fit and the number of
#          passes of the full fit and of every fold fit.
# Inputs:  data/liver/mixOmicsData.csv.
# Outputs: None (prints the differences; exits non-zero when a check fails).
#          Also checks the ridge path against its closed form, without
#          running R.
# Usage:   nix develop -c python -m pytest python/tests/test_glmnet.py          (no full liver grids)
#          nix develop -c python -m pytest python/tests/test_glmnet.py --full   (everything)
#          nix develop -c python python/tests/test_glmnet.py [quick]            (everything, or as the first)
# Source:  glmnet 5.0 (cv.glmnet, elnet; gaussian family, alpha = 0).
# -----------------------------------------------------------------------------

import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from adamant import glmnet  # noqa: E402

LIVER = ROOT / "data" / "liver" / "mixOmicsData.csv"
CVM_TOLERANCE = 1e-8  # largest relative difference of cvm allowed
COEF_TOLERANCE = 1e-8  # largest difference of a coefficient, relative to the largest one

# The data go to R as raw doubles, so both sides start from the same bits.
R_SCRIPT = r"""
args <- commandArgs(TRUE)
dir <- args[1]
fold_passes <- args[2] == "1"
rd <- function(name) {
  f <- file.path(dir, name)
  readBin(f, "double", n = file.size(f) / 8)
}
dims <- rd("dims.bin"); n <- dims[1]; p <- dims[2]; q <- dims[3]
X <- matrix(rd("X.bin"), n, p)
Y <- matrix(rd("Y.bin"), n, q)
lam <- rd("lam.bin")
folds <- matrix(as.integer(rd("fold.bin")), n, q)
suppressPackageStartupMessages(library(glmnet))
out <- file(file.path(dir, "out.bin"), "wb")
elapsed <- 0
for (j in seq_len(q)) {
  f <- folds[, j]
  elapsed <- elapsed + system.time(
    fit <- suppressWarnings(cv.glmnet(X, Y[, j], alpha = 0, nfolds = max(f),
                                      lambda = lam, foldid = f)))[["elapsed"]]
  passes <- rep(-1, max(f))
  if (fold_passes) passes <- sapply(seq_len(max(f)), function(i)
    suppressWarnings(glmnet(X[f != i, , drop = FALSE], Y[f != i, j], alpha = 0,
                            lambda = lam))$npasses)
  writeBin(as.double(c(length(fit$cvm), fit$glmnet.fit$npasses, fit$lambda.min,
                       length(passes), passes, fit$cvm, fit$lambda,
                       as.matrix(coef(fit$glmnet.fit)))), out)
}
writeBin(as.double(elapsed), out)
close(out)
"""


def run_r(X, Y, lambdas, folds, fold_passes=True):
    """cv.glmnet(X, Y[, j], alpha = 0, lambda = lambdas, foldid = folds[, j]) per column j.

    Returns: list of dicts (cvm, lambda, lambda_min, npasses, fold_npasses,
    a0, beta) and the seconds R spent in cv.glmnet.
    """
    n, p = X.shape
    q = Y.shape[1]
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        np.array([n, p, q], dtype="<f8").tofile(tmp / "dims.bin")
        np.asfortranarray(X, dtype="<f8").T.tofile(tmp / "X.bin")
        np.asfortranarray(Y, dtype="<f8").T.tofile(tmp / "Y.bin")
        np.asarray(lambdas, dtype="<f8").tofile(tmp / "lam.bin")
        np.asfortranarray(folds, dtype="<f8").T.tofile(tmp / "fold.bin")
        (tmp / "run.R").write_text(R_SCRIPT)
        subprocess.run(["Rscript", "--vanilla", str(tmp / "run.R"), str(tmp),
                        "1" if fold_passes else "0"], check=True)
        raw = np.fromfile(tmp / "out.bin", dtype="<f8")
    results, pos = [], 0
    for _ in range(q):
        nl, npasses, lambda_min, nf = int(raw[pos]), int(raw[pos + 1]), raw[pos + 2], int(raw[pos + 3])
        pos += 4
        fold_npasses = raw[pos:pos + nf].astype(int); pos += nf
        cvm = raw[pos:pos + nl]; pos += nl
        lam = raw[pos:pos + nl]; pos += nl
        coef = raw[pos:pos + (p + 1) * nl].reshape((p + 1, nl), order="F"); pos += (p + 1) * nl
        results.append(dict(cvm=cvm, lambda_min=lambda_min, npasses=npasses,
                            fold_npasses=fold_npasses, a0=coef[0], beta=coef[1:],
                            **{"lambda": lam}))
    return results, raw[pos]


def compare(name, X, Y, lambdas, folds, fold_passes=True):
    """Run both implementations on every column of Y, print and check the differences.

    Returns: dict of the largest differences and the two timings.
    """
    r_results, r_seconds = run_r(X, Y, lambdas, folds, fold_passes)
    worst = dict(cvm=0.0, coef=0.0, lam=0.0)
    lambda_min_equal = passes_equal = fold_passes_equal = bitwise_cvm = True
    start = time.perf_counter()
    py_results = [(glmnet.cv_glmnet_ridge(X, Y[:, j], lambdas, folds[:, j]),
                   glmnet.glmnet_ridge(X, Y[:, j], lambdas)) for j in range(Y.shape[1])]
    py_seconds = time.perf_counter() - start
    for (cv, (a0, beta)), r in zip(py_results, r_results):
        assert cv["cvm"].shape == r["cvm"].shape, name
        worst["cvm"] = max(worst["cvm"], np.max(np.abs(cv["cvm"] - r["cvm"]) / np.abs(r["cvm"])))
        scale = max(np.max(np.abs(r["beta"])), np.max(np.abs(r["a0"])))
        worst["coef"] = max(worst["coef"], np.max(np.abs(beta - r["beta"])) / scale,
                            np.max(np.abs(a0 - r["a0"])) / scale)
        worst["lam"] = max(worst["lam"], np.max(np.abs(cv["lambda"] - r["lambda"])))
        lambda_min_equal &= cv["lambda_min"] == r["lambda_min"]
        passes_equal &= cv["npasses"] == r["npasses"]
        bitwise_cvm &= bool(np.array_equal(cv["cvm"], r["cvm"]))
        if fold_passes:
            fold_passes_equal &= bool(np.array_equal(cv["fold_npasses"], r["fold_npasses"]))
    print(f"{name}: n={X.shape[0]} p={X.shape[1]} nlam={len(lambdas)} responses={Y.shape[1]}\n"
          f"    cvm max rel diff {worst['cvm']:.3g} (bitwise equal: {bitwise_cvm}); "
          f"coef max diff / max|coef| {worst['coef']:.3g}; lambda max diff {worst['lam']:.3g}\n"
          f"    lambda.min identical: {lambda_min_equal}; passes of full fit equal: {passes_equal}; "
          f"passes of fold fits equal: {fold_passes_equal if fold_passes else 'not compared'}\n"
          f"    seconds: R cv.glmnet {r_seconds:.2f}, Python cv + full fit {py_seconds:.2f}",
          flush=True)
    assert worst["cvm"] < CVM_TOLERANCE, name
    assert worst["coef"] < COEF_TOLERANCE, name
    assert worst["lam"] == 0.0, name
    assert lambda_min_equal and passes_equal and fold_passes_equal, name
    return dict(worst, r_seconds=r_seconds, py_seconds=py_seconds)


def make_folds(n, nfolds, q, seed):
    """One random fold assignment (rep(1..nfolds, length n), permuted) per response."""
    rng = np.random.default_rng(seed)
    base = np.arange(n) % nfolds + 1
    return np.column_stack([rng.permutation(base) for _ in range(q)])


def random_problem(n, p, q, seed, signal=5):
    """Correlated gaussian predictors and responses with a few non-zero coefficients."""
    rng = np.random.default_rng(seed)
    X = rng.standard_normal((n, p)) + 0.7 * rng.standard_normal((n, 1))
    X *= rng.uniform(0.5, 3.0, size=p)
    B = np.zeros((p, q))
    B[:signal] = rng.standard_normal((signal, q))
    Y = 1.5 + X @ B + rng.standard_normal((n, q))
    return X, Y


def read_liver():
    """Liver toxicity data, columns centred and scaled as R's scale() does."""
    data = np.genfromtxt(LIVER, delimiter=",", skip_header=1,
                         usecols=range(1, len(LIVER.open().readline().split(","))))
    scaled = (data - data.mean(axis=0)) / data.std(axis=0, ddof=1)
    return scaled[:, 10:], scaled[:, :10]


@pytest.mark.against_r
def test_random_covariance():
    X, Y = random_problem(40, 12, 2, seed=1)
    compare("random, covariance, n > p", X, Y, 10.0 ** np.linspace(2, -3, 60), make_folds(40, 5, 2, 11))
    X, Y = random_problem(25, 80, 2, seed=2)
    compare("random, covariance, n < p", X, Y, 10.0 ** np.linspace(3, -2, 80), make_folds(25, 5, 2, 12))


@pytest.mark.against_r
def test_random_naive():
    X, Y = random_problem(30, 600, 2, seed=3)
    compare("random, naive, p = 600", X, Y, 10.0 ** np.linspace(4, -2, 40), make_folds(30, 6, 2, 13))


@pytest.mark.against_r
def test_random_special():
    X, Y = random_problem(21, 9, 1, seed=4)
    X[:, 3] = 2.0  # a constant column, which glmnet leaves out
    compare("random, constant column, unsorted lambdas, fewer than 3 per fold",
            X, Y, np.array([0.5, 30.0, 0.01, 2.0, 7.0, 0.1]), make_folds(21, 10, 1, 14))


@pytest.mark.against_r
def test_liver_subset_thinned():
    X, Y = read_liver()
    genes = np.arange(0, 3100, 62)  # 50 fixed gene columns
    lambdas = 10.0 ** np.arange(4, -2.5 - 1e-9, -0.001)[::10]
    compare("liver A, 50 genes, every 10th lambda", X[:, genes], Y[:, :3], lambdas,
            make_folds(X.shape[0], 10, 3, 21))


@pytest.mark.against_r
@pytest.mark.full
def test_liver_subset_full_grid():
    X, Y = read_liver()
    genes = np.arange(0, 3100, 62)
    lambdas = 10.0 ** np.arange(4, -2.5 - 1e-9, -0.001)
    assert len(lambdas) == 6501
    compare("liver A, 50 genes, 6501 lambdas", X[:, genes], Y, lambdas,
            make_folds(X.shape[0], 10, 10, 22), fold_passes=False)


@pytest.mark.against_r
@pytest.mark.full
def test_liver_all_genes():
    X, Y = read_liver()
    lambdas = 10.0 ** np.arange(9, -5 - 1e-9, -0.1)
    assert len(lambdas) == 141
    compare("liver B, 3116 genes, 141 lambdas", X, Y, lambdas,
            make_folds(X.shape[0], 10, 10, 23), fold_passes=False)


def test_glmnet_ridge_is_the_closed_form_of_ridge_regression():
    # glmnet minimizes sum((y - a - X b)^2) / (2 n) + lambda * |b|^2 / 2 on
    # predictors and response scaled to variance one (divisor n): on centred
    # data b = (Z'Z + n lambda sd(y) I)^-1 Z'y per sd of each predictor
    X, Y = random_problem(40, 5, 1, seed=5)
    y = Y[:, 0]
    lambdas = np.array([0.01, 5.0, 0.5])              # used in decreasing order
    a0, beta = glmnet.glmnet_ridge(X, y, lambdas)
    assert a0.shape == (3,) and beta.shape == (5, 3)
    centre, spread = X.mean(axis=0), X.std(axis=0)
    Z = (X - centre) / spread
    for column, lam in enumerate([5.0, 0.5, 0.01]):
        b = np.linalg.solve(Z.T @ Z + 40 * lam / y.std() * np.eye(5), Z.T @ (y - y.mean()))
        # coordinate descent stops at thresh = 1e-7, short of the exact solution
        np.testing.assert_allclose(beta[:, column], b / spread, rtol=5e-3, atol=5e-4)
        np.testing.assert_allclose(a0[column], y.mean() - centre @ (b / spread), rtol=5e-3)


def test_cv_glmnet_ridge_returns_the_curve_and_its_minimum():
    X, Y = random_problem(40, 5, 1, seed=5)
    lambdas = np.array([0.01, 5.0, 0.5])
    folds = make_folds(40, 5, 1, 15)[:, 0]
    cv = glmnet.cv_glmnet_ridge(X, Y[:, 0], lambdas, folds)
    assert cv["lambda"].tolist() == [5.0, 0.5, 0.01]
    assert cv["cvm"].shape == (3,) and (cv["cvm"] > 0).all()
    assert cv["cvm"][0] > cv["cvm"][2]                # the largest penalty fits worst
    assert cv["lambda_min"] == cv["lambda"][np.argmin(cv["cvm"])]
    assert cv["index_min"] == int(np.argmin(cv["cvm"]))
    assert len(cv["fold_npasses"]) == 5
    # with two folds the curve is the mean squared error of each half
    # predicted from the fit on the other, weighted by the sizes of the halves
    halves = np.where(np.arange(40) < 15, 1, 2)
    two = glmnet.cv_glmnet_ridge(X, Y[:, 0], lambdas, halves)
    squared = np.empty((40, 3))
    for fold in (1, 2):
        a0, beta = glmnet.glmnet_ridge(X[halves != fold], Y[halves != fold, 0], lambdas)
        squared[halves == fold] = (Y[halves == fold, 0][:, None]
                                   - (a0 + X[halves == fold] @ beta)) ** 2
    np.testing.assert_allclose(two["cvm"], squared.mean(axis=0), rtol=1e-10)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q", "-s"] + ([] if "quick" in sys.argv[1:] else ["--full"])))
