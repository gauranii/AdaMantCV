# -----------------------------------------------------------------------------
# File:    python/adamant/rbase.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. All figures and tables.
# Purpose: What base R, stats and MASS supply to the code in R/, with R's
#          order of operations: matrix products, svd(), eigen() and ginv()
#          through the BLAS and LAPACK library R uses, sums and means that
#          add in R's order, scale(), quantile(), seq(), var(), pgamma() and
#          mvrnorm().
# Inputs:  The environment variable ADAMANT_BLAS names the library (OpenBLAS,
#          set by flake.nix to the one R is linked to). Without it the
#          products and decompositions come from numpy. The library may be
#          built with 32-bit or with 64-bit integers; the calls pass the
#          width it reports.
#          R adds its sums and means in the C type long double, which is a
#          double on arm64 macOS, 80 bits wide on x86_64 and 128 bits wide on
#          arm64 Linux; the sums here are added in numpy's longdouble, the
#          same type, so they round as R's do on each of them.
# Outputs: None.
# Usage:   from adamant.rbase import matprod, crossprod, svd, scale
# Source:  R 4.6.1: src/main/array.c (matprod, crossprod),
#          src/modules/lapack/Lapack.c (La_svd, La_rs), base::scale,
#          stats::quantile (type 7), MASS::ginv, MASS::mvrnorm.
# -----------------------------------------------------------------------------
import ctypes
import os

import numpy as np
from numba import njit
from . import rmath, rng
from .rng import muladd

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")  # as flake.nix sets it for R
_LIBRARY = os.environ.get("ADAMANT_BLAS")
_blas = ctypes.CDLL(_LIBRARY) if _LIBRARY else None


def _integer_type(library):
    """The Fortran INTEGER of the library: 64 bits when OpenBLAS was built
    with INTERFACE64 (its configuration then says USE64BITINT), else 32."""
    if library is None or not hasattr(library, "openblas_get_config"):
        return ctypes.c_int
    library.openblas_get_config.restype = ctypes.c_char_p
    library.openblas_get_config.argtypes = []
    wide = b"USE64BITINT" in library.openblas_get_config()
    return ctypes.c_int64 if wide else ctypes.c_int


_int = _integer_type(_blas)
_int_dtype = np.dtype(_int)
_dbl = ctypes.c_double
_ref = ctypes.byref

# Every argument of the Fortran routines is passed by reference; the lengths
# of the character arguments follow at the end, one size_t each, by value.
_chr, _ptr, _len = ctypes.c_char_p, ctypes.c_void_p, ctypes.c_size_t
if _blas is not None:
    for _name, _characters, _others in [("dgemv_", 1, 10), ("dgemm_", 2, 11),
                                        ("dgesdd_", 1, 13), ("dsyevr_", 3, 18)]:
        _routine = getattr(_blas, _name)
        _routine.argtypes = [_chr] * _characters + [_ptr] * _others + [_len] * _characters
        _routine.restype = None
SQRT_EPS = float(np.sqrt(np.finfo(float).eps))  # sqrt(.Machine$double.eps)


# ----- Matrix products -----

def as_matrix(x):
    """A float matrix in column-major storage; a vector becomes one column.

    Parameters: x (array). Returns: 2-d Fortran-ordered float64 array.
    """
    x = np.asarray(x, dtype=np.float64)
    if x.ndim == 1:
        x = x.reshape(-1, 1)
    return np.asfortranarray(x)


def _pointer(a):
    return a.ctypes.data_as(ctypes.c_void_p)


def _dgemv(trans, a, x):
    """BLAS dgemv: a %*% x (trans "N") or t(a) %*% x (trans "T")."""
    m, n = a.shape
    y = np.empty(m if trans == b"N" else n)
    _blas.dgemv_(trans, _ref(_int(m)), _ref(_int(n)), _ref(_dbl(1.0)), _pointer(a),
                 _ref(_int(m)), _pointer(x), _ref(_int(1)), _ref(_dbl(0.0)), _pointer(y),
                 _ref(_int(1)), ctypes.c_size_t(1))
    return y


def _dgemm(trans_a, a, b, rows, cols, inner):
    """BLAS dgemm: op(a) %*% b, a rows x cols matrix."""
    c = np.empty((rows, cols), order="F")
    _blas.dgemm_(trans_a, b"N", _ref(_int(rows)), _ref(_int(cols)), _ref(_int(inner)),
                 _ref(_dbl(1.0)), _pointer(a), _ref(_int(a.shape[0])), _pointer(b),
                 _ref(_int(b.shape[0])), _ref(_dbl(0.0)), _pointer(c), _ref(_int(rows)),
                 ctypes.c_size_t(1), ctypes.c_size_t(1))
    return c


def _all_finite(*arrays):
    return all(np.isfinite(a).all() for a in arrays)


def matprod(x, y):
    """`x %*% y` as R computes it.

    dgemv when y has one column or x one row, dgemm otherwise; sums in a
    loop when a value is not finite.

    Parameters: x, y (matrices; a vector is one column).
    Returns: Fortran-ordered matrix.
    """
    x, y = as_matrix(x), as_matrix(y)
    if _blas is None or not _all_finite(x, y):
        return np.asfortranarray(_simple_matprod(x, y, rng.FUSED)) if _blas is not None else \
            np.asfortranarray(x @ y)
    if y.shape[1] == 1:
        return _dgemv(b"N", x, y).reshape(-1, 1, order="F")
    if x.shape[0] == 1:
        return _dgemv(b"T", y, x).reshape(1, -1, order="F")
    return _dgemm(b"N", x, y, x.shape[0], y.shape[1], x.shape[1])


def crossprod(x, y):
    """`crossprod(x, y)`, t(x) %*% y, as R computes it.

    Parameters: x, y (matrices; a vector is one column).
    Returns: Fortran-ordered matrix.
    """
    x, y = as_matrix(x), as_matrix(y)
    if _blas is None or not _all_finite(x, y):
        return np.asfortranarray(_simple_matprod(np.asfortranarray(x.T), y, rng.FUSED)) \
            if _blas is not None else np.asfortranarray(x.T @ y)
    if y.shape[1] == 1:
        return _dgemv(b"T", x, y).reshape(-1, 1, order="F")
    if x.shape[1] == 1:
        return _dgemv(b"T", y, x).reshape(1, -1, order="F")
    return _dgemm(b"T", x, y, x.shape[1], y.shape[1], x.shape[0])


def crossprod_columns(X, Y):
    """`apply(Y, 2, function(y) crossprod(X, y))`: one product per column.

    Parameters: X (matrix), Y (matrix with as many rows).
    Returns: Fortran-ordered matrix, ncol(X) x ncol(Y).
    """
    X, Y = as_matrix(X), as_matrix(Y)
    if _blas is None:
        return np.asfortranarray(X.T @ Y)
    out = np.empty((X.shape[1], Y.shape[1]), order="F")
    for b in range(Y.shape[1]):
        out[:, b] = crossprod(X, Y[:, b])[:, 0]
    return out


@njit(cache=True)
def _simple_matprod(x, y, fused):
    """x %*% y with every entry summed in order (R's simple_matprod);
    `fused` adds each product with one rounding."""
    z = np.empty((x.shape[0], y.shape[1]))
    for i in range(x.shape[0]):
        for k in range(y.shape[1]):
            total = 0.0
            for j in range(x.shape[1]):
                if fused:
                    total = muladd(x[i, j], y[j, k], total)
                else:
                    total += x[i, j] * y[j, k]
            z[i, k] = total
    return z


# ----- Decompositions -----

def svd(X):
    """`svd(X)`: LAPACK dgesdd with the workspace R asks for.

    Parameters: X (n x p matrix).
    Returns: dict with `d` (singular values), `u` (n x min(n, p)) and `v`
    (p x min(n, p)).
    """
    X = as_matrix(X)
    n, p = X.shape
    k = min(n, p)
    if _blas is None:
        u, d, vt = np.linalg.svd(X, full_matrices=False)
        return {"d": d, "u": np.asfortranarray(u), "v": np.asfortranarray(vt.T)}
    a = X.copy(order="F")
    d = np.empty(k)
    u = np.empty((n, k), order="F")
    vt = np.empty((k, p), order="F")
    iwork = np.empty(8 * k, dtype=_int_dtype)
    info = _int(0)

    def dgesdd(work, lwork):
        _blas.dgesdd_(b"S", _ref(_int(n)), _ref(_int(p)), _pointer(a), _ref(_int(n)),
                      _pointer(d), _pointer(u), _ref(_int(n)), _pointer(vt), _ref(_int(k)),
                      work, _ref(_int(lwork)), _pointer(iwork), _ref(info), ctypes.c_size_t(1))
        if info.value != 0:
            raise RuntimeError(f"error code {info.value} from Lapack routine 'dgesdd'")

    size = _dbl(0.0)
    dgesdd(_ref(size), -1)
    work = np.empty(int(size.value))
    dgesdd(_pointer(work), len(work))
    return {"d": d, "u": u, "v": np.asfortranarray(vt.T)}


def eigen_symmetric(S):
    """`eigen(S, symmetric = TRUE)`: LAPACK dsyevr as R calls it.

    Parameters: S (symmetric matrix).
    Returns: dict with `values` (decreasing) and `vectors` (columns).
    """
    S = as_matrix(S)
    n = S.shape[0]
    if _blas is None:
        values, vectors = np.linalg.eigh(S)
        return {"values": values[::-1].copy(), "vectors": np.asfortranarray(vectors[:, ::-1])}
    a = S.copy(order="F")
    values = np.empty(n)
    z = np.empty((n, n), order="F")
    isuppz = np.empty(2 * n, dtype=_int_dtype)
    m, info = _int(0), _int(0)

    def dsyevr(work, lwork, iwork, liwork):
        _blas.dsyevr_(b"V", b"A", b"L", _ref(_int(n)), _pointer(a), _ref(_int(n)),
                      _ref(_dbl(0.0)), _ref(_dbl(0.0)), _ref(_int(0)), _ref(_int(0)),
                      _ref(_dbl(0.0)), _ref(m), _pointer(values), _pointer(z), _ref(_int(n)),
                      _pointer(isuppz), work, _ref(_int(lwork)), iwork, _ref(_int(liwork)),
                      _ref(info), ctypes.c_size_t(1), ctypes.c_size_t(1), ctypes.c_size_t(1))
        if info.value != 0:
            raise RuntimeError(f"error code {info.value} from Lapack routine 'dsyevr'")

    size, isize = _dbl(0.0), _int(0)
    dsyevr(_ref(size), -1, _ref(isize), -1)
    work = np.empty(int(size.value))
    iwork = np.empty(isize.value, dtype=_int_dtype)
    dsyevr(_pointer(work), len(work), _pointer(iwork), len(iwork))
    return {"values": values[::-1].copy(), "vectors": np.asfortranarray(z[:, ::-1])}


def ginv(X, tol=SQRT_EPS):
    """`MASS::ginv(X, tol)`: generalized inverse from the SVD.

    Parameters: X (matrix), tol (relative size below which a singular value
    counts as zero). Returns: Fortran-ordered matrix.
    """
    X = as_matrix(X)
    X_svd = svd(X)
    d, u, v = X_svd["d"], X_svd["u"], X_svd["v"]
    positive = d > max(tol * d[0], 0.0)
    if positive.all():
        return matprod(v, (1.0 / d)[:, None] * u.T)
    if not positive.any():
        return np.zeros(X.shape[::-1], order="F")
    return matprod(v[:, positive], (1.0 / d[positive])[:, None] * u[:, positive].T)


# ----- Sums, means and scaling -----

_long = np.longdouble  # the type R adds in (LDOUBLE)
WIDE = np.finfo(_long).nmant > np.finfo(np.float64).nmant  # wider than a double


def _long_sums(A):
    """Every column of A added from the first row down, in long double.

    Parameters: A (matrix, float or long double). Returns: long double array.
    """
    if A.shape[0] == 0:
        return np.zeros(A.shape[1], dtype=_long)
    return np.add.accumulate(np.asarray(A, dtype=_long), axis=0)[-1]


def _long_mean(x):
    """R's two-pass mean of a vector, before it is rounded to a double: the
    sum over the count, plus the mean of the deviations from it.

    Parameters: x (float vector, not empty). Returns: long double.
    """
    x = np.asarray(x, dtype=_long).reshape(-1, 1)
    n = x.shape[0]
    s = _long_sums(x)[0] / n
    if np.isfinite(np.float64(s)):
        s = s + _long_sums(x - s)[0] / n
    return s


def _long_var(x):
    """R's variance of a vector with more than one value: the deviations
    from the mean (a double) and their squares in long double.

    Parameters: x (float vector). Returns: float.
    """
    centre = np.float64(_long_mean(x))
    deviations = np.asarray(x, dtype=_long).reshape(-1, 1) - _long(centre)
    return float(_long_sums(deviations * deviations)[0] / (len(x) - 1))


def column_moments(A):
    """Mean and variance of every column as R's mean() and var() with
    na.rm = TRUE, added in long double.

    Parameters: A (matrix; NaN is skipped).
    Returns: means, variances (float arrays; the mean is NaN for a column
    without values, the variance for a column with fewer than two).
    """
    centre = np.full(A.shape[1], np.nan)
    spread = np.full(A.shape[1], np.nan)
    for j in range(A.shape[1]):
        x = A[:, j]
        x = x[~np.isnan(x)]
        if len(x) > 0:
            centre[j] = np.float64(_long_mean(x))
        if len(x) > 1:
            spread[j] = _long_var(x)
    return centre, spread


@njit(cache=True)
def _col_sums(A):
    out = np.empty(A.shape[1])
    for j in range(A.shape[1]):
        total = 0.0
        for i in range(A.shape[0]):
            total += A[i, j]
        out[j] = total
    return out


def col_sums(A):
    """`colSums(A)`: every column added from the first row down.

    Parameters: A (matrix). Returns: float array of length ncol(A).
    """
    A = np.asarray(A, dtype=np.float64)
    return _long_sums(A).astype(np.float64) if WIDE else _col_sums(A)


def col_means(A):
    """`colMeans(A)`: the sum of every column over the number of rows.

    Parameters: A (matrix). Returns: float array of length ncol(A).
    """
    A = np.asarray(A, dtype=np.float64)
    if WIDE:
        return (_long_sums(A) / A.shape[0]).astype(np.float64)
    return _col_sums(A) / A.shape[0]


def cum_sums(x):
    """`cumsum(x)`: the running sum, rounded to a double at every element.

    Parameters: x (vector). Returns: float array.
    """
    return np.add.accumulate(np.asarray(x, dtype=_long)).astype(np.float64)


def r_sum(x):
    """`sum(x)`: the elements added in order. Returns: float."""
    return float(col_sums(np.asarray(x, dtype=np.float64).reshape(-1, 1))[0])


@njit(cache=True)
def _mean(x):
    n = len(x)
    s = 0.0
    for v in x:
        s += v
    s /= n
    if np.isfinite(s):
        t = 0.0
        for v in x:
            t += v - s
        s += t / n
    return s


def r_mean(x, na_rm=False):
    """`mean(x, na.rm = na_rm)`: the mean with R's second pass.

    Parameters: x (vector), na_rm (drop missing values). Returns: float.
    """
    x = np.asarray(x, dtype=np.float64)
    if na_rm:
        x = x[~np.isnan(x)]
    if len(x) == 0:
        return float("nan")
    return float(_long_mean(x)) if WIDE else float(_mean(x))


@njit(cache=True)
def _var(x, fused):
    n = len(x)
    centre = _mean(x)
    total = 0.0
    for v in x:
        total = muladd(v - centre, v - centre, total) if fused else \
            total + (v - centre) * (v - centre)
    return total / (n - 1)


def r_var(x, na_rm=False):
    """`var(x, na.rm = na_rm)`. Returns: float."""
    x = np.asarray(x, dtype=np.float64)
    if na_rm:
        x = x[~np.isnan(x)]
    if len(x) <= 1:
        return float("nan")
    return _long_var(x) if WIDE else float(_var(x, rng.FUSED))


def r_sd(x, na_rm=False):
    """`sd(x, na.rm = na_rm)`. Returns: float."""
    return float(np.sqrt(r_var(x, na_rm)))


def scale(x):
    """`scale(x)`: columns centred by their mean and divided by their sd.

    Parameters: x (matrix or vector). Returns: Fortran-ordered matrix.
    """
    x = as_matrix(x)
    n = x.shape[0]
    x = x - col_means(x)[None, :]
    spread = np.sqrt(col_sums(x * x) / max(1, n - 1))
    return np.asfortranarray(x / spread[None, :])


def quantile(x, prob):
    """`quantile(x, probs = prob)[[1]]`, type 7. Returns: float."""
    x = np.sort(np.asarray(x, dtype=np.float64))
    index = 1 + max(len(x) - 1, 0) * prob
    lo, hi = int(np.floor(index)), int(np.ceil(index))
    q = x[lo - 1]
    if index > lo and x[hi - 1] != q:
        h = index - lo
        q = (1 - h) * q + h * x[hi - 1]
    return float(q)


def seq_length(start, end, length_out):
    """`seq(start, end, length.out = length_out)`. Returns: float array."""
    if length_out <= 2:
        return np.array([start, end], dtype=np.float64)[:length_out]
    if start == end:
        return np.full(length_out, float(start))
    step = (end - start) / (length_out - 1)
    inner = start + np.arange(1, length_out - 1) * step
    return np.concatenate(([start], inner, [end])).astype(np.float64)


def seq_by(start, end, by):
    """`seq(start, end, by = by)` for numbers that are not integers.

    Returns: float array.
    """
    count = int((end - start) / by + 1e-10)
    x = start + np.arange(count + 1) * by
    return np.minimum(x, end) if by > 0 else np.maximum(x, end)


def which_min_rows(A):
    """`apply(A, 1, which.min)`, zero-based: first smallest value per row.

    Parameters: A (matrix; NaN is skipped). Returns: int array.
    """
    return np.nanargmin(A, axis=1)


def pgamma_upper(x, shape, rate):
    """`pgamma(x, shape, rate, lower.tail = FALSE)`. Returns: float array."""
    return rmath.pgamma_upper(x, shape, rate)


# ----- Multivariate normal -----

def mvrnorm(n, Sigma):
    """`MASS::mvrnorm(n, rep(0, p), Sigma)` from the stream of rng.

    Parameters: n (number of rows), Sigma (p x p covariance).
    Returns: Fortran-ordered matrix, n x p.
    """
    p = Sigma.shape[0]
    eig = eigen_symmetric(Sigma)
    values = eig["values"]
    if not (values >= -1e-6 * abs(values[0])).all():
        raise ValueError("'Sigma' is not positive definite")
    Z = np.asfortranarray(rng.rnorm(p * n).reshape((n, p), order="F"))
    root = np.asfortranarray(np.diag(np.sqrt(np.maximum(values, 0.0))))
    X = matprod(matprod(eig["vectors"], root), np.asfortranarray(Z.T))
    return np.asfortranarray((0.0 + X).T)
