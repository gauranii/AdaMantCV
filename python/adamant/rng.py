# -----------------------------------------------------------------------------
# File:    python/adamant/rng.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. All figures and tables.
# Purpose: The random numbers of R (>= 3.6.0, default kinds): set.seed() with
#          its initial scrambling, the Mersenne-Twister, unif_rand(),
#          sample.int() with "Rejection" sampling and rnorm() by "Inversion",
#          so that a seed gives the same permutations, folds and simulated
#          data as in R. One stream per process, as in R.
# Inputs:  None.
# Outputs: None.
# Usage:   from adamant import rng
#          rng.set_seed(8891); rng.sample_int(10); rng.runif(3); rng.rnorm(3)
# Source:  R 4.6.1: src/main/RNG.c, src/main/random.c (do_sample),
#          src/nmath/snorm.c, src/nmath/qnorm.c.
# -----------------------------------------------------------------------------
import platform

import numpy as np
from llvmlite import ir
from numba import njit
from numba.core import types
from numba.extending import intrinsic



@intrinsic
def muladd(typingctx, a, b, c):
    """a * b + c with a single rounding (fused multiply-add)."""
    sig = types.float64(types.float64, types.float64, types.float64)

    def codegen(context, builder, signature, args):
        double = ir.DoubleType()
        fnty = ir.FunctionType(double, [double, double, double])
        fn = builder.module.declare_intrinsic("llvm.fma", [double], fnty)
        return builder.call(fn, args)

    return sig, codegen


N = 624  # words of the Mersenne-Twister state
M = 397

# The stream: N state words, then the position of the next word.
_state = np.zeros(N + 1, dtype=np.uint32)
_state[N] = N + 1

# R's C code is compiled with one rounding for a product followed by a sum
# on arm64 (fused multiply-add); the polynomials of qnorm() follow it.
FUSED = platform.machine().lower() in ("arm64", "aarch64")


@njit(cache=True)
def _seed(state, seed):
    """Fill the state from a seed as `set.seed()` does.

    Parameters: state (uint32 array of length N + 1), seed (integer).
    Returns: None; `state` is overwritten.
    """
    s = np.uint32(seed)
    for _ in range(50):  # initial scrambling
        s = np.uint32(69069) * s + np.uint32(1)
    for j in range(N + 1):
        s = np.uint32(69069) * s + np.uint32(1)
        if j > 0:
            state[j - 1] = s
    state[N] = N  # the first word is replaced by the position


@njit(cache=True)
def _unif(state):
    """One uniform number in (0, 1): `unif_rand()`.

    Parameters: state (the stream). Returns: float.
    """
    mti = state[N]
    if mti >= N:
        for kk in range(N):
            y = (state[kk] & np.uint32(0x80000000)) | (state[(kk + 1) % N] & np.uint32(0x7FFFFFFF))
            v = state[(kk + M) % N] ^ (y >> np.uint32(1))
            if y & np.uint32(1):
                v ^= np.uint32(0x9908B0DF)
            state[kk] = v
        mti = 0
    y = state[mti]
    state[N] = mti + 1
    y ^= y >> np.uint32(11)
    y ^= (y << np.uint32(7)) & np.uint32(0x9D2C5680)
    y ^= (y << np.uint32(15)) & np.uint32(0xEFC60000)
    y ^= y >> np.uint32(18)
    x = float(y) * 2.3283064365386963e-10
    if x <= 0.0:
        return 0.5 * 2.328306437080797e-10
    if 1.0 - x <= 0.0:
        return 1.0 - 0.5 * 2.328306437080797e-10
    return x


@njit(cache=True)
def _unif_index(state, dn):
    """A uniform index in 0, ..., dn - 1 by rejection: `R_unif_index()`.

    Parameters: state (the stream), dn (number of values). Returns: int.
    """
    if dn <= 0:
        return 0
    bits = int(np.ceil(np.log2(float(dn))))
    while True:
        v = 0
        n = 0
        while n <= bits:
            v = 65536 * v + int(np.floor(_unif(state) * 65536))
            n += 16
        if bits < 64:
            v &= (1 << bits) - 1
        if v < dn:
            return v


@njit(cache=True)
def _sample_int(state, n, size):
    """`sample.int(n, size)` without replacement, values 1 to n.

    Parameters: state (the stream), n, size. Returns: int64 array.
    """
    pool = np.arange(n)
    out = np.empty(size, dtype=np.int64)
    for i in range(size):
        j = _unif_index(state, n)
        out[i] = pool[j] + 1
        n -= 1
        pool[j] = pool[n]
    return out


@njit(cache=True)
def _horner(r, coefficients, fused):
    """Polynomial in r, highest power first, evaluated as qnorm() does.

    Parameters: r, coefficients (array), fused (one rounding per step).
    Returns: float.
    """
    value = coefficients[0]
    for c in coefficients[1:]:
        value = muladd(value, r, c) if fused else value * r + c
    return value


_CENTRE_NUM = np.array([2509.0809287301226727, 33430.575583588128105, 67265.770927008700853,
                        45921.953931549871457, 13731.693765509461125, 1971.5909503065514427,
                        133.14166789178437745, 3.387132872796366608])
_CENTRE_DEN = np.array([5226.495278852854561, 28729.085735721942674, 39307.89580009271061,
                        21213.794301586595867, 5394.1960214247511077, 687.1870074920579083,
                        42.313330701600911252, 1.0])
_MIDDLE_NUM = np.array([7.7454501427834140764e-4, .0227238449892691845833, .24178072517745061177,
                        1.27045825245236838258, 3.64784832476320460504, 5.7694972214606914055,
                        4.6303378461565452959, 1.42343711074968357734])
_MIDDLE_DEN = np.array([1.05075007164441684324e-9, 5.475938084995344946e-4,
                        .0151986665636164571966, .14810397642748007459, .68976733498510000455,
                        1.6763848301838038494, 2.05319162663775882187, 1.0])
_TAIL_NUM = np.array([2.01033439929228813265e-7, 2.71155556874348757815e-5,
                      .0012426609473880784386, .026532189526576123093, .29656057182850489123,
                      1.7848265399172913358, 5.4637849111641143699, 6.6579046435011037772])
_TAIL_DEN = np.array([2.04426310338993978564e-15, 1.4215117583164458887e-7,
                      1.8463183175100546818e-5, 7.868691311456132591e-4,
                      .0148753612908506148525, .13692988092273580531, .59983220655588793769,
                      1.0])


@njit(cache=True)
def _qnorm(p, fused):
    """Standard normal quantile of a probability in (0, 1): `qnorm()`.

    Parameters: p, fused (as `_horner()`). Returns: float.
    """
    q = p - 0.5
    if abs(q) <= 0.425:
        r = muladd(-q, q, 0.180625) if fused else 0.180625 - q * q
        return q * _horner(r, _CENTRE_NUM, fused) / _horner(r, _CENTRE_DEN, fused)
    r = np.sqrt(-np.log(1.0 - p if q > 0 else p))
    if r <= 5.0:
        r += -1.6
        value = _horner(r, _MIDDLE_NUM, fused) / _horner(r, _MIDDLE_DEN, fused)
    else:
        r += -5.0
        value = _horner(r, _TAIL_NUM, fused) / _horner(r, _TAIL_DEN, fused)
    return -value if q < 0.0 else value


@njit(cache=True)
def _rnorm(state, n, mean, sd, fused):
    """`rnorm(n, mean, sd)` by inversion, two uniform numbers per value.

    Parameters: state (the stream), n, mean, sd, fused. Returns: float array.
    """
    big = 134217728.0  # 2^27
    out = np.empty(n)
    for i in range(n):
        u = _unif(state)
        u = int(big * u) + _unif(state)
        z = _qnorm(u / big, fused)
        out[i] = mean + sd * z
    return out


@njit(cache=True)
def _runif(state, n):
    """`runif(n)`. Parameters: state, n. Returns: float array."""
    out = np.empty(n)
    for i in range(n):
        out[i] = _unif(state)
    return out


def set_seed(seed):
    """`set.seed(seed)`: start the stream of this process at a seed.

    Parameters: seed (integer). Returns: None.
    """
    _seed(_state, int(seed) & 0xFFFFFFFF)


def runif(n):
    """`runif(n)`: n uniform numbers in (0, 1). Returns: float array."""
    return _runif(_state, n)


def rnorm(n, mean=0.0, sd=1.0):
    """`rnorm(n, mean, sd)`: n normal numbers. Returns: float array."""
    return _rnorm(_state, n, float(mean), float(sd), FUSED)


def sample_int(n, size=None):
    """`sample.int(n, size)`: `size` of 1, ..., n without replacement.

    Parameters: n, size (n when None). Returns: int64 array, values 1 to n.
    """
    return _sample_int(_state, n, n if size is None else size)


def sample(x):
    """`sample(x)` for a vector of length above one: x in random order.

    Parameters: x (array). Returns: array of the same length.
    """
    x = np.asarray(x)
    return x[sample_int(len(x)) - 1]


def permutation_rows(n, count):
    """`count` consecutive calls of `sample.int(n)`, zero-based, one per row.

    Parameters: n, count. Returns: int64 array, count x n.
    """
    return _permutation_rows(_state, n, count)


@njit(cache=True)
def _permutation_rows(state, n, count):
    out = np.empty((count, n), dtype=np.int64)
    for b in range(count):
        out[b] = _sample_int(state, n, n) - 1
    return out
