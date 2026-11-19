# -----------------------------------------------------------------------------
# File:    python/adamant/glmnet.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Figure 1.
# Purpose: Ridge regression along a given decreasing sequence of penalties and
#          its K-fold cross-validated MSE, computed step for step as
#          glmnet::glmnet() and glmnet::cv.glmnet() do for the gaussian
#          family with alpha = 0 and default options (standardize, intercept,
#          thresh = 1e-7, maxit = 1e5, covariance updating for fewer than 500
#          predictors and naive updating otherwise, folds given by the caller).
# Inputs:  None.
# Outputs: None.
# Usage:   from adamant.glmnet import glmnet_ridge, cv_glmnet_ridge
#          a0, beta = glmnet_ridge(X, y, lambdas)
#          cv = cv_glmnet_ridge(X, y, lambdas, foldid)   # cv["cvm"], cv["lambda_min"]
# Source:  glmnet 5.0 (cv.glmnet, elnet; gaussian family, alpha = 0).
# -----------------------------------------------------------------------------

import numpy as np
from llvmlite import ir
from numba import njit
from numba.core import types
from numba.extending import intrinsic

THRESH = 1e-7  # convergence threshold of glmnet (thresh)
MAXIT = 100000  # largest number of passes over the data for a whole path (maxit)
COVARIANCE_BELOW = 500  # covariance updating for fewer predictors, naive otherwise


@intrinsic
def _muladd(typingctx, a, b, c):
    """a * b + c with a single rounding (LLVM `fma`).

    The compiled glmnet evaluates `g -= c * d` and `a0 -= b * m` in this way
    on machines whose C++ compiler contracts a product and a sum (arm64).

    Parameters: a, b, c floats.  Returns: float.
    """
    sig = types.float64(types.float64, types.float64, types.float64)

    def codegen(context, builder, signature, args):
        double = args[0].type
        fnty = ir.FunctionType(double, [double, double, double])
        fn = builder.module.declare_intrinsic("llvm.fma", [double], fnty)
        return builder.call(fn, args)

    return sig, codegen


@njit(cache=True)
def _dot_columns(A, j, B, k):
    """Inner product of column j of A and column k of B.

    The terms are added in the order of Eigen's vectorised sum with packets
    of two doubles (four running sums, combined at the end, then the odd
    last term), which is what glmnet's C++ code compiles to.

    Parameters: A, B Fortran-ordered matrices with the same number of rows;
    j, k column numbers.  Returns: float.
    """
    n = A.shape[0]
    if n == 0:
        return 0.0
    n2 = (n // 2) * 2
    n4 = (n // 4) * 4
    if n2 == 0:
        return A[0, j] * B[0, k]
    p00 = A[0, j] * B[0, k]
    p01 = A[1, j] * B[1, k]
    if n2 > 2:
        p10 = A[2, j] * B[2, k]
        p11 = A[3, j] * B[3, k]
        for i in range(4, n4, 4):
            p00 += A[i, j] * B[i, k]
            p01 += A[i + 1, j] * B[i + 1, k]
            p10 += A[i + 2, j] * B[i + 2, k]
            p11 += A[i + 3, j] * B[i + 3, k]
        p00 += p10
        p01 += p11
        if n2 > n4:
            p00 += A[n4, j] * B[n4, k]
            p01 += A[n4 + 1, j] * B[n4 + 1, k]
    res = p00 + p01
    for i in range(n2, n):
        res += A[i, j] * B[i, k]
    return res


@njit(cache=True)
def _standardize(x, y):
    """Centre and scale the predictors and the response with weights 1/n.

    glmnet's Standardize1: each column of x and y is centred at its weighted
    mean, multiplied by sqrt(1/n) and divided by its norm, in place.  A
    column that is constant is left out of the fit (ju false).

    Parameters: x (n, p) Fortran-ordered and y (n, 1), both overwritten.
    Returns: ju (p,) bool, xm (p,) means, xs (p,) scales, ym, ys.
    """
    n, p = x.shape
    ju = np.zeros(p, dtype=np.bool_)
    for j in range(p):
        t = x[0, j]
        for i in range(1, n):
            if x[i, j] != t:
                ju[j] = True
                break
    if not ju.any():
        raise ValueError("all predictors are constant")
    w = np.empty((n, 1))
    wsum = 0.0
    for i in range(n):
        wsum += 1.0
    for i in range(n):
        w[i, 0] = 1.0 / wsum
    v = np.sqrt(w)
    xm = np.zeros(p)
    xs = np.zeros(p)
    for j in range(p):
        if not ju[j]:
            continue
        xm[j] = _dot_columns(x, j, w, 0)
        for i in range(n):
            x[i, j] = v[i, 0] * (x[i, j] - xm[j])
        xs[j] = np.sqrt(_dot_columns(x, j, x, j))
    for j in range(p):
        if not ju[j]:
            continue
        for i in range(n):
            x[i, j] /= xs[j]
    ym = _dot_columns(y, 0, w, 0)
    for i in range(n):
        y[i, 0] = v[i, 0] * (y[i, 0] - ym)
    ys = np.sqrt(_dot_columns(y, 0, y, 0))
    for i in range(n):
        y[i, 0] /= ys
    return ju, xm, xs, ym, ys


@njit(cache=True)
def _path_covariance(x, y, ju, vlam, thr, maxit, fused):
    """Coordinate descent along the penalties, covariance updating.

    glmnet's ElnetPath / ElnetPoint for the gaussian family, mode cov, with
    alpha = 0, penalty factors 1, no limits and a given penalty sequence:
    per penalty, passes over the active set from the previous solution, then
    full passes, each followed by passes over the active set, until a full
    pass changes no coefficient by more than the threshold.  The gradient
    X'r is kept up to date through the inner products of the active columns.

    Parameters: x (n, p) and y (n, 1) standardized; ju (p,) columns in the
    fit; vlam penalties on the standardized scale, decreasing; thr, maxit;
    fused whether a product and a sum are rounded once.
    Returns: ao (p, nlam) coefficients in order of entry, ia (p,) the
    columns in order of entry (numbered from 1), kin (nlam,) number entered
    per penalty, lmu number of penalties solved, nlp number of passes, jerr
    (0, or -m when the passes ran out at the m-th penalty).
    """
    n, p = x.shape
    nlam = vlam.shape[0]
    nx = p
    ao = np.zeros((nx, nlam))
    ia = np.zeros(nx, dtype=np.int64)
    kin = np.zeros(nlam, dtype=np.int64)
    mm = np.zeros(p, dtype=np.int64)
    a = np.zeros(p)
    da = np.zeros(p)
    g = np.zeros(p)
    for j in range(p):
        if ju[j]:
            g[j] = _dot_columns(x, j, y, 0)
    ct = np.zeros((nx, p))  # ct[l, j]: inner product of column j and the l-th entered column
    nin = 0
    nlp = 0
    lmu = 0
    iz = False
    for m in range(nlam):
        alm = vlam[m]
        ab = alm * 0.0
        dem = alm * 1.0
        first = True
        while True:
            if (not first) or iz:
                # passes over the active set
                for l in range(nin):
                    da[l] = a[ia[l] - 1]
                iz = True
                while True:
                    nlp += 1
                    dlx = 0.0
                    for l in range(nin):
                        k = ia[l] - 1
                        ak = a[k]
                        u = g[k] + ak * 1.0
                        v = abs(u) - ab
                        a[k] = 0.0
                        if v > 0.0:
                            a[k] = np.copysign(v, u) / (1.0 + dem)
                        if a[k] == ak:
                            continue
                        d = a[k] - ak
                        dlx = max(1.0 * d * d, dlx)
                        row = ct[mm[k] - 1]
                        if fused:
                            for q in range(nin):
                                j = ia[q] - 1
                                g[j] = _muladd(-row[j], d, g[j])
                        else:
                            for q in range(nin):
                                j = ia[q] - 1
                                g[j] -= row[j] * d
                    if dlx < thr:
                        break
                    if nlp > maxit:
                        return ao, ia, kin, lmu, nlp, -(m + 1)
                for l in range(nin):
                    da[l] -= a[ia[l] - 1]
                for j in range(p):
                    if mm[j] != 0 or not ju[j]:
                        continue
                    s = 0.0
                    for l in range(nin):
                        if l == 0:
                            s = da[0] * ct[0, j]
                        else:
                            s += da[l] * ct[l, j]
                    g[j] += s
            first = False
            # full pass
            nlp += 1
            dlx = 0.0
            for k in range(p):
                if not ju[k]:
                    continue
                ak = a[k]
                u = g[k] + ak * 1.0
                v = abs(u) - ab
                a[k] = 0.0
                if v > 0.0:
                    a[k] = np.copysign(v, u) / (1.0 + dem)
                if a[k] == ak:
                    continue
                if mm[k] == 0:
                    nin += 1
                    mm[k] = nin
                    ia[nin - 1] = k + 1
                    row = ct[nin - 1]
                    for j in range(p):
                        if not ju[j]:
                            continue
                        if j == k:
                            row[j] = 1.0
                        elif mm[j] != 0:
                            row[j] = ct[mm[j] - 1, k]
                        else:
                            row[j] = _dot_columns(x, j, x, k)
                d = a[k] - ak
                dlx = max(1.0 * d * d, dlx)
                row = ct[mm[k] - 1]
                if fused:
                    for j in range(p):
                        if ju[j]:
                            g[j] = _muladd(-row[j], d, g[j])
                else:
                    for j in range(p):
                        if ju[j]:
                            g[j] -= row[j] * d
            if dlx < thr:
                break
            if nlp > maxit:
                return ao, ia, kin, lmu, nlp, -(m + 1)
        for l in range(nin):
            ao[l, m] = a[ia[l] - 1]
        kin[m] = nin
        lmu = m + 1
    return ao, ia, kin, lmu, nlp, 0


@njit(cache=True)
def _path_naive(x, y, ju, vlam, thr, maxit):
    """Coordinate descent along the penalties, naive updating.

    glmnet's ElnetPath / ElnetPoint for the gaussian family, mode naive,
    with alpha = 0, penalty factors 1, no limits and a given penalty
    sequence.  The residual (y, overwritten) is updated after every change
    of a coefficient and the gradient of a coordinate is its inner product
    with the residual.  The strong set holds the columns whose absolute
    gradient is positive; a converged full pass is followed by a check of
    the columns outside it.

    Parameters and return value as `_path_covariance` (without `fused`).
    """
    n, p = x.shape
    nlam = vlam.shape[0]
    nx = p
    ao = np.zeros((nx, nlam))
    ia = np.zeros(nx, dtype=np.int64)
    kin = np.zeros(nlam, dtype=np.int64)
    mm = np.zeros(p, dtype=np.int64)
    a = np.zeros(p)
    g = np.zeros(p)
    ix = np.zeros(p, dtype=np.bool_)
    for j in range(p):
        if ju[j]:
            g[j] = abs(_dot_columns(y, 0, x, j))
    nin = 0
    nlp = 0
    lmu = 0
    iz = False
    alm0 = 0.0
    for m in range(nlam):
        alm = vlam[m]
        ab = alm * 0.0
        dem = alm * 1.0
        tlam = 0.0 * (2.0 * alm - alm0)
        alm0 = alm
        for k in range(p):
            if ix[k] or not ju[k]:
                continue
            if g[k] > tlam * 1.0:
                ix[k] = True
        partial = iz
        while True:
            if partial:
                iz = True
                while True:
                    nlp += 1
                    dlx = 0.0
                    for l in range(nin):
                        k = ia[l] - 1
                        ak = a[k]
                        gk = _dot_columns(y, 0, x, k)
                        u = gk + ak * 1.0
                        v = abs(u) - ab
                        a[k] = 0.0
                        if v > 0.0:
                            a[k] = np.copysign(v, u) / (1.0 + dem)
                        if a[k] == ak:
                            continue
                        d = a[k] - ak
                        dlx = max(1.0 * d * d, dlx)
                        for i in range(n):
                            y[i, 0] -= d * x[i, k]
                    if dlx < thr:
                        break
                    if nlp > maxit:
                        return ao, ia, kin, lmu, nlp, -(m + 1)
            # full passes over the strong set until one fails to converge or
            # converges with no column to add to the strong set
            done = False
            while True:
                if nlp > maxit:
                    return ao, ia, kin, lmu, nlp, -(m + 1)
                nlp += 1
                dlx = 0.0
                for k in range(p):
                    if not ix[k]:
                        continue
                    ak = a[k]
                    gk = _dot_columns(y, 0, x, k)
                    u = gk + ak * 1.0
                    v = abs(u) - ab
                    a[k] = 0.0
                    if v > 0.0:
                        a[k] = np.copysign(v, u) / (1.0 + dem)
                    if a[k] == ak:
                        continue
                    if mm[k] == 0:
                        nin += 1
                        mm[k] = nin
                        ia[nin - 1] = k + 1
                    d = a[k] - ak
                    dlx = max(1.0 * d * d, dlx)
                    for i in range(n):
                        y[i, 0] -= d * x[i, k]
                if dlx < thr:
                    updated = False
                    for k in range(p):
                        if ix[k] or not ju[k]:
                            continue
                        g[k] = abs(_dot_columns(y, 0, x, k))
                        if g[k] > ab * 1.0:
                            ix[k] = True
                            updated = True
                    if not updated:
                        done = True
                        break
                    continue
                if nlp > maxit:
                    return ao, ia, kin, lmu, nlp, -(m + 1)
                break
            if done:
                break
            partial = True
        for l in range(nin):
            ao[l, m] = a[ia[l] - 1]
        kin[m] = nin
        lmu = m + 1
    return ao, ia, kin, lmu, nlp, 0


@njit(cache=True)
def _unstandardize(ao, ia, kin, lmu, xm, xs, ym, ys, p, fused):
    """Coefficients and intercepts on the original scale.

    Parameters: output of a path function, the means and scales of
    `_standardize`, the number of predictors, and `fused` as in
    `_path_covariance`.  Returns: a0 (lmu,), beta (p, lmu).
    """
    a0 = np.zeros(lmu)
    beta = np.zeros((p, lmu))
    for k in range(lmu):
        nk = kin[k]
        for l in range(nk):
            ao[l, k] *= ys / xs[ia[l] - 1]
        s = 0.0
        for l in range(nk):
            if fused:
                s = _muladd(-ao[l, k], xm[ia[l] - 1], s)
            else:
                s -= ao[l, k] * xm[ia[l] - 1]
        s += ym
        a0[k] = s
        for l in range(nk):
            beta[ia[l] - 1, k] = ao[l, k]
    return a0, beta


@njit(cache=True)
def _predict(x, rows, a0, beta, fused):
    """Linear predictor a0 + x beta of some rows for every penalty.

    The terms are added in the order of the predictors, as the product of a
    dense and a sparse matrix in R's Matrix package does.

    Parameters: x (n, p); rows the row numbers; a0 (nlam,); beta (p, nlam);
    fused as in `_path_covariance`.  Returns: (len(rows), nlam).
    """
    p, nlam = beta.shape
    out = np.empty((rows.shape[0], nlam))
    for r in range(rows.shape[0]):
        i = rows[r]
        for l in range(nlam):
            s = a0[l]
            if fused:
                for j in range(p):
                    s = _muladd(beta[j, l], x[i, j], s)
            else:
                for j in range(p):
                    s += x[i, j] * beta[j, l]
            out[r, l] = s
    return out


def _fused_default():
    """Whether glmnet's compiler rounds a product and a sum once (arm64)."""
    import platform

    return platform.machine().lower() in ("arm64", "aarch64")


def _fit(X, y, lambdas, thresh=THRESH, maxit=MAXIT, fused=None):
    """One glmnet fit; see `glmnet_ridge`.

    Returns: dict with a0 (lmu,), beta (p, lmu), lambda (lmu,) (the
    penalties as glmnet reports them, that is divided and multiplied back by
    the scale of the response), npasses and jerr.  lmu is the number of
    penalties, fewer than given only when the passes ran out (jerr < 0).
    """
    if fused is None:
        fused = _fused_default()
    x = np.array(X, dtype=np.float64, order="F")
    yy = np.array(y, dtype=np.float64).reshape(-1, 1).copy(order="F")
    n, p = x.shape
    if yy.shape[0] != n:
        raise ValueError("X and y differ in the number of observations")
    ulam = np.sort(np.asarray(lambdas, dtype=np.float64).ravel())[::-1].copy()
    if np.any(ulam < 0):
        raise ValueError("lambdas should be non-negative")
    if np.all(yy == yy[0, 0]):
        raise ValueError("y is constant")
    ju, xm, xs, ym, ys = _standardize(x, yy)
    vlam = ulam / ys
    if p < COVARIANCE_BELOW:
        ao, ia, kin, lmu, nlp, jerr = _path_covariance(x, yy, ju, vlam, thresh, maxit, fused)
    else:
        ao, ia, kin, lmu, nlp, jerr = _path_naive(x, yy, ju, vlam, thresh, maxit)
    if lmu < 1:
        raise RuntimeError("no penalty converged within maxit passes")
    a0, beta = _unstandardize(ao, ia, kin, lmu, xm, xs, ym, ys, p, fused)
    return dict(a0=a0, beta=beta, **{"lambda": vlam[:lmu] * ys}, npasses=nlp, jerr=jerr)


def glmnet_ridge(X, y, lambdas):
    """Ridge regression path, as coef(glmnet(X, y, alpha = 0, lambda = lambdas)).

    Gaussian family, standardized predictors, intercept, thresh = 1e-7,
    maxit = 1e5, covariance updating when X has fewer than 500 columns and
    naive updating otherwise.

    Parameters
    ----------
    X : (n, p) array of predictors.
    y : (n,) response.
    lambdas : penalties; used in decreasing order.

    Returns
    -------
    a0 : (nlam,) intercepts, one per penalty in decreasing order.
    beta : (p, nlam) coefficients on the original scale of X and y.
    Both have fewer columns when the passes ran out before the last penalty
    (glmnet then returns the solutions of the larger penalties only).
    """
    fit = _fit(X, y, lambdas)
    return fit["a0"], fit["beta"]


def _approx(x, y, xout):
    """Linear interpolation as R's approx(x, y, xout) for increasing x.

    Tied x are replaced by one point with the mean of their y.
    """
    ux, inverse = np.unique(x, return_inverse=True)
    uy = np.bincount(inverse, weights=y) / np.bincount(inverse)
    n = ux.shape[0]
    out = np.empty(xout.shape[0])
    for q, v in enumerate(xout):
        i, j = 0, n - 1
        if v < ux[i] or v > ux[j]:
            out[q] = np.nan
            continue
        while i < j - 1:
            ij = (i + j) // 2
            if v < ux[ij]:
                j = ij
            else:
                i = ij
        if v == ux[j]:
            out[q] = uy[j]
        elif v == ux[i]:
            out[q] = uy[i]
        else:
            out[q] = uy[i] + (uy[j] - uy[i]) * ((v - ux[i]) / (ux[j] - ux[i]))
    return out


def _lambda_interp(lam, s):
    """Position of the penalties s in the sequence lam (glmnet's lambda.interp).

    Parameters: lam decreasing penalties of a fit; s penalties wanted.
    Returns: left, right (indices from 0) and frac, so that the coefficients
    at s are frac * coef[left] + (1 - frac) * coef[right].
    """
    if lam.shape[0] == 1:
        z = np.zeros(s.shape[0], dtype=np.int64)
        return z, z, np.ones(s.shape[0])
    k = lam.shape[0]
    sfrac = (lam[0] - s) / (lam[0] - lam[k - 1])
    lamn = (lam[0] - lam) / (lam[0] - lam[k - 1])
    sfrac = np.where(sfrac < lamn.min(), lamn.min(), sfrac)
    sfrac = np.where(sfrac > lamn.max(), lamn.max(), sfrac)
    coord = _approx(lamn, np.arange(1.0, k + 1.0), sfrac)
    left = np.floor(coord).astype(np.int64) - 1
    right = np.ceil(coord).astype(np.int64) - 1
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = (sfrac - lamn[right]) / (lamn[left] - lamn[right])
    frac[left == right] = 1.0
    frac[np.abs(lamn[left] - lamn[right]) < np.finfo(np.float64).eps] = 1.0
    return left, right, frac


def _weighted_mean(values, w):
    """Column-wise weighted.mean(values, w) of R: sum(values * w) / sum(w).

    The sums run down the rows in R's extended precision (long double).
    """
    ld = np.longdouble
    wsum = ld(0)
    for wi in w:
        wsum += ld(wi)
    acc = np.zeros(values.shape[1], dtype=ld)
    for i in range(values.shape[0]):
        acc += (values[i] * w[i]).astype(ld)
    return acc.astype(np.float64) / np.float64(wsum)


def cv_glmnet_ridge(X, y, lambdas, foldid):
    """Cross-validated MSE of ridge regression, as cv.glmnet(X, y, alpha = 0,
    lambda = lambdas, foldid = foldid).

    One fit on all observations and one per fold on the others, all with
    the given penalties; each observation is predicted from the fit without
    its fold at the penalties of the full fit; the squared errors are
    averaged within folds and the fold means averaged with the fold sizes
    as weights (over all observations when a fold would have fewer than
    three on average, as cv.glmnet does).

    Parameters
    ----------
    X : (n, p) array of predictors.
    y : (n,) response.
    lambdas : penalties; used in decreasing order.
    foldid : (n,) integers 1..nfolds, the fold of each observation.

    Returns
    -------
    dict with
      cvm : (nlam,) cross-validated MSE per penalty (decreasing order);
      lambda_min : the largest penalty at which cvm is smallest, as reported
        by glmnet (the given value divided and multiplied back by the scale
        of the response, so it may differ from it in the last digit);
      lambda : (nlam,) the penalties as reported by glmnet;
      index_min : position of lambda_min in `lambda` (from 0);
      npasses : passes over the data of the full fit;
      fold_npasses : (nfolds,) passes over the data of each fold fit.
    """
    X = np.asarray(X, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64).ravel()
    foldid = np.asarray(foldid).astype(np.int64).ravel()
    n = X.shape[0]
    if foldid.shape[0] != n or foldid.min() < 1:
        raise ValueError("foldid must give a fold 1..nfolds for every observation")
    nfolds = int(foldid.max())
    fused = _fused_default()
    full = _fit(X, y, lambdas)
    lam = full["lambda"]
    nlam = lam.shape[0]
    xc = np.ascontiguousarray(X)
    predmat = np.full((n, nlam), np.nan)
    fold_npasses = []
    for i in range(1, nfolds + 1):
        which = foldid == i
        fit = _fit(X[~which], y[~which], lambdas)
        fold_npasses.append(fit["npasses"])
        left, right, frac = _lambda_interp(fit["lambda"], lam)
        a0 = fit["a0"][left] * frac + fit["a0"][right] * (1.0 - frac)
        beta = fit["beta"][:, left] * frac + fit["beta"][:, right] * (1.0 - frac)
        predmat[which] = _predict(xc, np.flatnonzero(which), a0, np.ascontiguousarray(beta),
                                  fused)
    cvraw = (y[:, None] - predmat) ** 2
    weights = np.ones(n)
    if n / nfolds < 3:
        cvm = _weighted_mean(cvraw, weights)
    else:
        folds = [i for i in range(1, nfolds + 1) if np.any(foldid == i)]
        outmat = np.empty((len(folds), nlam))
        wisum = np.empty(len(folds))
        for r, i in enumerate(folds):
            which = foldid == i
            outmat[r] = _weighted_mean(cvraw[which], weights[which])
            wisum[r] = weights[which].sum()
        cvm = _weighted_mean(outmat, wisum)
    idmin = cvm <= np.nanmin(cvm)
    lambda_min = float(np.max(lam[idmin]))
    index_min = int(np.flatnonzero(lam == lambda_min)[0])
    return dict(cvm=cvm, lambda_min=lambda_min, **{"lambda": lam},
                index_min=index_min, npasses=full["npasses"],
                fold_npasses=np.array(fold_npasses))
