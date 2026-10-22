"""Model points of the fitted state models by component-stratified scrambled Sobol sequences (derivation.md
Section 4.3, eq. (26)), and weighted distribution metrics of the rebuilt electrical quantities (eq. (29))."""

import warnings

import numpy as np
from scipy.special import gammaincinv, ndtri
from scipy.stats import qmc

from . import mixture as mx

EPS = np.finfo(float).eps


def allocation(weights, count):
    """N_k^s = N^s pi_k rounded by largest remainder, at least one point per component."""
    quota = count * np.asarray(weights)
    alloc = np.maximum(1, np.floor(quota)).astype(int)
    while alloc.sum() < count:
        alloc[np.argmax(quota - alloc)] += 1
    while alloc.sum() > count:
        c = np.flatnonzero(alloc > 1)
        alloc[c[np.argmin(quota[c] - alloc[c])]] -= 1
    return alloc


def _sobol(dim, n, rng):
    eng = qmc.Sobol(dim, scramble=True, seed=rng)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        u = eng.random(n)
    return np.clip(u, EPS, 1 - EPS)


def _strength(u, nu):
    return np.ones(len(u)) if np.isinf(nu) else 2.0 * gammaincinv(nu / 2, u) / nu


def sample_component(model, k, n, rng):
    """Standardized points z of component k, eq. (26)."""
    if isinstance(model, mx.Mixture):
        d = model.means.shape[1]
        heavy = not np.isinf(model.df[k])
        u = _sobol(d + int(heavy), n, rng)
        L = np.linalg.cholesky(model.scatter[k])
        g = ndtri(u[:, :d]) @ L.T
        if heavy:
            g /= np.sqrt(_strength(u[:, d], model.df[k]))[:, None]
        return model.means[k] + g
    if isinstance(model, mx.BlockTMM):
        r, dw = model.means.shape[1], model.intercept.shape[1]
        u = _sobol(r + dw + 1, n, rng)
        U = model.means[k] + (ndtri(u[:, :r]) @ np.linalg.cholesky(model.scatter[k]).T) / np.sqrt(_strength(u[:, -1], model.df[k]))[:, None]
        W = model.intercept[k] + U @ model.coefficient[k].T + ndtri(u[:, r:r + dw]) @ np.linalg.cholesky(model.residual[k]).T
        return U @ model.Qu.T + W @ model.Qw.T
    if isinstance(model, mx.TwoBlockTMM):
        r1, r2, dw = model.means.shape[1], model.intercept2.shape[1], model.intercept.shape[1]
        u = _sobol(r1 + r2 + dw + 2, n, rng)
        U1 = model.means[k] + (ndtri(u[:, :r1]) @ np.linalg.cholesky(model.scatter[k]).T) / np.sqrt(_strength(u[:, -2], model.df[k]))[:, None]
        U2 = (model.intercept2[k] + U1 @ model.coefficient2[k].T
              + (ndtri(u[:, r1:r1 + r2]) @ np.linalg.cholesky(model.scatter2[k]).T) / np.sqrt(_strength(u[:, -1], model.df2[k]))[:, None])
        U = np.hstack([U1, U2])
        W = model.intercept[k] + U @ model.coefficient[k].T + ndtri(u[:, r1 + r2:r1 + r2 + dw]) @ np.linalg.cholesky(model.residual[k]).T
        return U1 @ model.Q1.T + U2 @ model.Q2.T + W @ model.Qw.T
    raise TypeError(type(model))


def model_points(model, count, entropy):
    """All points z (count, d) and weights pi_k / N_k^s. Component k uses its own scrambled Sobol sequence
    seeded by SeedSequence(entropy).spawn(K)[k]."""
    alloc = allocation(model.weights, count)
    children = np.random.SeedSequence(entropy).spawn(model.K)
    Z, w = [], []
    for k in range(model.K):
        Z.append(sample_component(model, k, int(alloc[k]), np.random.default_rng(children[k])))
        w.append(np.full(alloc[k], model.weights[k] / alloc[k]))
    return np.vstack(Z), np.concatenate(w), alloc


def evaluate_points(Z, coords, outmap, keys, chunk=1 << 16):
    """Electrical quantities of model points, computed in chunks; returns dict of (count, n_obj) arrays."""
    res = {}
    for s in range(0, len(Z), chunk):
        q = outmap.evaluate(coords.to_x(Z[s:s + chunk]))
        for k in keys:
            if k not in res:
                res[k] = np.empty((len(Z), q[k].shape[1]))
            res[k][s:s + chunk] = q[k]
    return res


# ---------------------------------------------------------------- distribution metrics

PROBS = (0.001, 0.01, 0.99, 0.999)


def weighted_quantiles(Y, w, probs):
    """Y (n, m), weights w (n,) summing to 1. Returns (len(probs), m) with the lower quantile convention."""
    order = np.argsort(Y, axis=0)
    Ys = np.take_along_axis(Y, order, axis=0)
    cw = np.cumsum(w[order], axis=0)
    out = np.empty((len(probs), Y.shape[1]))
    for i, p in enumerate(probs):
        idx = np.minimum((cw < p).sum(axis=0), len(Y) - 1)
        out[i] = Ys[idx, np.arange(Y.shape[1])]
    return out, Ys, cw


def ks_distance(Ys, cw, R):
    """sup |F_model - F_ref|: Ys, cw sorted model values and cumulative weights per column; R reference sample."""
    n = len(R)
    out = np.empty(R.shape[1])
    for j in range(R.shape[1]):
        r = np.sort(R[:, j])
        Fm = np.concatenate([[0.0], cw[:, j]])
        lo = Fm[np.searchsorted(Ys[:, j], r, side="left")]
        hi = Fm[np.searchsorted(Ys[:, j], r, side="right")]
        Fr_hi = np.arange(1, n + 1) / n
        Fr_lo = np.arange(0, n) / n
        out[j] = max(np.abs(hi - Fr_hi).max(), np.abs(lo - Fr_lo).max())
    return out


def distribution_metrics(Y, w, R):
    """Mean error, std relative error, quantile errors (physical and in reference std) and KS distance of the
    weighted sample (Y, w) against the reference sample R (equal weights)."""
    w = w / w.sum()
    mu = w @ Y
    sd = np.sqrt(np.maximum(w @ (Y - mu) ** 2, 0))
    q, Ys, cw = weighted_quantiles(Y, w, PROBS)
    Rm, Rs = R.mean(axis=0), R.std(axis=0)
    Rq = np.quantile(R, PROBS, axis=0, method="inverted_cdf")
    out = {"mean_ref": Rm, "std_ref": Rs, "mean_error": mu - Rm, "std_rel_error": sd / Rs - 1, "ks": ks_distance(Ys, cw, R)}
    for p, a, b in zip(PROBS, q, Rq):
        tag = f"q{p * 100:g}".replace(".", "p")
        out[f"{tag}_ref"] = b
        out[f"{tag}_error"] = a - b
        out[f"{tag}_error_in_std"] = (a - b) / Rs
    return out
