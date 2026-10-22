"""Joint probability models of the power-flow state (01_IDEA/derivation.md Section 3).

Coordinates: eq. (11) standardization by the eigen-decomposition of the modeling-scenario covariance.
Models: GMM, standard TMM eqs. (12)-(14) fitted by EM eqs. (19)-(20), and the block TMM eqs. (15)-(18)
fitted by ECM eqs. (19)-(21). Settings follow Section 3.4 and the MATLAB reference
+TMM_Transmission_core/+statistics/fit_block_tmm.m: covariance eigenvalue floor 1e-7, nu in [4.1, 1e5],
initial splits by the projection on principal directions 1, 2, 3, GMM first, TMM and block TMM started
from the GMM with nu0 in {12, 200}; a fit is accepted when it converged, the mean log-likelihood never
decreased, and every component keeps N_k > d + 2.
"""

import time
from dataclasses import dataclass, field

import numpy as np
from scipy import linalg, optimize
from scipy.special import digamma, gammaln

COV_FLOOR = 1e-7
DF_BOUNDS = (4.1, 1e5)
LL_TOL = 1e-6
MAX_ITER = 3000
DECREASE_TOL = 1e-8
INITIAL_DF = (12.0, 200.0)
RESTART_DIRECTIONS = (0, 1, 2)


@dataclass
class Coordinates:
    center: np.ndarray
    basis: np.ndarray
    scale: np.ndarray
    eigenvalues: np.ndarray
    loss_tolerance: float
    resolvable_tolerance: float

    @property
    def d(self):
        return self.basis.shape[1]

    def to_z(self, X):
        return ((X - self.center) @ self.basis) / self.scale

    def to_x(self, Z):
        return (Z * self.scale) @ self.basis.T + self.center


def fit_coordinates(X, loss_tolerance=1e-10):
    """Eq. (11). d = number of leading directions such that the dropped share of total variance is at most
    loss_tolerance, limited to eigenvalues above p * eps * lambda_1 (numerically resolvable)."""
    center = X.mean(axis=0)
    C = np.cov(X, rowvar=False)
    w, E = linalg.eigh(C)
    order = np.argsort(w)[::-1]
    w, E = w[order], E[:, order]
    wpos = np.clip(w, 0.0, None)
    dropped = np.concatenate([np.cumsum(wpos[::-1])[::-1], [0.0]]) / wpos.sum()
    d_loss = int(np.argmax(dropped <= loss_tolerance))
    res_tol = X.shape[1] * np.finfo(float).eps * w[0]
    d_res = int((w > res_tol).sum())
    d = min(d_loss, d_res)
    return Coordinates(center, E[:, :d], np.sqrt(w[:d]), w, loss_tolerance, res_tol)


def floor_cov(A, floor):
    A = 0.5 * (A + A.T)
    w, V = linalg.eigh(A)
    return (V * np.maximum(w, floor)) @ V.T


def logsumexp_rows(L):
    peak = L.max(axis=1, keepdims=True)
    return (peak + np.log(np.exp(L - peak).sum(axis=1, keepdims=True)))[:, 0]


@dataclass
class Mixture:
    """GMM (df = inf) or standard TMM in standardized coordinates."""
    weights: np.ndarray
    means: np.ndarray
    scatter: np.ndarray
    df: np.ndarray
    family: str = "tmm"
    audit: dict = field(default_factory=dict)

    @property
    def K(self):
        return len(self.weights)


def mixture_terms(Z, model):
    """Per-component log(pi_k) + log density and Mahalanobis distances delta_ik."""
    n, d = Z.shape
    logs = np.empty((n, model.K))
    delta = np.empty((n, model.K))
    for k in range(model.K):
        L = linalg.cholesky(model.scatter[k], lower=True)
        y = linalg.solve_triangular(L, (Z - model.means[k]).T, lower=True)
        dk = np.einsum("ij,ij->j", y, y)
        logdet = np.log(np.diag(L)).sum()
        nu = model.df[k]
        if np.isinf(nu):
            lp = -0.5 * (d * np.log(2 * np.pi) + dk)
        else:
            lp = gammaln((nu + d) / 2) - gammaln(nu / 2) - d / 2 * np.log(nu * np.pi) - (nu + d) / 2 * np.log1p(dk / nu)
        logs[:, k] = np.log(model.weights[k]) + lp - logdet
        delta[:, k] = dk
    return logs, delta


def log_density(Z, model):
    if isinstance(model, BlockTMM):
        return block_terms(Z, model)[0]
    if isinstance(model, TwoBlockTMM):
        return two_block_terms(Z, model)[0]
    logs, _ = mixture_terms(Z, model)
    return logsumexp_rows(logs)


def _solve_df(gamma_k, tau_k, delta_k, nu_old, dim, Nk):
    stat = gamma_k @ (digamma((nu_old + dim) / 2) - np.log((nu_old + delta_k) / 2) - tau_k) / Nk

    def score(v):
        return np.log(v / 2) - digamma(v / 2) + 1 + stat

    lo, hi = DF_BOUNDS
    if score(lo) <= 0:
        return lo
    if score(hi) >= 0:
        return hi
    return optimize.brentq(score, lo, hi, xtol=1e-10, rtol=1e-12)


def _em(Z, model, update_df):
    """EM for GMM (update_df False, df inf) or TMM. Returns model with audit."""
    t0 = time.time()
    n, d = Z.shape
    logs, delta = mixture_terms(Z, model)
    ll = logsumexp_rows(logs)
    hist = [ll.mean()]
    gamma = np.exp(logs - ll[:, None])
    status = "maximum-iterations"
    it = 0
    for it in range(1, MAX_ITER + 1):
        Nk = gamma.sum(axis=0)
        if Nk.min() <= d + 2:
            status = "collapsed-component"
            break
        nxt = Mixture(Nk / n, model.means.copy(), model.scatter.copy(), model.df.copy(), model.family)
        for k in range(model.K):
            nu = model.df[k]
            tau = np.ones(n) if np.isinf(nu) else (nu + d) / (nu + delta[:, k])
            wgt = gamma[:, k] * tau
            mu = wgt @ Z / wgt.sum()
            D = Z - mu
            S = (D * wgt[:, None]).T @ D / Nk[k]
            floor = COV_FLOOR if np.isinf(nu) else COV_FLOOR * (nu - 2) / nu
            nxt.means[k] = mu
            nxt.scatter[k] = floor_cov(S, floor)
            if update_df:
                nxt.df[k] = _solve_df(gamma[:, k], tau, delta[:, k], nu, d, Nk[k])
        logs_n, delta_n = mixture_terms(Z, nxt)
        ll_n = logsumexp_rows(logs_n)
        value = ll_n.mean()
        gain = value - hist[-1]
        if gain < -DECREASE_TOL:
            status = "rejected-decreasing-update"
            break
        model, delta, gamma = nxt, delta_n, np.exp(logs_n - ll_n[:, None])
        hist.append(value)
        if abs(gain) < LL_TOL:
            status = "converged"
            break
    Nk = gamma.sum(axis=0)
    model.audit = {"status": status, "iterations": it, "mean_log_likelihood": float(hist[-1]),
                   "effective_counts": Nk.tolist(), "seconds": time.time() - t0,
                   "accepted": status == "converged" and Nk.min() > d + 2}
    return model


def initial_gmm(Z, K, direction):
    """Sort the modeling scenarios by their projection on principal direction `direction` (the z coordinate)
    and split them into K equal-count groups."""
    n, d = Z.shape
    order = np.argsort(Z[:, direction], kind="stable")
    groups = np.array_split(order, K)
    means = np.stack([Z[g].mean(axis=0) for g in groups])
    scatter = np.stack([floor_cov(np.cov(Z[g], rowvar=False, bias=True), COV_FLOOR) for g in groups])
    weights = np.array([len(g) for g in groups], float) / n
    return Mixture(weights, means, scatter, np.full(K, np.inf), "gmm")


def fit_gmm(Z, K, direction):
    m = _em(Z, initial_gmm(Z, K, direction), update_df=False)
    m.audit.update({"method": "GMM", "K": K, "restart_direction": direction + 1})
    return m


def fit_tmm(Z, gmm, nu0):
    start = Mixture(gmm.weights.copy(), gmm.means.copy(), gmm.scatter * (nu0 - 2) / nu0, np.full(gmm.K, nu0), "tmm")
    m = _em(Z, start, update_df=True)
    m.audit.update({"method": "TMM", "K": gmm.K, "restart_direction": gmm.audit["restart_direction"], "initial_df": nu0})
    return m


# ---------------------------------------------------------------- block TMM, eqs. (15)-(18), (21)

@dataclass
class BlockTMM:
    weights: np.ndarray
    means: np.ndarray        # (K, r)
    scatter: np.ndarray      # (K, r, r)
    df: np.ndarray
    intercept: np.ndarray    # (K, d - r)
    coefficient: np.ndarray  # (K, d - r, r)
    residual: np.ndarray     # (K, d - r, d - r)
    Qu: np.ndarray
    Qw: np.ndarray
    family: str = "block-tmm"
    audit: dict = field(default_factory=dict)

    @property
    def K(self):
        return len(self.weights)


def block_terms(Z, model):
    U, W = Z @ model.Qu, Z @ model.Qw
    n, r = U.shape
    dw = W.shape[1]
    logs = np.empty((n, model.K))
    delta = np.empty((n, model.K))
    for k in range(model.K):
        L = linalg.cholesky(model.scatter[k], lower=True)
        y = linalg.solve_triangular(L, (U - model.means[k]).T, lower=True)
        dk = np.einsum("ij,ij->j", y, y)
        nu = model.df[k]
        if np.isinf(nu):
            tl = -0.5 * (r * np.log(2 * np.pi) + dk)
        else:
            tl = gammaln((nu + r) / 2) - gammaln(nu / 2) - r / 2 * np.log(nu * np.pi) - (nu + r) / 2 * np.log1p(dk / nu)
        R = linalg.cholesky(model.residual[k], lower=True)
        e = linalg.solve_triangular(R, (W - model.intercept[k] - U @ model.coefficient[k].T).T, lower=True)
        gl = -0.5 * (dw * np.log(2 * np.pi) + np.einsum("ij,ij->j", e, e)) - np.log(np.diag(R)).sum()
        logs[:, k] = np.log(model.weights[k]) + tl - np.log(np.diag(L)).sum() + gl
        delta[:, k] = dk
    ll = logsumexp_rows(logs)
    return ll, np.exp(logs - ll[:, None]), delta


def fit_block_tmm(Z, gmm, nu0, Qu, Qw):
    t0 = time.time()
    U, W = Z @ Qu, Z @ Qw
    n, r = U.shape
    dw = W.shape[1]
    K = gmm.K
    m = BlockTMM(gmm.weights.copy(), gmm.means @ Qu, np.zeros((K, r, r)), np.full(K, nu0), np.zeros((K, dw)),
                 np.zeros((K, dw, r)), np.zeros((K, dw, dw)), Qu, Qw)
    for k in range(K):
        C = gmm.scatter[k]
        Suu, Swu, Sww = Qu.T @ C @ Qu, Qw.T @ C @ Qu, Qw.T @ C @ Qw
        G = linalg.solve(Suu, Swu.T, assume_a="pos").T
        m.scatter[k] = floor_cov(Suu * (nu0 - 2) / nu0, COV_FLOOR * (nu0 - 2) / nu0)
        m.coefficient[k] = G
        m.intercept[k] = gmm.means[k] @ Qw - G @ (gmm.means[k] @ Qu)
        m.residual[k] = floor_cov(Sww - G @ Suu @ G.T, COV_FLOOR)
    ll, gamma, delta = block_terms(Z, m)
    hist = [ll.mean()]
    design = np.hstack([np.ones((n, 1)), U])
    status = "maximum-iterations"
    it = 0
    for it in range(1, MAX_ITER + 1):
        Nk = gamma.sum(axis=0)
        if Nk.min() <= dw + r + 2:
            status = "collapsed-component"
            break
        nxt = BlockTMM(Nk / n, m.means.copy(), m.scatter.copy(), m.df.copy(), m.intercept.copy(), m.coefficient.copy(),
                       m.residual.copy(), Qu, Qw)
        for k in range(K):
            nu = m.df[k]
            tau = (nu + r) / (nu + delta[:, k])
            wgt = gamma[:, k] * tau
            mu = wgt @ U / wgt.sum()
            D = U - mu
            nxt.means[k] = mu
            nxt.scatter[k] = floor_cov((D * wgt[:, None]).T @ D / Nk[k], COV_FLOOR * (nu - 2) / nu)
            nxt.df[k] = _solve_df(gamma[:, k], tau, delta[:, k], nu, r, Nk[k])
            Xg = design * gamma[:, k][:, None]
            b = linalg.solve(design.T @ Xg, Xg.T @ W, assume_a="pos")
            nxt.intercept[k] = b[0]
            nxt.coefficient[k] = b[1:].T
            e = W - design @ b
            nxt.residual[k] = floor_cov((e * gamma[:, k][:, None]).T @ e / Nk[k], COV_FLOOR)
        ll_n, gamma_n, delta_n = block_terms(Z, nxt)
        value = ll_n.mean()
        gain = value - hist[-1]
        if gain < -DECREASE_TOL:
            status = "rejected-decreasing-update"
            break
        m, gamma, delta = nxt, gamma_n, delta_n
        hist.append(value)
        if abs(gain) < LL_TOL:
            status = "converged"
            break
    Nk = gamma.sum(axis=0)
    m.audit = {"method": "block TMM", "K": K, "initial_df": nu0, "restart_direction": gmm.audit["restart_direction"],
               "status": status, "iterations": it, "mean_log_likelihood": float(hist[-1]),
               "effective_counts": Nk.tolist(), "seconds": time.time() - t0,
               "accepted": status == "converged" and Nk.min() > dw + r + 2}
    return m


# ---------------------------------------------------------------- two-strength block TMM (wind and PV)

@dataclass
class TwoBlockTMM:
    """u1 (wind-driven directions) | k ~ t_nu1(m1, S1); u2 (PV-driven directions) | u1, k ~ t_nu2(a2 + G2 u1, O2)
    with an independent strength; w | u1, u2, k ~ N(a + G [u1, u2], Omega)."""
    weights: np.ndarray
    means: np.ndarray        # (K, r1)
    scatter: np.ndarray      # (K, r1, r1)
    df: np.ndarray           # (K,) nu1
    intercept2: np.ndarray   # (K, r2)
    coefficient2: np.ndarray # (K, r2, r1)
    scatter2: np.ndarray     # (K, r2, r2)
    df2: np.ndarray          # (K,) nu2
    intercept: np.ndarray    # (K, dw)
    coefficient: np.ndarray  # (K, dw, r1 + r2)
    residual: np.ndarray     # (K, dw, dw)
    Q1: np.ndarray
    Q2: np.ndarray
    Qw: np.ndarray
    family: str = "two-block-tmm"
    audit: dict = field(default_factory=dict)

    @property
    def K(self):
        return len(self.weights)


def _t_log(dk, nu, dim, logdet):
    return gammaln((nu + dim) / 2) - gammaln(nu / 2) - dim / 2 * np.log(nu * np.pi) - (nu + dim) / 2 * np.log1p(dk / nu) - logdet


def _maha(R, Y):
    L = linalg.cholesky(R, lower=True)
    y = linalg.solve_triangular(L, Y.T, lower=True)
    return np.einsum("ij,ij->j", y, y), np.log(np.diag(L)).sum()


def two_block_terms(Z, m):
    U1, U2, W = Z @ m.Q1, Z @ m.Q2, Z @ m.Qw
    U = np.hstack([U1, U2])
    n, r1 = U1.shape
    r2, dw = U2.shape[1], W.shape[1]
    logs = np.empty((n, m.K))
    d1 = np.empty((n, m.K))
    d2 = np.empty((n, m.K))
    for k in range(m.K):
        a, ld1 = _maha(m.scatter[k], U1 - m.means[k])
        b, ld2 = _maha(m.scatter2[k], U2 - m.intercept2[k] - U1 @ m.coefficient2[k].T)
        c, ldw = _maha(m.residual[k], W - m.intercept[k] - U @ m.coefficient[k].T)
        logs[:, k] = (np.log(m.weights[k]) + _t_log(a, m.df[k], r1, ld1) + _t_log(b, m.df2[k], r2, ld2)
                      - 0.5 * (dw * np.log(2 * np.pi) + c) - ldw)
        d1[:, k], d2[:, k] = a, b
    ll = logsumexp_rows(logs)
    return ll, np.exp(logs - ll[:, None]), d1, d2


def _wls(design, wgt, Y):
    Xg = design * wgt[:, None]
    b = linalg.solve(design.T @ Xg, Xg.T @ Y, assume_a="pos")
    return b, Y - design @ b


def fit_two_block_tmm(Z, gmm, nu0, Q1, Q2, Qw):
    t0 = time.time()
    U1, U2, W = Z @ Q1, Z @ Q2, Z @ Qw
    U = np.hstack([U1, U2])
    n, r1 = U1.shape
    r2, dw = U2.shape[1], W.shape[1]
    K = gmm.K
    s0 = (nu0 - 2) / nu0
    m = TwoBlockTMM(gmm.weights.copy(), gmm.means @ Q1, np.zeros((K, r1, r1)), np.full(K, nu0), np.zeros((K, r2)),
                    np.zeros((K, r2, r1)), np.zeros((K, r2, r2)), np.full(K, nu0), np.zeros((K, dw)),
                    np.zeros((K, dw, r1 + r2)), np.zeros((K, dw, dw)), Q1, Q2, Qw)
    Qu = np.hstack([Q1, Q2])
    for k in range(K):
        C, mu = gmm.scatter[k], gmm.means[k]
        S11, S21, S22 = Q1.T @ C @ Q1, Q2.T @ C @ Q1, Q2.T @ C @ Q2
        G2 = linalg.solve(S11, S21.T, assume_a="pos").T
        Suu, Swu, Sww = Qu.T @ C @ Qu, Qw.T @ C @ Qu, Qw.T @ C @ Qw
        G = linalg.solve(Suu, Swu.T, assume_a="pos").T
        m.scatter[k] = floor_cov(S11 * s0, COV_FLOOR * s0)
        m.coefficient2[k] = G2
        m.intercept2[k] = mu @ Q2 - G2 @ (mu @ Q1)
        m.scatter2[k] = floor_cov((S22 - G2 @ S11 @ G2.T) * s0, COV_FLOOR * s0)
        m.coefficient[k] = G
        m.intercept[k] = mu @ Qw - G @ (mu @ Qu)
        m.residual[k] = floor_cov(Sww - G @ Suu @ G.T, COV_FLOOR)
    ll, gamma, d1, d2 = two_block_terms(Z, m)
    hist = [ll.mean()]
    des1 = np.hstack([np.ones((n, 1)), U1])
    desu = np.hstack([np.ones((n, 1)), U])
    status, it = "maximum-iterations", 0
    for it in range(1, MAX_ITER + 1):
        Nk = gamma.sum(axis=0)
        if Nk.min() <= Z.shape[1] + 2:
            status = "collapsed-component"
            break
        x = TwoBlockTMM(Nk / n, *(np.copy(getattr(m, f)) for f in ("means", "scatter", "df", "intercept2", "coefficient2",
                                                                    "scatter2", "df2", "intercept", "coefficient", "residual")), Q1, Q2, Qw)
        for k in range(K):
            g = gamma[:, k]
            nu1, nu2 = m.df[k], m.df2[k]
            t1 = (nu1 + r1) / (nu1 + d1[:, k])
            w1 = g * t1
            mu = w1 @ U1 / w1.sum()
            D = U1 - mu
            x.means[k] = mu
            x.scatter[k] = floor_cov((D * w1[:, None]).T @ D / Nk[k], COV_FLOOR * (nu1 - 2) / nu1)
            x.df[k] = _solve_df(g, t1, d1[:, k], nu1, r1, Nk[k])
            t2 = (nu2 + r2) / (nu2 + d2[:, k])
            b, e = _wls(des1, g * t2, U2)
            x.intercept2[k], x.coefficient2[k] = b[0], b[1:].T
            x.scatter2[k] = floor_cov((e * (g * t2)[:, None]).T @ e / Nk[k], COV_FLOOR * (nu2 - 2) / nu2)
            x.df2[k] = _solve_df(g, t2, d2[:, k], nu2, r2, Nk[k])
            b, e = _wls(desu, g, W)
            x.intercept[k], x.coefficient[k] = b[0], b[1:].T
            x.residual[k] = floor_cov((e * g[:, None]).T @ e / Nk[k], COV_FLOOR)
        ll_n, gamma_n, d1n, d2n = two_block_terms(Z, x)
        value = ll_n.mean()
        gain = value - hist[-1]
        if gain < -DECREASE_TOL:
            status = "rejected-decreasing-update"
            break
        m, gamma, d1, d2 = x, gamma_n, d1n, d2n
        hist.append(value)
        if abs(gain) < LL_TOL:
            status = "converged"
            break
    Nk = gamma.sum(axis=0)
    m.audit = {"method": "two-strength block TMM", "K": K, "initial_df": nu0, "restart_direction": gmm.audit["restart_direction"],
               "status": status, "iterations": it, "mean_log_likelihood": float(hist[-1]),
               "effective_counts": Nk.tolist(), "seconds": time.time() - t0,
               "accepted": status == "converged" and Nk.min() > Z.shape[1] + 2}
    return m


def two_strength_directions(H1, H2, coords, rank_tolerance=1e-10):
    """Q1 spans the standardized wind-driven columns H1; Q2 spans the PV-driven columns H2 after removing
    their Q1 part; Qw is the orthonormal complement of [Q1, Q2]."""
    def orth(A):
        Ua, s, _ = linalg.svd(A, full_matrices=False)
        return Ua[:, : int((s > rank_tolerance * s[0]).sum())]
    A1 = (coords.basis.T @ H1) / coords.scale[:, None]
    A2 = (coords.basis.T @ H2) / coords.scale[:, None]
    Q1 = orth(A1)
    Q2 = orth(A2 - Q1 @ (Q1.T @ A2))
    Qu = np.hstack([Q1, Q2])
    Uf, _, _ = linalg.svd(Qu, full_matrices=True)
    Qw = Uf[:, Qu.shape[1]:]
    assert np.linalg.norm(Qu.T @ Qw) < 1e-8 and np.linalg.norm(Q1.T @ Q2) < 1e-8
    return Q1, Q2, Qw


def strength_directions(H, coords, rank_tolerance=1e-10):
    """Eq. (16): map the sensitivity columns H to standardized coordinates, keep the left singular vectors
    with singular value above rank_tolerance * the largest; Qw is the orthonormal complement."""
    Ht = (coords.basis.T @ H) / coords.scale[:, None]
    Uh, s, _ = linalg.svd(Ht, full_matrices=True)
    r = int((s > rank_tolerance * s[0]).sum())
    Qu, Qw = Uh[:, :r], Uh[:, r:]
    assert np.linalg.norm(Qu.T @ Qw) < 1e-8
    return Qu, Qw, s


# ---------------------------------------------------------------- order selection

def fit_one_order(Z, Zv, K, Qu=None, Qw=None, extra_block=None, verbose=False):
    """One component count. GMM from 3 initial splits, TMM from each accepted GMM with nu0 in {12, 200},
    block TMM from each accepted GMM with nu0 in {12, 200}; the accepted fit with the highest modeling mean
    log-likelihood is kept per method. extra_block: optional {method name: fit function(Z, gmm, nu0)}."""
    n, d = Z.shape
    rows, best = [], {}
    if n / K <= d + 2:
        return best, [{"method": "all", "K": K, "status": "skipped: N/K <= d + 2"}]
    gmms = [fit_gmm(Z, K, j) for j in RESTART_DIRECTIONS if j < d]
    cands = {"GMM": gmms, "TMM": [], "block TMM": []}
    for name in (extra_block or {}):
        cands[name] = []
    for g in gmms:
        if not g.audit["accepted"]:
            continue
        for nu0 in INITIAL_DF:
            cands["TMM"].append(fit_tmm(Z, g, nu0))
            if Qu is not None:
                cands["block TMM"].append(fit_block_tmm(Z, g, nu0, Qu, Qw))
            for name, fn in (extra_block or {}).items():
                cands[name].append(fn(Z, g, nu0))
    for method, lst in cands.items():
        for c in lst:
            a = c.audit
            rows.append({"method": method, "K": K, "restart_direction": a["restart_direction"], "initial_df": a.get("initial_df"),
                         "status": a["status"], "accepted": a["accepted"], "iterations": a["iterations"],
                         "modeling_mean_log_likelihood": a["mean_log_likelihood"], "min_effective_count": float(min(a["effective_counts"])),
                         "seconds": a["seconds"]})
            if verbose:
                print(method, K, a["restart_direction"], a.get("initial_df"), a["status"], a["iterations"], round(a["mean_log_likelihood"], 4), f"{a['seconds']:.1f}s", flush=True)
        acc = [c for c in lst if c.audit["accepted"]]
        if acc:
            m = max(acc, key=lambda c: c.audit["mean_log_likelihood"])
            m.audit["determination_mean_log_likelihood"] = float(log_density(Zv, m).mean())
            best[(method, K)] = m
    return best, rows


def _fit_one_order_args(args):
    return fit_one_order(*args)


def fit_orders(Z, Zv, Ks, Qu=None, Qw=None, workers=1, extra_block=None, verbose=False):
    """All candidate component counts (eq. (22) selection on the determination scenarios Zv).
    workers > 1 runs the component counts in separate processes (spawn); set OPENBLAS_NUM_THREADS=1 in the
    calling entry before numpy is imported."""
    order = sorted(Ks, reverse=True)
    args = [(Z, Zv, K, Qu, Qw, extra_block, verbose) for K in order]
    if workers > 1:
        import multiprocessing as mp
        with mp.get_context("spawn").Pool(min(workers, len(order))) as pool:
            parts = pool.map(_fit_one_order_args, args, chunksize=1)
    else:
        parts = [fit_one_order(*a) for a in args]
    best, rows = {}, []
    for (b, r), K in sorted(zip(parts, order), key=lambda t: t[1]):
        best.update(b)
        rows.extend(r)
    return best, select_orders(best), rows


def select_orders(best):
    selected = {}
    for method in sorted({mm for (mm, _) in best}):
        ks = [K for (mm, K) in best if mm == method]
        selected[method] = max(ks, key=lambda K: best[(method, K)].audit["determination_mean_log_likelihood"])
    return selected


# ---------------------------------------------------------------- storage

TWO_BLOCK_FIELDS = ("weights", "means", "scatter", "df", "intercept2", "coefficient2", "scatter2", "df2",
                    "intercept", "coefficient", "residual", "Q1", "Q2", "Qw")


def to_npz_dict(model, prefix=""):
    d = {f"{prefix}weights": model.weights, f"{prefix}means": model.means, f"{prefix}scatter": model.scatter,
         f"{prefix}df": model.df, f"{prefix}family": np.array(model.family)}
    if isinstance(model, BlockTMM):
        d.update({f"{prefix}intercept": model.intercept, f"{prefix}coefficient": model.coefficient,
                  f"{prefix}residual": model.residual, f"{prefix}Qu": model.Qu, f"{prefix}Qw": model.Qw})
    if isinstance(model, TwoBlockTMM):
        for k in TWO_BLOCK_FIELDS[4:]:
            d[f"{prefix}{k}"] = getattr(model, k)
    return d


def from_npz(z, prefix=""):
    fam = str(z[f"{prefix}family"])
    if fam == "two-block-tmm":
        return TwoBlockTMM(*(z[f"{prefix}{k}"] for k in TWO_BLOCK_FIELDS))
    if fam == "block-tmm":
        return BlockTMM(*(z[f"{prefix}{k}"] for k in ("weights", "means", "scatter", "df", "intercept", "coefficient", "residual", "Qu", "Qw")))
    return Mixture(z[f"{prefix}weights"], z[f"{prefix}means"], z[f"{prefix}scatter"], z[f"{prefix}df"], fam)
