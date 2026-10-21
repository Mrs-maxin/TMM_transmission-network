"""Renewable forecast-error model at an operating point.

Historical hours whose capacity-normalised day-ahead forecast vector is closest to the operating point
(k nearest neighbours) supply the errors. Each site has a kernel-smoothed empirical marginal; values
outside the physical range forecast + error in [0, capacity] are set to the nearer bound (output at 0 or
at capacity). Dependence between sites is a t copula whose correlation matrix comes from Kendall's tau
and whose degrees of freedom is estimated by maximum likelihood; the shared chi-square variable of the t
copula is the common fluctuation strength. A Gaussian copula with the same marginals serves as the
counterpart without common fluctuation strength.

When wind and photovoltaic sites are both active, a grouped t copula is also fitted: each technology
group has its own fluctuation strength (independent chi-square variables, degrees of freedom nu_g) while
the Gaussian part keeps the full correlation matrix. Its likelihood is evaluated by Gauss-Laguerre
quadrature over the group strengths. The heavy-tailed model is the one with the lower AIC.
"""

from dataclasses import dataclass

import numpy as np
from scipy import optimize, special, stats

GRID_POINTS = 8001


@dataclass
class ErrorModel:
    sites: np.ndarray          # all mapped sites
    active: np.ndarray         # bool, sites with random error at this operating point
    forecast: np.ndarray       # MW (IEEE scale), all sites
    capacity: np.ndarray       # MW
    neighbours: np.ndarray     # (k,) historical hour indices
    distance: np.ndarray       # (k,) normalised forecast distance
    errors: np.ndarray         # (k, n_active) neighbour errors, MW (IEEE scale)
    bandwidth: np.ndarray      # (n_active,) MW
    grid: np.ndarray           # (n_active, GRID_POINTS) MW
    cdf: np.ndarray            # (n_active, GRID_POINTS)
    R_t: np.ndarray
    nu: float
    R_g: np.ndarray
    nu_profile: np.ndarray     # (m, 2): nu, log-likelihood
    nu_ci95: tuple
    group: np.ndarray = None   # (n_active,) technology group index, 0 wind, 1 photovoltaic
    nu_groups: np.ndarray = None
    heavy_kind: str = "t"      # "t" single strength or "grouped_t"

    def lower(self):
        return -self.forecast[self.active]

    def upper(self):
        return (self.capacity - self.forecast)[self.active]


def knn_hours(features, target, k):
    d = np.sqrt(((features - target[None, :]) ** 2).sum(axis=1))
    idx = np.argsort(d, kind="stable")[:k]
    return idx, d[idx]


def silverman_bandwidth(x):
    n = len(x)
    s = x.std(ddof=1)
    iqr = np.subtract(*np.quantile(x, [0.75, 0.25]))
    return 0.9 * min(s, iqr / 1.34) * n ** (-0.2)


def kde_cdf(x, h, lo, hi):
    grid = np.linspace(min(x.min(), lo) - 6 * h, max(x.max(), hi) + 6 * h, GRID_POINTS)
    cdf = np.zeros_like(grid)
    for chunk in np.array_split(x, max(1, len(x) // 200)):
        cdf += stats.norm.cdf((grid[:, None] - chunk[None, :]) / h).sum(axis=1)
    cdf /= len(x)
    return grid, cdf


def inverse_marginal(model, u):
    """Errors (MW) for uniform scores u of shape (n, n_active), set to the physical bounds when outside."""
    out = np.empty_like(u)
    lo, hi = model.lower(), model.upper()
    for j in range(u.shape[1]):
        c = model.cdf[j]
        keep = np.concatenate([[True], np.diff(c) > 0])
        out[:, j] = np.interp(u[:, j], c[keep], model.grid[j][keep])
        out[:, j] = np.clip(out[:, j], lo[j], hi[j])
    return out


def pseudo_obs(E):
    k = E.shape[0]
    return stats.rankdata(E, axis=0, method="average") / (k + 1)


def nearest_corr(R):
    w, V = np.linalg.eigh((R + R.T) / 2)
    w = np.maximum(w, 1e-6)
    R = (V * w) @ V.T
    d = np.sqrt(np.diag(R))
    return R / d[:, None] / d[None, :]


def kendall_corr(U):
    m = U.shape[1]
    R = np.eye(m)
    for i in range(m):
        for j in range(i + 1, m):
            tau = stats.kendalltau(U[:, i], U[:, j]).statistic
            R[i, j] = R[j, i] = np.sin(np.pi * tau / 2)
    return nearest_corr(R)


def t_copula_loglik(U, R, nu):
    x = stats.t.ppf(U, nu)
    m = U.shape[1]
    L = np.linalg.cholesky(R)
    logdet = 2 * np.log(np.diag(L)).sum()
    y = np.linalg.solve(L, x.T)
    q = (y ** 2).sum(axis=0)
    log_joint = (special.gammaln((nu + m) / 2) - special.gammaln(nu / 2) - m / 2 * np.log(nu * np.pi)
                 - 0.5 * logdet - (nu + m) / 2 * np.log1p(q / nu))
    log_marg = stats.t.logpdf(x, nu).sum(axis=1)
    return float((log_joint - log_marg).sum())


def gauss_copula_loglik(U, R):
    x = stats.norm.ppf(U)
    m = U.shape[1]
    L = np.linalg.cholesky(R)
    logdet = 2 * np.log(np.diag(L)).sum()
    y = np.linalg.solve(L, x.T)
    q = (y ** 2).sum(axis=0)
    return float((-0.5 * logdet - 0.5 * (q - (x ** 2).sum(axis=1))).sum())


def fit_nu(U, R, lo=2.05, hi=500.0):
    f = lambda lnu: -t_copula_loglik(U, R, np.exp(lnu))
    res = optimize.minimize_scalar(f, bounds=(np.log(lo), np.log(hi)), method="bounded", options={"xatol": 1e-6})
    nu = float(np.exp(res.x))
    grid = np.exp(np.linspace(np.log(lo), np.log(hi), 80))
    prof = np.array([[g, t_copula_loglik(U, R, g)] for g in grid])
    lmax = -res.fun
    inside = prof[prof[:, 1] >= lmax - 0.5 * stats.chi2.ppf(0.95, 1)]
    ci = (float(inside[:, 0].min()), float(inside[:, 0].max())) if len(inside) else (nu, nu)
    return nu, lmax, prof, ci


def grouped_t_loglik(U, R, group, nus, n_quad=48):
    """Log-likelihood of a t copula with independent fluctuation strengths per group (quadrature)."""
    groups = np.unique(group)
    x = np.empty_like(U)
    log_marg = np.zeros(U.shape[0])
    for g, nu in zip(groups, nus):
        sel = group == g
        x[:, sel] = stats.t.ppf(U[:, sel], nu)
        log_marg += stats.t.logpdf(x[:, sel], nu).sum(axis=1)
    m = U.shape[1]
    Rinv = np.linalg.inv(R)
    logdet = np.linalg.slogdet(R)[1]
    nodes, logw = [], []
    for nu in nus:
        t, w = special.roots_genlaguerre(n_quad, nu / 2 - 1)
        nodes.append(2 * t / nu)
        logw.append(np.log(w) - special.gammaln(nu / 2))
    sizes = [int((group == g).sum()) for g in groups]
    grids = np.meshgrid(*nodes, indexing="ij")
    lw = sum(np.meshgrid(*logw, indexing="ij"))
    A = np.stack([gr.ravel() for gr in grids], axis=1)
    lw = lw.ravel()
    scale = np.sqrt(A[:, np.searchsorted(groups, group)])
    out = np.empty(U.shape[0])
    for i0 in range(0, U.shape[0], 50):
        xs = x[i0:i0 + 50]
        y = xs[:, None, :] * scale[None, :, :]
        q = np.einsum("pqi,ij,pqj->pq", y, Rinv, y)
        logdens = (-0.5 * q - 0.5 * logdet - m / 2 * np.log(2 * np.pi)
                   + sum(sz / 2 * np.log(A[:, j]) for j, sz in enumerate(sizes))[None, :] + lw[None, :])
        out[i0:i0 + 50] = special.logsumexp(logdens, axis=1)
    return float((out - log_marg).sum())


def fit_grouped_nu(U, R, group, start):
    f = lambda lnus: -grouped_t_loglik(U, R, group, np.exp(np.clip(lnus, np.log(2.05), np.log(500))))
    res = optimize.minimize(f, np.log(start), method="Nelder-Mead", options={"xatol": 1e-3, "fatol": 1e-3, "maxiter": 200})
    return np.exp(np.clip(res.x, np.log(2.05), np.log(500))), -res.fun


def fit(sites, forecast_norm_hist, errors_hist_MW, forecast_op_MW, capacity_MW, active, k, group=None):
    """forecast_norm_hist: (T, n_sites) day-ahead / capacity; errors_hist_MW: (T, n_sites) IEEE-scale MW."""
    target = forecast_op_MW / capacity_MW
    idx, dist = knn_hours(forecast_norm_hist, target, k)
    E = errors_hist_MW[idx][:, active]
    if not np.isfinite(E).all():
        raise ValueError("non-finite neighbour errors")
    lo = -forecast_op_MW[active]
    hi = (capacity_MW - forecast_op_MW)[active]
    h = np.array([silverman_bandwidth(E[:, j]) for j in range(E.shape[1])])
    grids, cdfs = zip(*[kde_cdf(E[:, j], h[j], lo[j], hi[j]) for j in range(E.shape[1])])
    U = pseudo_obs(E)
    R_t = kendall_corr(U)
    nu, ll_t, prof, ci = fit_nu(U, R_t)
    R_g = nearest_corr(np.corrcoef(stats.norm.ppf(U), rowvar=False))
    ll_g = gauss_copula_loglik(U, R_g)
    m = E.shape[1]
    n_corr = m * (m - 1) // 2
    group = np.zeros(m, int) if group is None else np.asarray(group)
    ng = len(np.unique(group))
    info = {"loglik_t_copula": ll_t, "loglik_gauss_copula": ll_g, "lr_statistic_t_vs_gauss": 2 * (ll_t - ll_g),
            "aic_gauss": 2 * n_corr - 2 * ll_g, "aic_t": 2 * (n_corr + 1) - 2 * ll_t}
    heavy_kind, nus = "t", np.array([nu])
    if ng > 1:
        group_nu_separate = []
        for g in np.unique(group):
            sel = group == g
            group_nu_separate.append(fit_nu(U[:, sel], kendall_corr(U[:, sel]))[0])
        nus_g, ll_gt = fit_grouped_nu(U, R_t, group, np.array(group_nu_separate))
        info.update({"group_nu_from_group_subsets": group_nu_separate, "grouped_t_nu": nus_g.tolist(),
                     "loglik_grouped_t_copula": ll_gt, "aic_grouped_t": 2 * (n_corr + ng) - 2 * ll_gt,
                     "loglik_t_copula_quadrature_check": grouped_t_loglik(U, R_t, np.zeros(m, int), np.array([nu]))})
        if info["aic_grouped_t"] < info["aic_t"]:
            heavy_kind, nus = "grouped_t", nus_g
    info["heavy_tailed_model"] = heavy_kind
    model = ErrorModel(sites=np.asarray(sites, dtype=str), active=np.asarray(active), forecast=forecast_op_MW,
                       capacity=capacity_MW, neighbours=idx, distance=dist, errors=E, bandwidth=h,
                       grid=np.array(grids), cdf=np.array(cdfs), R_t=R_t, nu=nu, R_g=R_g, nu_profile=prof, nu_ci95=ci,
                       group=group, nu_groups=nus if heavy_kind == "grouped_t" else np.array([nu] * ng),
                       heavy_kind=heavy_kind)
    return model, info


def sample_uniform(model, n, rng, kind):
    """kind: "t" single-strength t copula, "grouped_t", "heavy" (the selected heavy-tailed model) or "gauss".
    Returns uniform scores and the fluctuation strength variables w (n, n_groups); w = 1 for "gauss"."""
    if kind == "heavy":
        kind = model.heavy_kind
    m = int(model.active.sum())
    group = model.group if model.group is not None else np.zeros(m, int)
    groups = np.unique(group)
    R = model.R_g if kind == "gauss" else model.R_t
    L = np.linalg.cholesky(R)
    z = rng.standard_normal((n, m)) @ L.T
    if kind == "gauss":
        return stats.norm.cdf(z), np.ones((n, len(groups)))
    if kind == "t":
        w1 = rng.chisquare(model.nu, n) / model.nu
        return stats.t.cdf(z / np.sqrt(w1)[:, None], model.nu), np.repeat(w1[:, None], len(groups), axis=1)
    u = np.empty_like(z)
    w = np.empty((n, len(groups)))
    for j, (g, nu) in enumerate(zip(groups, model.nu_groups)):
        sel = group == g
        w[:, j] = rng.chisquare(nu, n) / nu
        u[:, sel] = stats.t.cdf(z[:, sel] / np.sqrt(w[:, j])[:, None], nu)
    return u, w


def sample_errors(model, n, rng, kind):
    """Errors in MW for all sites (zeros at inactive sites) and the common strength variable."""
    u, w = sample_uniform(model, n, rng, kind)
    u = np.clip(u, 1e-12, 1 - 1e-12)
    e = np.zeros((n, len(model.sites)))
    e[:, model.active] = inverse_marginal(model, u)
    return e, w


def to_npz_dict(model):
    d = {k: getattr(model, k) for k in model.__dataclass_fields__ if k not in ("nu_ci95", "heavy_kind")}
    return d | {"nu_ci95": np.array(model.nu_ci95), "heavy_kind": np.array(model.heavy_kind)}


def from_npz(z):
    d = {k: z[k] for k in z.files}
    d["nu"] = float(d["nu"])
    d["nu_ci95"] = tuple(d["nu_ci95"])
    d["heavy_kind"] = str(d["heavy_kind"])
    return ErrorModel(**d)


def dependence_statistics(E):
    """Excess kurtosis per column, pairwise Pearson correlation of errors and of squared standardised errors."""
    Z = (E - E.mean(axis=0)) / E.std(axis=0)
    kurt = (Z ** 4).mean(axis=0) - 3
    C = np.corrcoef(E, rowvar=False)
    C2 = np.corrcoef(Z ** 2, rowvar=False)
    return kurt, C, C2
