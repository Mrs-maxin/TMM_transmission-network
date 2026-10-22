"""Distribution statistics, tail indices, binomial intervals and power-flow sensitivities for Monte Carlo results."""

import numpy as np
from scipy import stats

QUANTILES = [0.005, 0.01, 0.25, 0.75, 0.99, 0.995]


def column_statistics(X):
    """Mean, std, excess kurtosis, quantiles and tail spread index (q99.5 - q0.5) / IQR per column."""
    X = np.asarray(X, np.float64)
    mu = X.mean(axis=0)
    sd = X.std(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        z = (X - mu) / sd
        kurt = (z ** 4).mean(axis=0) - 3
    q = np.quantile(X, QUANTILES, axis=0)
    iqr = q[3] - q[2]
    with np.errstate(invalid="ignore", divide="ignore"):
        tail = (q[5] - q[0]) / iqr
    return {"mean": mu, "std": sd, "excess_kurtosis": kurt, "q005": q[0], "q01": q[1], "q25": q[2],
            "q75": q[3], "q99": q[4], "q995": q[5], "tail_index": tail}


def clopper_pearson(k, n, level=0.95):
    a = (1 - level) / 2
    k = np.asarray(k)
    lo = np.where(k > 0, stats.beta.ppf(a, k, n - k + 1), 0.0)
    hi = np.where(k < n, stats.beta.ppf(1 - a, k + 1, n - k), 1.0)
    return lo, hi


def probability_table(events, names):
    """events: bool (n, m). Returns list of dicts with count, probability and 95% interval."""
    n = events.shape[0]
    k = events.sum(axis=0)
    lo, hi = clopper_pearson(k, n)
    return [{"object": nm, "count": int(c), "n": int(n), "probability": float(c / n), "ci95_low": float(l), "ci95_high": float(h)}
            for nm, c, l, h in zip(names, k, lo, hi)]


def renewable_sensitivity(solver, V0, lam0, isPV, Qsp, P0):
    """Linear sensitivity of bus voltage magnitude, angle and branch-end P/Q to renewable active injections.

    Returns dict of arrays (n_objects, n_ren) per MW of renewable injection: dVm (p.u./MW), dVa (deg/MW),
    dPf, dQf, dPt, dQt (MW/MW or Mvar/MW)."""
    net = solver.net
    n = net.n
    m = n - 1
    _, I, _ = solver._mismatch(V0[None], np.array([lam0]), P0[None], Qsp[None], isPV[None])
    J = solver._jacobian(V0[None], I, isPV[None])[0]
    rows = solver.op.ren_bus
    rhs = np.zeros((2 * m + 1, len(rows)))
    rhs[rows, np.arange(len(rows))] = 1.0 / net.baseMVA
    dx = np.linalg.solve(J, rhs)
    dV = np.zeros((n, len(rows)), complex)
    dV[solver.nr] = dx[:m] + 1j * dx[m:2 * m]
    Vm = np.abs(V0)
    dVm = (V0.real[:, None] * dV.real + V0.imag[:, None] * dV.imag) / Vm[:, None]
    dVa = np.rad2deg((V0.real[:, None] * dV.imag - V0.imag[:, None] * dV.real) / Vm[:, None] ** 2)
    Vf, Vt = V0[net.f][:, None], V0[net.t][:, None]
    dVf, dVt = dV[net.f], dV[net.t]
    If = net.Yff[:, None] * Vf + net.Yft[:, None] * Vt
    It = net.Ytf[:, None] * Vf + net.Ytt[:, None] * Vt
    dIf = net.Yff[:, None] * dVf + net.Yft[:, None] * dVt
    dIt = net.Ytf[:, None] * dVf + net.Ytt[:, None] * dVt
    dSf = (dVf * np.conj(If) + Vf * np.conj(dIf)) * net.baseMVA
    dSt = (dVt * np.conj(It) + Vt * np.conj(dIt)) * net.baseMVA
    return {"dVm": dVm, "dVa": dVa, "dPf": dSf.real, "dQf": dSf.imag, "dPt": dSt.real, "dQt": dSt.imag,
            "dlam": dx[2 * m] * net.baseMVA}
