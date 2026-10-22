"""Expansion-point batch power flow (01_IDEA/derivation.md Sections 1-2) with the distributed active-power
imbalance as an additional unknown.

Unknowns x = [e_nr, f_nr, lam]. The equations P_i(V) - alpha_i * lam (all buses), Q_i(V) (PQ buses) and
|V_i|^2 (PV buses) are quadratic in x, so F(x0 + d) = F(x0) + J0 d + h(d) with
h(d) = [F(x0 + d) + F(x0 - d)] / 2 - F(x0) holds exactly.

For one unit-mode group: the first expansion point is the scenario whose specified quantities are closest
to the group mean; it is solved by Newton-Raphson. Each scenario starts from the second-order value
eq. (8), then takes at most 10 fixed-Jacobian corrections eq. (9) with step 1, 1/2, ..., 1/2^8 (largest
step that lowers the residual 2-norm). Convergence follows eq. (10). Unconverged scenarios supply the next
expansion point (largest ||J0^-1||_2 ||r||_2), at most 12 expansion points per group; the remaining
scenarios are solved point by point by Newton-Raphson. Unit reactive limits follow Section 2.3.
"""

import numpy as np
from scipy import linalg

from . import powerflow as pf

EPS_Y = 1e-8
EPS_X = 1e-7
MAX_CORR = 10
MAX_EXPANSION = 12
STEPS = 0.5 ** np.arange(9)


class BatchExpansionSolver:
    def __init__(self, net, op, eps_q=5e-8):
        self.pf = pf.BatchSolver(net, op, eps_q=eps_q)
        self.net = net
        self.stats = {"expansion_points": 0, "newton_point_solutions": 0, "mode_groups": 0, "group_sizes": []}

    def _x_to_V(self, x):
        s = self.pf
        m = s.n - 1
        V = np.empty((x.shape[0], s.n), complex)
        V[:, s.nr] = x[:, :m] + 1j * x[:, m:2 * m]
        V[:, self.net.ref] = s.Vref
        return V, x[:, 2 * m]

    def _V_to_x(self, V, lam):
        s = self.pf
        return np.concatenate([V[:, s.nr].real, V[:, s.nr].imag, lam[:, None]], axis=1)

    def _mis(self, x, P0, Qsp, isPV):
        V, lam = self._x_to_V(x)
        F, _, _ = self.pf._mismatch(V, lam, P0, Qsp, isPV)
        return F

    def _converged(self, r, norm_inv, x):
        V, _ = self._x_to_V(x)
        Vm = np.abs(V)
        ok_rng = (Vm.min(axis=1) >= 0.5) & (Vm.max(axis=1) <= 1.5)
        return (np.abs(r).max(axis=1) <= EPS_Y) & (norm_inv * np.linalg.norm(r, axis=1) <= EPS_X) & ok_rng

    def _newton_points(self, P0, Qsp, isPV):
        s = self.pf
        B = P0.shape[0]
        V = np.tile(s.Vset_bus.astype(complex), (B, 1))
        V[:, self.net.ref] = s.Vref
        V, lam, it, fl = s.newton(V, np.zeros(B), P0, Qsp, isPV)
        self.stats["newton_point_solutions"] += B
        return self._V_to_x(V, lam), fl

    def solve_group(self, P0, Qsp, isPV):
        """Solve scenarios that share one bus-type pattern. Returns x, converged flag, method code
        (0 batch, 1 expansion point itself, 2 point-wise Newton), number of expansion points."""
        B = P0.shape[0]
        y = np.concatenate([P0, np.where(isPV, 0.0, Qsp)[:, self.pf.nr]], axis=1)
        x = np.full((B, 2 * (self.pf.n - 1) + 1), np.nan)
        done = np.zeros(B, bool)
        method = np.full(B, -1)
        pending = np.arange(B)
        first = int(np.argmin(np.linalg.norm(y - y.mean(axis=0), axis=1)))
        n_exp = 0
        while len(pending) and n_exp < MAX_EXPANSION:
            k = first if n_exp == 0 else next_point
            x0, fl = self._newton_points(P0[k:k + 1], Qsp[k:k + 1], isPV[k:k + 1])
            n_exp += 1
            if fl[0] != pf.FAIL_OK:
                break
            x[k], done[k], method[k] = x0[0], True, 1
            pending = pending[pending != k]
            if not len(pending):
                break
            V0, lam0 = self._x_to_V(x0)
            _, I0, _ = self.pf._mismatch(V0, lam0, P0[k:k + 1], Qsp[k:k + 1], isPV[k:k + 1])
            J0 = self.pf._jacobian(V0, I0, isPV[k:k + 1])[0]
            lu = linalg.lu_factor(J0)
            norm_inv = 1.0 / linalg.svdvals(J0).min()
            Pp, Qp, Vp = P0[pending], Qsp[pending], isPV[pending]
            X0 = np.repeat(x0, len(pending), axis=0)
            m0 = self._mis(X0, Pp, Qp, Vp)
            d1 = linalg.lu_solve(lu, -m0.T).T
            h = 0.5 * (self._mis(X0 + d1, Pp, Qp, Vp) + self._mis(X0 - d1, Pp, Qp, Vp)) - m0
            xs = X0 + d1 - linalg.lu_solve(lu, h.T).T
            r = -self._mis(xs, Pp, Qp, Vp)
            conv = self._converged(r, norm_inv, xs)
            for _ in range(MAX_CORR):
                act = np.flatnonzero(~conv & np.isfinite(r).all(axis=1))
                if not len(act):
                    break
                dx = linalg.lu_solve(lu, r[act].T).T
                nr0 = np.linalg.norm(r[act], axis=1)
                best_x = xs[act].copy()
                best_r = r[act].copy()
                accepted = np.zeros(len(act), bool)
                for a in STEPS:
                    cand = ~accepted
                    if not cand.any():
                        break
                    xt = xs[act[cand]] + a * dx[cand]
                    rt = -self._mis(xt, Pp[act[cand]], Qp[act[cand]], Vp[act[cand]])
                    better = np.linalg.norm(rt, axis=1) < nr0[cand]
                    ci = np.flatnonzero(cand)[better]
                    best_x[ci], best_r[ci] = xt[better], rt[better]
                    accepted[ci] = True
                xs[act], r[act] = best_x, best_r
                conv[act] = self._converged(r[act], norm_inv, xs[act])
                if not accepted.any():
                    break
            ok = conv
            x[pending[ok]], done[pending[ok]], method[pending[ok]] = xs[ok], True, 0
            rest = ~ok
            if not rest.any():
                pending = pending[:0]
                break
            score = norm_inv * np.linalg.norm(np.nan_to_num(r[rest], nan=np.inf), axis=1)
            next_point = int(pending[rest][np.argmax(score)])
            pending = pending[rest]
        if len(pending):
            xn, fl = self._newton_points(P0[pending], Qsp[pending], isPV[pending])
            good = fl == pf.FAIL_OK
            x[pending[good]], done[pending[good]], method[pending[good]] = xn[good], True, 2
        self.stats["expansion_points"] += n_exp
        return x, done, method, n_exp

    def solve(self, Pd, Qd, Pren, max_rounds=30):
        s = self.pf
        net = self.net
        B = Pd.shape[0]
        P0, Q0 = s.injections(Pd, Qd, Pren)
        mode = np.zeros((B, len(s.reg_gen)), int)
        isPV = np.zeros((B, s.n), bool)
        isPV[:, s.reg_bus] = True
        Qsp = Q0.copy()
        x = np.full((B, 2 * (s.n - 1) + 1), np.nan)
        method = np.full(B, -1)
        fail = np.full(B, -1)
        active = np.arange(B)
        for rnd in range(max_rounds + 1):
            keys, inv = np.unique(mode[active], axis=0, return_inverse=True)
            inv = inv.ravel()
            for g in range(len(keys)):
                idx = active[inv == g]
                self.stats["mode_groups"] += 1
                self.stats["group_sizes"].append(len(idx))
                xg, ok, mg, _ = self.solve_group(P0[idx], Qsp[idx], isPV[idx])
                x[idx], method[idx] = xg, mg
                fail[idx[~ok]] = pf.FAIL_MAXITER
            active = active[fail[active] == -1]
            if not len(active):
                break
            V, lam = self._x_to_V(x[active])
            S = V * np.conj(V @ s.Y.T)
            Qg = S.imag[:, s.reg_bus] + Qd[active][:, s.reg_bus]
            reg = mode[active] == 0
            up = reg & (Qg > s.Qmax_bus[s.reg_bus] + s.eps_q)
            lo = reg & (Qg < s.Qmin_bus[s.reg_bus] - s.eps_q)
            new = up | lo
            has = new.any(axis=1)
            fail[active[~has]] = pf.FAIL_OK
            if rnd == max_rounds:
                fail[active[has]] = pf.FAIL_MAXROUNDS
                break
            sel = active[has]
            ms = mode[sel]
            ms[up[has]] = 1
            ms[lo[has]] = -1
            mode[sel] = ms
            rows, cols = np.nonzero(new[has])
            gb = s.reg_bus[cols]
            lim = np.where(up[has][rows, cols], s.Qmax_bus[gb], s.Qmin_bus[gb])
            Qsp[sel[rows], gb] = lim - Qd[sel[rows], gb]
            isPV[sel[rows], gb] = False
            active = sel
        V, lam = self._x_to_V(np.nan_to_num(x))
        return {"V": V, "lam_MW": lam * net.baseMVA, "mode": mode, "fail": fail, "method": method}
