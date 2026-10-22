"""Batched rectangular-coordinate Newton-Raphson AC power flow with distributed active-power balancing
and unit reactive-power limits.

Unknowns per scenario: real and imaginary voltage parts of every non-reference bus plus the scalar
system imbalance lam (p.u.). Unit g produces Pg0_g + alpha_g * lam. Equations: active power at every
bus (reference bus included), reactive power at PQ buses, squared voltage magnitude at PV buses.

Reactive limits: start with all non-reference units regulating; after convergence every regulating
unit whose reactive output leaves [Qmin - eps_q, Qmax + eps_q] is fixed at the violated limit and its
bus becomes PQ. All violating units switch together, a fixed unit never returns to regulation, and the
procedure repeats until no new violation appears. The reference unit has no reactive limit.
"""

from dataclasses import dataclass

import numpy as np

FAIL_OK = 0
FAIL_MAXITER = 1
FAIL_VOLTAGE_RANGE = 2
FAIL_NONFINITE = 3
FAIL_MAXROUNDS = 4
FAIL_SINGULAR = 5
FAIL_REASONS = {
    FAIL_OK: "converged",
    FAIL_MAXITER: "Newton iterations exceeded limit",
    FAIL_VOLTAGE_RANGE: "voltage magnitude left numerical range 0.5-1.5 p.u.",
    FAIL_NONFINITE: "non-finite mismatch",
    FAIL_MAXROUNDS: "reactive-limit switching rounds exceeded limit",
    FAIL_SINGULAR: "singular Jacobian",
}


@dataclass
class OperatingPoint:
    Pg0: np.ndarray      # unit base active output, p.u.
    Vset: np.ndarray     # unit voltage setpoint, p.u.
    alpha: np.ndarray    # unit participation factors, sum 1
    Pd: np.ndarray       # bus active load at the operating point, p.u.
    Qd: np.ndarray       # bus reactive load, p.u.
    ren_bus: np.ndarray  # host bus index of each renewable injection
    Pren: np.ndarray     # renewable day-ahead forecast, p.u.
    ren_cap: np.ndarray  # renewable installed capacity, p.u.


class BatchSolver:
    def __init__(self, net, op, tol=1e-10, max_iter=30, max_rounds=30, eps_q=5e-8):
        self.net = net
        self.op = op
        self.tol = tol
        self.max_iter = max_iter
        self.max_rounds = max_rounds
        self.eps_q = eps_q
        n = net.n
        self.n = n
        self.Y = net.Y
        self.Yc = np.conj(net.Y)
        self.nr = np.setdiff1d(np.arange(n), [net.ref])
        self.alpha_bus = np.zeros(n)
        np.add.at(self.alpha_bus, net.gen_bus, op.alpha)
        self.Pg0_bus = np.zeros(n)
        np.add.at(self.Pg0_bus, net.gen_bus, op.Pg0)
        self.reg_gen = np.flatnonzero(net.gen_bus != net.ref)
        self.reg_bus = net.gen_bus[self.reg_gen]
        self.ref_gen = np.flatnonzero(net.gen_bus == net.ref)
        self.Vref = op.Vset[self.ref_gen[0]] if len(self.ref_gen) else 1.0
        self.Vset_bus = np.ones(n)
        self.Vset_bus[net.gen_bus] = op.Vset
        self.Qmax_bus = np.full(n, np.inf)
        self.Qmin_bus = np.full(n, -np.inf)
        self.Qmax_bus[self.reg_bus] = net.Qmax[self.reg_gen]
        self.Qmin_bus[self.reg_bus] = net.Qmin[self.reg_gen]
        self.pos_in_nr = -np.ones(n, int)
        self.pos_in_nr[self.nr] = np.arange(n - 1)

    def injections(self, Pd, Qd, Pren):
        """Specified bus injections without the shared imbalance term."""
        P = self.Pg0_bus[None, :] - Pd
        np.add.at(P, (slice(None), self.op.ren_bus), Pren)
        Q = -Qd.copy()
        return P, Q

    def _mismatch(self, V, lam, P0, Qsp, isPV):
        I = V @ self.Y.T
        S = V * np.conj(I)
        mP = S.real - (P0 + lam[:, None] * self.alpha_bus[None, :])
        Vm2 = V.real ** 2 + V.imag ** 2
        mQ = np.where(isPV, Vm2 - self.Vset_bus[None, :] ** 2, S.imag - Qsp)
        F = np.concatenate([mP, mQ[:, self.nr]], axis=1)
        return F, I, S

    def _jacobian(self, V, I, isPV):
        B, n = V.shape
        nr = self.nr
        VY = V[:, :, None] * self.Yc[None, :, :]
        dSe = VY.copy()
        idx = np.arange(n)
        dSe[:, idx, idx] += np.conj(I)
        dSf = -1j * VY
        dSf[:, idx, idx] += 1j * np.conj(I)
        m = n - 1
        J = np.zeros((B, n + m, 2 * m + 1))
        J[:, :n, :m] = dSe.real[:, :, nr]
        J[:, :n, m:2 * m] = dSf.real[:, :, nr]
        J[:, :n, 2 * m] = -self.alpha_bus[None, :]
        Jq_e = dSe.imag[:, nr][:, :, nr]
        Jq_f = dSf.imag[:, nr][:, :, nr]
        pv = isPV[:, nr]
        k = np.arange(m)
        Jv_e = np.zeros((B, m, m))
        Jv_f = np.zeros((B, m, m))
        Jv_e[:, k, k] = 2 * V.real[:, nr]
        Jv_f[:, k, k] = 2 * V.imag[:, nr]
        J[:, n:, :m] = np.where(pv[:, :, None], Jv_e, Jq_e)
        J[:, n:, m:2 * m] = np.where(pv[:, :, None], Jv_f, Jq_f)
        return J

    def newton(self, V, lam, P0, Qsp, isPV):
        """Newton-Raphson for a batch with fixed bus types. Returns V, lam, iterations, fail code."""
        B = V.shape[0]
        m = self.n - 1
        it = np.zeros(B, int)
        fail = np.full(B, -1)
        active = np.arange(B)
        for k in range(self.max_iter + 1):
            F, I, _ = self._mismatch(V[active], lam[active], P0[active], Qsp[active], isPV[active])
            err = np.abs(F).max(axis=1)
            Vm = np.abs(V[active])
            bad_num = ~np.isfinite(err)
            bad_rng = (Vm.min(axis=1) < 0.5) | (Vm.max(axis=1) > 1.5)
            done = err < self.tol
            fail[active[done & ~bad_rng]] = FAIL_OK
            fail[active[bad_num]] = FAIL_NONFINITE
            fail[active[bad_rng & ~bad_num]] = FAIL_VOLTAGE_RANGE
            keep = ~(done | bad_num | bad_rng)
            if k == self.max_iter:
                fail[active[keep]] = FAIL_MAXITER
                break
            active, F, I = active[keep], F[keep], I[keep]
            if len(active) == 0:
                break
            J = self._jacobian(V[active], I, isPV[active])
            try:
                dx = np.linalg.solve(J, -F[:, :, None])[:, :, 0]
            except np.linalg.LinAlgError:
                dx = np.empty_like(F)
                for j in range(len(active)):
                    try:
                        dx[j] = np.linalg.solve(J[j], -F[j])
                    except np.linalg.LinAlgError:
                        dx[j] = np.nan
                        fail[active[j]] = FAIL_SINGULAR
            ok = np.isfinite(dx).all(axis=1)
            a = active[ok]
            V[a[:, None], self.nr[None, :]] += dx[ok, :m] + 1j * dx[ok, m:2 * m]
            lam[a] += dx[ok, 2 * m]
            it[active] += 1
            active = a
        return V, lam, it, fail

    def solve(self, Pd, Qd, Pren, V0=None):
        """Solve a batch of scenarios. Pd, Qd: (B, n) p.u.; Pren: (B, n_ren) p.u."""
        net = self.net
        B = Pd.shape[0]
        P0, Q0 = self.injections(Pd, Qd, Pren)
        if V0 is None:
            V = np.tile(self.Vset_bus.astype(complex), (B, 1))
            V[:, net.ref] = self.Vref
        else:
            V = V0.copy()
        lam = np.zeros(B)
        mode = np.zeros((B, len(self.reg_gen)), int)
        isPV = np.zeros((B, self.n), bool)
        isPV[:, self.reg_bus] = True
        Qsp = Q0.copy()
        iters = np.zeros(B, int)
        rounds = np.zeros(B, int)
        fail = np.full(B, -1)
        active = np.arange(B)
        for r in range(self.max_rounds + 1):
            Va, la, it, fl = self.newton(V[active], lam[active], P0[active], Qsp[active], isPV[active])
            V[active], lam[active] = Va, la
            iters[active] += it
            rounds[active] += 1
            bad = fl != FAIL_OK
            fail[active[bad]] = fl[bad]
            active = active[~bad]
            if len(active) == 0:
                break
            S = V[active] * np.conj(V[active] @ self.Y.T)
            Qg = S.imag[:, self.reg_bus] + Qd[active][:, self.reg_bus]
            reg = mode[active] == 0
            up = reg & (Qg > self.Qmax_bus[self.reg_bus] + self.eps_q)
            lo = reg & (Qg < self.Qmin_bus[self.reg_bus] - self.eps_q)
            new = up | lo
            has_new = new.any(axis=1)
            fail[active[~has_new]] = FAIL_OK
            if r == self.max_rounds:
                fail[active[has_new]] = FAIL_MAXROUNDS
                break
            sel = active[has_new]
            m_sel = mode[sel]
            m_sel[up[has_new]] = 1
            m_sel[lo[has_new]] = -1
            mode[sel] = m_sel
            rows, cols = np.nonzero(new[has_new])
            gb = self.reg_bus[cols]
            lim = np.where(up[has_new][rows, cols], self.Qmax_bus[gb], self.Qmin_bus[gb])
            Qsp[sel[rows], gb] = lim - Qd[sel[rows], gb]
            isPV[sel[rows], gb] = False
            active = sel
        return self.results(V, lam, mode, Pd, Qd, Pren, iters, rounds, fail)

    def results(self, V, lam, mode, Pd, Qd, Pren, iters, rounds, fail):
        net = self.net
        base = net.baseMVA
        S = V * np.conj(V @ self.Y.T)
        Pg = (self.op.Pg0[None, :] + lam[:, None] * self.op.alpha[None, :]) * base
        Qg = (S.imag[:, net.gen_bus] + Qd[:, net.gen_bus]) * base
        Sf, St = net.branch_flows(V)
        Pres = S.real - (self.injections(Pd, Qd, Pren)[0] + lam[:, None] * self.alpha_bus[None, :])
        return {
            "V": V,
            "lam_MW": lam * base,
            "Pg_MW": Pg,
            "Qg_Mvar": Qg,
            "mode": mode,
            "Sf_MVA": Sf,
            "St_MVA": St,
            "iterations": iters,
            "rounds": rounds,
            "fail": fail,
            "max_P_residual_pu": np.abs(Pres).max(axis=1),
        }
