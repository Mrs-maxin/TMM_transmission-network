"""Power-flow state vector x = [e of non-reference buses, f of non-reference buses, system imbalance lam (p.u.)],
the distributed-slack Jacobian used for the common-strength directions (derivation.md eq. (15)), and the
electrical quantities rebuilt from x by eqs. (23)-(24)."""

from pathlib import Path

import numpy as np

from . import montecarlo as mcarlo
from . import powerflow as pf


def scenario_state(paths, net):
    """Return (x, data) for the stage D chunks in paths. Non-converged scenarios raise: phase-1 sets have none."""
    d = mcarlo.load_set(sorted(Path(p) for p in paths))
    if not d["converged"].all():
        raise RuntimeError(f"{int((~d['converged']).sum())} non-converged scenarios in {paths[0]}")
    nr = np.setdiff1d(np.arange(net.n), [net.ref])
    V = d["V"].astype(np.complex128)
    x = np.hstack([V.real[:, nr], V.imag[:, nr], d["imbalance_MW"][:, None] / net.baseMVA])
    return x, d


def most_frequent_mode(mode):
    u, cnt = np.unique(mode, axis=0, return_counts=True)
    j = int(np.argmax(cnt))
    return u[j].astype(int), float(cnt[j] / len(mode))


def jacobian_at_mode(net, op, mode_vec, Pd, Qd, Pren):
    """Solve the power flow with the unit reactive-limit mode fixed (0 regulating, +1 at Qmax, -1 at Qmin) at the
    given injections and return the Jacobian of eq. (7) with the imbalance column, plus V and lam."""
    s = pf.BatchSolver(net, op)
    P0, Q0 = s.injections(Pd[None], Qd[None], Pren[None])
    isPV = np.zeros((1, net.n), bool)
    isPV[0, s.reg_bus] = mode_vec == 0
    Qsp = Q0.copy()
    lim = np.where(mode_vec == 1, s.Qmax_bus[s.reg_bus], s.Qmin_bus[s.reg_bus])
    sel = mode_vec != 0
    Qsp[0, s.reg_bus[sel]] = lim[sel] - Qd[s.reg_bus[sel]]
    V = s.Vset_bus.astype(complex)[None].copy()
    V[0, net.ref] = s.Vref
    V, lam, it, fl = s.newton(V, np.zeros(1), P0, Qsp, isPV)
    if fl[0] != pf.FAIL_OK:
        raise RuntimeError(f"Jacobian operating point did not converge: {pf.FAIL_REASONS[int(fl[0])]}")
    _, I, _ = s._mismatch(V, lam, P0, Qsp, isPV)
    J = s._jacobian(V, I, isPV)[0]
    return J, V[0], float(lam[0]), int(it[0])


def injection_columns(J, rows):
    """Columns of J^-1 for unit injections at the given equation rows (active-power rows are bus indices)."""
    E = np.zeros((J.shape[0], len(rows)))
    E[np.asarray(rows), np.arange(len(rows))] = 1.0
    return np.linalg.solve(J, E)


def renewable_rows(net, op, mapping):
    """Active-power equation rows of the renewable host buses, read from the stage-B renewable mapping
    (rows in the order of the renewable injections), checked against the operating point."""
    rows = net.idx(mapping.host_bus.to_numpy())
    if not np.array_equal(rows, op.ren_bus):
        raise ValueError("renewable mapping and operating point host buses differ")
    return rows, mapping.kind.to_numpy()


def strength_setup(net, op, mapping, modeling, coords):
    """Eqs. (15)-(16) for the OP2 input: one common strength over all renewable active injections (main block
    model), and wind / PV columns separately for the two-strength variant."""
    from . import mixture as mx
    mode, share = most_frequent_mode(modeling["mode"])
    k = modeling["load_factor"].mean(axis=0)
    Pren = modeling["Pren_MW"].mean(axis=0) / net.baseMVA
    J, V0, lam0, it = jacobian_at_mode(net, op, mode, op.Pd * k, op.Qd * k, Pren)
    rows, kind = renewable_rows(net, op, mapping)
    H = injection_columns(J, rows)
    Qu, Qw, sv = mx.strength_directions(H, coords)
    Q1, Q2, Qw2 = mx.two_strength_directions(H[:, kind == "wind"], H[:, kind == "pv"], coords)
    info = {"mode_share": share, "mode_units_at_limit": int((mode != 0).sum()), "newton_iterations": it,
            "imbalance_at_point_MW": lam0 * net.baseMVA, "jacobian_condition": float(np.linalg.cond(J)),
            "r": int(Qu.shape[1]), "singular_values": sv.tolist(), "r_wind": int(Q1.shape[1]), "r_pv": int(Q2.shape[1]),
            "host_buses": net.bus_id[rows].tolist(), "host_kind": kind.tolist()}
    return {"Qu": Qu, "Qw": Qw, "Q1": Q1, "Q2": Q2, "Qw2": Qw2, "H": H, "J": J, "mode": mode, "V0": V0}, info


class OutputMap:
    """Electrical quantities of state points x (n, 2m+1), eqs. (23)-(24) and Pg = Pg0 + alpha * lam."""

    def __init__(self, net, op):
        self.net = net
        self.op = op
        self.nr = np.setdiff1d(np.arange(net.n), [net.ref])
        s = pf.BatchSolver(net, op)
        self.Vref = s.Vref
        self.share = net.Pmax > 0

    def voltages(self, x):
        m = len(self.nr)
        V = np.empty((x.shape[0], self.net.n), complex)
        V[:, self.nr] = x[:, :m] + 1j * x[:, m:2 * m]
        V[:, self.net.ref] = self.Vref
        return V

    def evaluate(self, x):
        net = self.net
        V = self.voltages(x)
        Sf, St = net.branch_flows(V)
        lam_MW = x[:, -1] * net.baseMVA
        Pg = self.op.Pg0[None] * net.baseMVA + lam_MW[:, None] * self.op.alpha[None]
        return {"Vm": np.abs(V[:, self.nr]), "Va": np.rad2deg(np.angle(V[:, self.nr])),
                "Pf": Sf.real, "Qf": Sf.imag, "Pt": St.real, "Qt": St.imag,
                "loading": np.maximum(np.abs(Sf), np.abs(St)) / net.rateA[None],
                "lam": lam_MW[:, None], "Pg": Pg[:, self.share]}

    def from_stage_d(self, d):
        """Same quantities from stored stage D power-flow results."""
        net = self.net
        Vm = np.abs(d["V"][:, self.nr]).astype(np.float64)
        Sf = d["Pf_MW"].astype(np.float64) + 1j * d["Qf_Mvar"].astype(np.float64)
        St = d["Pt_MW"].astype(np.float64) + 1j * d["Qt_Mvar"].astype(np.float64)
        return {"Vm": Vm, "Va": np.rad2deg(np.angle(d["V"][:, self.nr])).astype(np.float64),
                "Pf": Sf.real, "Qf": Sf.imag, "Pt": St.real, "Qt": St.imag,
                "loading": np.maximum(np.abs(Sf), np.abs(St)) / net.rateA[None],
                "lam": d["imbalance_MW"].astype(np.float64)[:, None],
                "Pg": d["Pg_MW"].astype(np.float64)[:, self.share]}
