"""Operating-limit violation events of the OP2 power-flow state and their probabilities from AC Monte Carlo
scenarios (counting, Clopper-Pearson interval) and from stratified Sobol points of a fitted state model."""

import numpy as np
from scipy import stats

from . import propagation as pg

CURVE_LEVELS = np.round(np.arange(1.00, 1.3001, 0.01), 2)
EVENTS = ["a_branch_100_103_above_1.05", "b_any_unit_below_Pmin", "c_any_unit_above_Pmax", "d_branch_77_82_above_1.00",
          "e_bus83_above_Vmax", "f_any_branch_above_1.05", "f_any_branch_above_1.10", "g_any_bus_above_Vmax_plus_0.0025"]


def thresholds(net, op):
    """Event thresholds from the network and operating-point data. Units with alpha > 0 move with the imbalance;
    unit i is below Pmin iff lam < (Pmin_i - Pg0_i)/alpha_i, so 'any unit below Pmin' iff lam < max_i of that ratio."""
    br = {f"{net.bus_id[f]}-{net.bus_id[t]}": j for j, (f, t) in enumerate(zip(net.f, net.t))}
    nr = np.setdiff1d(np.arange(net.n), [net.ref])
    part = op.alpha > 0
    low = (net.Pmin[part] - op.Pg0[part]) / op.alpha[part] * net.baseMVA
    high = (net.Pmax[part] - op.Pg0[part]) / op.alpha[part] * net.baseMVA
    b83 = int(net.idx([83])[0])
    return {"branch_100_103": br["100-103"], "branch_77_82": br["77-82"], "bus83_col": int(np.flatnonzero(nr == b83)[0]),
            "Vmax_bus83": float(net.Vmax[b83]), "Vmax_nonref": net.Vmax[nr].copy(),
            "lam_low_MW": float(low.max()), "lam_high_MW": float(high.min()),
            "lam_low_unit_bus": int(net.bus_id[net.gen_bus[part][np.argmax(low)]]),
            "lam_high_unit_bus": int(net.bus_id[net.gen_bus[part][np.argmin(high)]]),
            "Pmin_MW": net.Pmin[net.Pmax > 0] * net.baseMVA, "Pmax_MW": net.Pmax[net.Pmax > 0] * net.baseMVA,
            "rateA_100_103_MVA": float(net.rateA[br["100-103"]]), "rateA_77_82_MVA": float(net.rateA[br["77-82"]])}


def indicators(q, thr):
    """Boolean (n, events) array in EVENTS order and (n, levels) array of the branch 100-103 curve, from the
    quantities of OutputMap.evaluate / from_stage_d."""
    L, lam, Vm = q["loading"], q["lam"][:, 0], q["Vm"]
    l100 = L[:, thr["branch_100_103"]]
    Lmax = L.max(axis=1)
    ev = np.column_stack([
        l100 > 1.05,
        lam < thr["lam_low_MW"],
        lam > thr["lam_high_MW"],
        L[:, thr["branch_77_82"]] > 1.00,
        Vm[:, thr["bus83_col"]] > thr["Vmax_bus83"] + 1e-6,
        Lmax > 1.05,
        Lmax > 1.10,
        (Vm > thr["Vmax_nonref"][None] + 0.0025).any(axis=1),
    ])
    curve = l100[:, None] > CURVE_LEVELS[None]
    return ev, curve


def unit_limit_check(q, thr):
    """Direct unit-limit events from the stored unit outputs, for the equivalence with the imbalance thresholds."""
    return (q["Pg"] < thr["Pmin_MW"][None]).any(axis=1), (q["Pg"] > thr["Pmax_MW"][None]).any(axis=1)


def clopper_pearson(k, n, level=0.95):
    a = 1 - level
    lo = stats.beta.ppf(a / 2, k, n - k + 1) if k > 0 else 0.0
    hi = stats.beta.ppf(1 - a / 2, k + 1, n - k) if k < n else 1.0
    return float(lo), float(hi)


def counted(q, thr):
    ev, curve = indicators(q, thr)
    return ev.sum(axis=0), curve.sum(axis=0), len(ev)


def model_probabilities(model, coords, outmap, thr, count, entropy, chunk=1 << 16):
    """Weighted event probabilities of count stratified Sobol points of the model (weights pi_k / N_k^s), and
    per-component probabilities."""
    Z, w, alloc = pg.model_points(model, count, entropy)
    comp = np.repeat(np.arange(model.K), alloc)
    pe, pc = np.zeros(len(EVENTS)), np.zeros(len(CURVE_LEVELS))
    hits = np.zeros(len(EVENTS), int)
    per_comp = np.zeros((model.K, len(EVENTS)))
    for s in range(0, len(Z), chunk):
        q = outmap.evaluate(coords.to_x(Z[s:s + chunk]))
        ev, curve = indicators(q, thr)
        ws = w[s:s + chunk]
        pe += ws @ ev
        pc += ws @ curve
        hits += ev.sum(axis=0)
        np.add.at(per_comp, comp[s:s + chunk], ev * ws[:, None])
    return {"events": pe, "curve": pc, "point_hits": hits, "per_component": per_comp, "allocation": alloc}
