"""Stage I: marginal accuracy of the selected GMM, standard TMM and block TMM against the 100000 heavy-tailed
AC Monte Carlo check scenarios of OP2 (derivation.md Sections 4.1-4.3 and eq. (29)).

Usage: python cloud/tests/run_stage_i_marginals.py --result-stamp YYYYmmddTHHMMSS
"""

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from tmm_cloud import matpower_case as mc
from tmm_cloud import mixture as mx
from tmm_cloud import montecarlo as mcarlo
from tmm_cloud import network as nw
from tmm_cloud import propagation as pg
from tmm_cloud import run_record as rr
from tmm_cloud import state_model as sm

TASK = "phase1b_tmm_op2"
OP = "OP2"
N_POINTS = 1 << 20
SOBOL_ENTROPY = 2026102201
METHODS = ["GMM", "TMM", "block TMM", "two-strength block TMM"]
GROUPS = [["Vm", "Va", "lam", "Pg"], ["Pf", "Qf"], ["Pt", "Qt"], ["loading"]]
FIXED_TOL = {"Vm": 1e-6, "Va": 1e-5, "Pf": 1e-3, "Qf": 1e-3, "Pt": 1e-3, "Qt": 1e-3, "loading": 1e-6, "lam": 1e-6, "Pg": 1e-6}
UNITS = {"Vm": "p.u.", "Va": "deg", "Pf": "MW", "Qf": "Mvar", "Pt": "MW", "Qt": "Mvar", "loading": "fraction of rateA", "lam": "MW", "Pg": "MW"}
COL_BLOCK = 24
STAGE_D = rr.CLOUD_ROOT / "results" / "phase1_rts118" / "20261021T202200" / "stage_d_montecarlo" / OP
STAGE_B = rr.CLOUD_ROOT / "results" / "phase1_rts118" / "20261021T201240" / "stage_b_network"


def object_names(net):
    nr = np.setdiff1d(np.arange(net.n), [net.ref])
    br = [f"{net.bus_id[f]}-{net.bus_id[t]}" for f, t in zip(net.f, net.t)]
    return {"Vm": [f"bus{net.bus_id[i]}" for i in nr], "Va": [f"bus{net.bus_id[i]}" for i in nr],
            "Pf": br, "Qf": br, "Pt": br, "Qt": br, "loading": br, "lam": ["system"],
            "Pg": [f"unit_bus{net.bus_id[b]}" for b in net.gen_bus[net.Pmax > 0]]}


def load_stage_h(hdir):
    z = np.load(hdir / f"coordinates_{OP}.npz")
    coords = mx.Coordinates(z["center"], z["basis"], z["scale"], z["eigenvalues"], 1e-10, 0.0)
    models = {}
    for p in sorted(hdir.glob(f"model_*_{OP}.npz")):
        m = mx.from_npz(np.load(p))
        name = {"mixture": None}.get(m.family)
        if m.family == "gmm":
            name = "GMM"
        elif m.family == "tmm":
            name = "TMM"
        elif m.family == "block-tmm":
            name = "block TMM"
        elif m.family == "two-block-tmm":
            name = "two-strength block TMM"
        models[name] = (m, p)
    return coords, models


def metrics_rows(Y, w, R, key, names, method, live, unit_bus=None):
    """unit_bus: optional (is unit bus, share of check scenarios at the voltage setpoint) per object."""
    rows = []
    for s in range(0, Y.shape[1], COL_BLOCK):
        cols = [j for j in range(s, min(Y.shape[1], s + COL_BLOCK)) if live[j]]
        if not cols:
            continue
        met = pg.distribution_metrics(Y[:, cols], w, R[:, cols])
        for i, j in enumerate(cols):
            row = {"method": method, "quantity": key, "object": names[j], "unit": UNITS[key]}
            row.update({k: float(v[i]) for k, v in met.items()})
            if unit_bus is not None:
                row["unit_bus"] = bool(unit_bus[0][j])
                row["check_share_at_setpoint"] = float(unit_bus[1][j])
            rows.append(row)
    return rows


def summarise(df):
    out = []
    df = df.copy()
    df["group"] = df.quantity
    if "unit_bus" in df:
        df.loc[(df.quantity == "Vm") & (df.unit_bus == True), "group"] = "Vm unit buses"
        df.loc[(df.quantity == "Vm") & (df.unit_bus != True), "group"] = "Vm PQ buses"
    for (method, key), g in df.groupby(["method", "group"], sort=False):
        row = {"method": method, "quantity": key, "objects": len(g)}
        a = {"mean_error_in_std": (g.mean_error / g.std_ref).abs(), "std_rel_error": g.std_rel_error.abs(), "ks": g.ks}
        for p in ("q0p1", "q1", "q99", "q99p9"):
            a[f"{p}_error_in_std"] = g[f"{p}_error_in_std"].abs()
        for k, v in a.items():
            row[f"{k}_median"] = float(v.median())
            row[f"{k}_p95"] = float(v.quantile(0.95))
        w = g.loc[g.ks.idxmax()]
        row["worst_object_by_ks"] = w.object
        row["worst_ks"] = float(w.ks)
        out.append(row)
    return pd.DataFrame(out)


def plot_object(net, key, col, label, models, coords, outmap, R, Ymod, out):
    vals = {}
    for i, method in enumerate(METHODS):
        if method not in models:
            continue
        Z, w, _ = pg.model_points(models[method][0], N_POINTS, SOBOL_ENTROPY + i)
        vals[method] = (pg.evaluate_points(Z, coords, outmap, [key])[key][:, col], w / w.sum())
    r = R[:, col]
    lo, hi = np.quantile(r, [0.0005, 0.9995])
    span = hi - lo
    bins = np.linspace(lo - 0.15 * span, hi + 0.15 * span, 121)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(r, bins=bins, density=True, color="0.8", label="AC Monte Carlo check (100000)")
    styles = {"GMM": ("tab:blue", "--"), "TMM": ("tab:red", "-"), "block TMM": ("tab:green", "-."), "two-strength block TMM": ("tab:purple", ":")}
    for method, (y, w) in vals.items():
        h, e = np.histogram(y, bins=bins, weights=w, density=True)
        ax.plot(0.5 * (e[1:] + e[:-1]), h, color=styles[method][0], ls=styles[method][1], label=method)
    ax.set_xlabel(f"{label} ({UNITS[key]})")
    ax.set_ylabel("probability density")
    ax.legend(fontsize=7)
    fig.tight_layout()
    tag = label.replace(" ", "_").replace("/", "_")
    fig.savefig(out / f"pdf_{tag}.png", dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4))
    rs = np.sort(r)
    n = len(rs)
    ax.semilogy(rs, 1 - np.arange(1, n + 1) / n + 1e-12, color="k", label="check, P(Y > y)")
    ax.semilogy(rs, np.arange(1, n + 1) / n, color="k", ls="--", label="check, P(Y < y)")
    ms = np.sort(Ymod[:, col])
    nm = len(ms)
    ax.semilogy(ms, 1 - np.arange(1, nm + 1) / nm + 1e-12, color="0.5", lw=0.8, label="modeling scenarios (5256), P(Y > y)")
    ax.semilogy(ms, np.arange(1, nm + 1) / nm, color="0.5", lw=0.8, ls="--")
    for method, (y, w) in vals.items():
        o = np.argsort(y)
        cw = np.cumsum(w[o])
        c, ls = styles[method]
        ax.semilogy(y[o], np.clip(1 - cw, 1e-12, None), color=c, ls=ls, lw=1, label=f"{method}, P(Y > y)")
        ax.semilogy(y[o], cw, color=c, ls=ls, lw=1, alpha=0.6)
    ax.set_ylim(1e-6, 1)
    ax.set_xlim(rs[0] - 0.1 * span, rs[-1] + 0.1 * span)
    ax.set_xlabel(f"{label} ({UNITS[key]})")
    ax.set_ylabel("probability (upper tail solid, lower tail dashed)")
    ax.legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(out / f"exceedance_{tag}.png", dpi=130)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--result-stamp", required=True)
    a = ap.parse_args()
    t0 = time.time()
    out = rr.stamped_stage_dir(TASK, a.result_stamp, "stage_i_marginals")
    hdir = rr.CLOUD_ROOT / "results" / TASK / a.result_stamp / "stage_h_fit"
    net = nw.from_matpower(mc.PGLIB118)
    op = mcarlo.load_operating_point(STAGE_B / f"operating_point_{OP}.npz")
    outmap = sm.OutputMap(net, op)
    coords, models = load_stage_h(hdir)
    names = object_names(net)

    check = outmap.from_stage_d(mcarlo.load_set(sorted(STAGE_D.glob("check_heavy_chunk*.npz"))))
    _, dmod = sm.scenario_state(sorted(STAGE_D.glob("modeling_chunk*.npz")), net)
    mod = outmap.from_stage_d(dmod)
    live = {k: check[k].std(axis=0) >= FIXED_TOL[k] for k in check}
    fixed = [{"quantity": k, "object": names[k][j], "check_std": float(check[k][:, j].std())} for k in check for j in np.flatnonzero(~live[k])]
    pd.DataFrame(fixed).to_csv(out / f"fixed_objects_{OP}.csv", index=False)

    # Vm of a unit bus equals the unit setpoint while the unit regulates: a point mass the continuous models cannot
    # carry; such objects are summarised separately (reporting split, no change of the models)
    nr = np.setdiff1d(np.arange(net.n), [net.ref])
    vset = np.full(net.n, np.nan)
    vset[net.gen_bus] = op.Vset
    is_unit = np.isin(nr, net.gen_bus)
    at_set = (np.abs(check["Vm"] - vset[nr][None]) < 1e-6).mean(axis=0)
    ub = (is_unit, at_set)
    rows, timing, alloc_info = [], {}, {}
    nm = len(mod["lam"])
    for key in check:
        rows += metrics_rows(mod[key], np.full(nm, 1.0 / nm), check[key], key, names[key], "modeling scenarios (5256)", live[key], ub if key == "Vm" else None)
    for i, method in enumerate(METHODS):
        if method not in models:
            continue
        ts = time.time()
        Z, w, alloc = pg.model_points(models[method][0], N_POINTS, SOBOL_ENTROPY + i)
        alloc_info[method] = alloc.tolist()
        for group in GROUPS:
            Y = pg.evaluate_points(Z, coords, outmap, group)
            for key in group:
                rows += metrics_rows(Y[key], w, check[key], key, names[key], method, live[key], ub if key == "Vm" else None)
            del Y
        timing[method] = time.time() - ts
        print(method, f"{timing[method]:.0f} s", flush=True)
        del Z
    df = pd.DataFrame(rows)
    df.to_csv(out / f"marginal_metrics_{OP}.csv", index=False)
    summ = summarise(df)
    summ.to_csv(out / f"marginal_summary_{OP}.csv", index=False)

    br = {f"{net.bus_id[f]}-{net.bus_id[t]}": j for j, (f, t) in enumerate(zip(net.f, net.t))}
    nr = list(nr)
    plots = [("loading", br["100-103"], "branch 100-103 loading"), ("Pf", br["12-14"], "branch 12-14 P"),
             ("loading", br["77-82"], "branch 77-82 loading"), ("Vm", nr.index(int(net.idx([83])[0])), "bus 83 Vm"),
             ("lam", 0, "system imbalance")]
    top = df[df.method == "TMM"].sort_values("ks", ascending=False)
    top = top[~((top.quantity == "Vm") & (top.unit_bus == True))]
    chosen = {(k, names[k][c]) for k, c, _ in plots}
    extra = []
    for _, r in top.iterrows():
        if (r.quantity, r.object) not in chosen and len(extra) < 3:
            extra.append((r.quantity, names[r.quantity].index(r.object), f"{r.quantity} {r.object} (TMM KS {r.ks:.3f})"))
            chosen.add((r.quantity, r.object))
    for key, col, label in plots + extra:
        plot_object(net, key, col, label, models, coords, outmap, check[key], mod[key], out)

    summary = {"operating_point": OP, "points_per_model": N_POINTS, "component_allocation": alloc_info,
               "propagation_seconds": timing, "fixed_objects": len(fixed), "plotted_objects": [p[2] for p in plots + extra],
               "models": {m: rr.relpath(p) for m, (_, p) in models.items()}, "wall_time_s": time.time() - t0}
    rr.write_json(out / f"stage_i_summary_{OP}.json", summary)
    rr.write_manifest(
        out,
        command=f"python cloud/tests/run_stage_i_marginals.py --result-stamp {a.result_stamp}",
        data_files=[mc.PGLIB118, STAGE_B / f"operating_point_{OP}.npz", hdir / f"coordinates_{OP}.npz",
                    *[p for _, p in models.values()], *sorted(STAGE_D.glob("modeling_chunk*.npz"))],
        parameters={"reference": "check_heavy (100000, float32 storage), phase-1 stage D OP2; chunk SHA256 in stage_g_identity.json",
                    "floor": "empirical distribution of the 5256 modeling scenarios", "points_per_model": N_POINTS,
                    "quantities": "Vm, Va of 117 non-reference buses; P, Q at both ends of 186 branches; loading max(|Sf|,|St|)/rateA; imbalance (MW); Pg of 19 units with Pmax > 0",
                    "fixed_object_std_tolerance": FIXED_TOL, "quantile_probabilities": pg.PROBS,
                    "quantile_convention": "model: smallest value with cumulative weight >= p; check: inverted_cdf",
                    "ks": "sup over the check sample points of |F_model - F_check|, both one-sided limits",
                    "reporting_split": "Vm of unit buses (point mass at the setpoint while regulating) summarised apart from Vm of PQ buses; the 3 extra plotted objects are the largest-TMM-KS objects outside Vm of unit buses"},
        random_settings={"sobol": "scipy.stats.qmc.Sobol(scramble=True): linear matrix scramble with digital random shift; dimension d+1 per component (GMM d)",
                         "entropy": {m: SOBOL_ENTROPY + i for i, m in enumerate(METHODS)},
                         "per_component": "numpy default_rng(SeedSequence(entropy).spawn(K)[k])",
                         "allocation": "N_k = 2^20 pi_k by largest remainder, at least 1"},
    )
    print(rr.relpath(out), f"{time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
