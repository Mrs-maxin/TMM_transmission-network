"""Stage E: propagation of the heavy-tailed renewable input to bus voltages and branch flows, violation
probabilities with 95% binomial intervals, and conventional units leaving [Pmin, Pmax]."""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from tmm_cloud import matpower_case as mc
from tmm_cloud import montecarlo as mcarlo
from tmm_cloud import network as nw
from tmm_cloud import powerflow as pf
from tmm_cloud import run_record as rr
from tmm_cloud import tail_stats as ts

TASK = "phase1_rts118"
V_TOL = 1e-6
FIXED_TOL = {"Vm": 1e-6, "Va": 1e-5, "Pf": 1e-3, "Qf": 1e-3, "Pt": 1e-3, "Qt": 1e-3}
N_SUBSETS = 5
OVER_DEPTHS = [0.0, 0.0025, 0.005, 0.01, 0.02]
LOAD_DEPTHS = [0.0, 0.05, 0.10, 0.20, 0.30]


def find_stage_d():
    base = rr.CLOUD_ROOT / "results" / TASK
    cands = sorted(p for p in base.glob("*/stage_d_montecarlo") if all((p / o / "run_manifest.json").exists() for o in ("OP1", "OP2")))
    if not cands:
        raise FileNotFoundError("stage D results for OP1 and OP2 not found")
    return cands[-1]


def object_blocks(net, d):
    nr = np.setdiff1d(np.arange(net.n), [net.ref])
    br = [f"{net.bus_id[f]}-{net.bus_id[t]}#{k}" for k, (f, t) in enumerate(zip(net.f, net.t))]
    return [
        ("Vm", [f"bus{net.bus_id[i]}" for i in nr], lambda: np.abs(d["V"][:, nr]).astype(np.float64), "p.u."),
        ("Va", [f"bus{net.bus_id[i]}" for i in nr], lambda: np.rad2deg(np.angle(d["V"][:, nr])).astype(np.float64), "deg"),
        ("Pf", br, lambda: d["Pf_MW"].astype(np.float64), "MW"),
        ("Qf", br, lambda: d["Qf_Mvar"].astype(np.float64), "Mvar"),
        ("Pt", br, lambda: d["Pt_MW"].astype(np.float64), "MW"),
        ("Qt", br, lambda: d["Qt_Mvar"].astype(np.float64), "Mvar"),
    ]


def most_frequent_mode(mode):
    u, inv, cnt = np.unique(mode, axis=0, return_inverse=True, return_counts=True)
    order = np.argsort(-cnt)
    return u[order[0]], cnt[order] / len(mode), len(u)


def sensitivity_at_mode(net, op, mode_vec):
    solver = pf.BatchSolver(net, op)
    P0, Q0 = solver.injections(op.Pd[None], op.Qd[None], op.Pren[None])
    isPV = np.zeros((1, net.n), bool)
    isPV[0, solver.reg_bus] = mode_vec == 0
    Qsp = Q0.copy()
    lim = np.where(mode_vec == 1, solver.Qmax_bus[solver.reg_bus], solver.Qmin_bus[solver.reg_bus])
    sel = mode_vec != 0
    Qsp[0, solver.reg_bus[sel]] = lim[sel] - op.Qd[solver.reg_bus[sel]]
    V = np.tile(solver.Vset_bus.astype(complex), (1, 1))
    V[:, net.ref] = solver.Vref
    V, lam, it, fl = solver.newton(V, np.zeros(1), P0, Qsp, isPV)
    if fl[0] != pf.FAIL_OK:
        raise RuntimeError("sensitivity operating point did not converge")
    return ts.renewable_sensitivity(solver, V[0], lam[0], isPV[0], Qsp[0], P0[0]), V[0]


def analyse(net, op, name, ddir, out):
    t0 = time.time()
    H = mcarlo.load_set(sorted((ddir / name).glob("check_heavy_chunk*.npz")))
    L = mcarlo.load_set(sorted((ddir / name).glob("check_light_chunk*.npz")))
    res = {"n_heavy": int(len(H["converged"])), "n_light": int(len(L["converged"])),
           "non_converged_heavy": int((~H["converged"]).sum()), "non_converged_light": int((~L["converged"]).sum())}
    for D in (H, L):
        ok = D["converged"]
        for k in list(D):
            D[k] = D[k][ok]
    mode_vec, mode_share, n_modes = most_frequent_mode(H["mode"])
    res["distinct_modes_heavy"] = int(n_modes)
    res["top_mode_shares_heavy"] = mode_share[:5].round(4).tolist()
    reg_gen = np.flatnonzero(net.gen_bus != net.ref)
    at_lim = (H["mode"] != 0).mean(axis=0)
    res["units_at_Q_limit_frequency_heavy"] = {int(net.bus_id[net.gen_bus[g]]): round(float(f), 4) for g, f in zip(reg_gen, at_lim) if f > 0}

    sens, V0 = sensitivity_at_mode(net, op, mode_vec)
    Xren_h = H["Pren_MW"]
    Cov_ren = np.cov(Xren_h, rowvar=False)
    idx = np.arange(len(H["converged"]))
    subsets = np.array_split(idx, N_SUBSETS)

    rows = []
    for key, names, getter_h, unit in object_blocks(net, H):
        Xh = getter_h()
        Xl = dict((k, g) for k, _, g, _ in object_blocks(net, L))[key]()
        sh = ts.column_statistics(Xh)
        sl = ts.column_statistics(Xl)
        sub_tail = np.array([ts.column_statistics(Xh[s])["tail_index"] for s in subsets])
        sub_kurt = np.array([ts.column_statistics(Xh[s])["excess_kurtosis"] for s in subsets])
        s = sens["d" + key]
        if key in ("Vm", "Va"):
            s = s[np.setdiff1d(np.arange(net.n), [net.ref])]
        ren_var = np.einsum("ij,jk,ik->i", s, Cov_ren, s)
        with np.errstate(invalid="ignore", divide="ignore"):
            ren_share = ren_var / sh["std"] ** 2
        fixed = sh["std"] < FIXED_TOL[key]
        for j, nm in enumerate(names):
            row = {"quantity": key, "object": nm, "unit": unit, "fixed": bool(fixed[j])}
            for stat in sh:
                row[f"{stat}_heavy"] = sh[stat][j]
                row[f"{stat}_light"] = sl[stat][j]
            row["tail_index_subset_sd"] = sub_tail[:, j].std(ddof=1)
            row["kurtosis_subset_sd"] = sub_kurt[:, j].std(ddof=1)
            row["renewable_sensitivity_norm"] = float(np.linalg.norm(s[j]))
            row["renewable_variance_share_linear"] = ren_share[j]
            rows.append(row)
        del Xh, Xl
    df = pd.DataFrame(rows)
    df["tail_index_diff"] = df.tail_index_heavy - df.tail_index_light
    # light set has 20000 scenarios, one subset of the heavy set also 20000
    df["tail_index_diff_z"] = df.tail_index_diff / (df.tail_index_subset_sd * np.sqrt(1 + 1 / N_SUBSETS))
    df["kurtosis_diff"] = df.excess_kurtosis_heavy - df.excess_kurtosis_light
    df["kurtosis_diff_z"] = df.kurtosis_diff / (df.kurtosis_subset_sd * np.sqrt(1 + 1 / N_SUBSETS))
    df["q995_shift_in_light_std"] = (df.q995_heavy - df.q995_light) / df.std_light
    df["q005_shift_in_light_std"] = (df.q005_heavy - df.q005_light) / df.std_light
    df.to_csv(out / f"tail_statistics_{name}.csv", index=False)

    live = df[~df.fixed]
    res["objects"] = int(len(df))
    res["fixed_objects"] = int(df.fixed.sum())
    res["by_quantity"] = {}
    for key, g in live.groupby("quantity"):
        rs = stats.spearmanr(g.renewable_variance_share_linear, g.tail_index_diff, nan_policy="omit").statistic
        res["by_quantity"][key] = {
            "objects": int(len(g)),
            "median_excess_kurtosis_heavy": float(g.excess_kurtosis_heavy.median()),
            "median_excess_kurtosis_light": float(g.excess_kurtosis_light.median()),
            "max_excess_kurtosis_heavy": float(g.excess_kurtosis_heavy.max()),
            "median_tail_index_heavy": float(g.tail_index_heavy.median()),
            "median_tail_index_light": float(g.tail_index_light.median()),
            "objects_tail_diff_z_gt_3": int((g.tail_index_diff_z > 3).sum()),
            "objects_tail_diff_z_lt_-3": int((g.tail_index_diff_z < -3).sum()),
            "objects_kurtosis_diff_z_gt_3": int((g.kurtosis_diff_z > 3).sum()),
            "median_tail_index_diff_z": float(g.tail_index_diff_z.median()),
            "spearman_renewable_share_vs_tail_diff": float(rs),
            "median_renewable_variance_share": float(g.renewable_variance_share_linear.median()),
        }
    top = live.sort_values("tail_index_diff_z", ascending=False).groupby("quantity").head(10)
    top.to_csv(out / f"largest_tail_difference_{name}.csv", index=False)

    # violations
    viol = {}
    graded_all = []
    for tag, D in (("heavy", H), ("light", L)):
        Vm = np.abs(D["V"]).astype(np.float64)
        under = Vm < net.Vmin[None] - V_TOL
        over = Vm > net.Vmax[None] + V_TOL
        Sf = np.hypot(D["Pf_MW"].astype(np.float64), D["Qf_Mvar"].astype(np.float64))
        St = np.hypot(D["Pt_MW"].astype(np.float64), D["Qt_Mvar"].astype(np.float64))
        load = np.maximum(Sf, St) / net.rateA[None]
        ovl = load > 1.0
        Pg = D["Pg_MW"].astype(np.float64)
        share = net.Pmax > 0
        p_hi = (Pg > net.Pmax[None] * net.baseMVA + 1e-6) & share[None]
        p_lo = (Pg < net.Pmin[None] * net.baseMVA - 1e-6) & share[None]
        bus_names = [f"bus{b}" for b in net.bus_id]
        br_names = [f"{net.bus_id[f]}-{net.bus_id[t]}#{k}" for k, (f, t) in enumerate(zip(net.f, net.t))]
        tabs = {
            "undervoltage": ts.probability_table(under, bus_names),
            "overvoltage": ts.probability_table(over, bus_names),
            "branch_overload": ts.probability_table(ovl, br_names),
            "unit_above_Pmax": ts.probability_table(p_hi, [f"unit_bus{net.bus_id[b]}" for b in net.gen_bus]),
            "unit_below_Pmin": ts.probability_table(p_lo, [f"unit_bus{net.bus_id[b]}" for b in net.gen_bus]),
        }
        system = np.stack([under.any(1), over.any(1), ovl.any(1), (p_hi | p_lo).any(1),
                           under.any(1) | over.any(1) | ovl.any(1)], axis=1)
        tabs["system"] = ts.probability_table(system, ["any_undervoltage", "any_overvoltage", "any_branch_overload",
                                                       "any_unit_outside_P_limits", "any_voltage_or_branch_violation"])
        for k, t in tabs.items():
            pd.DataFrame(t).query("count > 0").to_csv(out / f"violation_{k}_{name}_{tag}.csv", index=False)
        graded = []
        dv_over = Vm - net.Vmax[None]
        dv_under = net.Vmin[None] - Vm
        for kind, depth, thresholds, unit in (("overvoltage", dv_over, OVER_DEPTHS, "p.u. above Vmax"),
                                               ("undervoltage", dv_under, OVER_DEPTHS, "p.u. below Vmin"),
                                               ("branch_overload", load - 1.0, LOAD_DEPTHS, "loading above rateA (fraction)")):
            for th in thresholds:
                ev = depth > th + (V_TOL if kind != "branch_overload" else 0.0)
                k_any = int(ev.any(axis=1).sum())
                lo_, hi_ = ts.clopper_pearson(k_any, ev.shape[0])
                per = ev.mean(axis=0)
                graded.append({"set": tag, "violation": kind, "depth_threshold": th, "unit": unit,
                               "system_count": k_any, "system_probability": k_any / ev.shape[0],
                               "system_ci95_low": float(lo_), "system_ci95_high": float(hi_),
                               "objects_with_probability_gt_0": int((per > 0).sum()),
                               "objects_with_probability_1e-3_to_1e-2": int(((per >= 1e-3) & (per <= 1e-2)).sum()),
                               "objects_with_probability_1e-4_to_1e-3": int(((per >= 1e-4) & (per < 1e-3)).sum())})
        graded_all.extend(graded)
        viol[tag] = {
            "system": tabs["system"],
            "objects_with_probability_1e-3_to_1e-2": {k: [r["object"] for r in tabs[k] if 1e-3 <= r["probability"] <= 1e-2]
                                                     for k in ("undervoltage", "overvoltage", "branch_overload")},
            "max_object_probability": {k: max(tabs[k], key=lambda r: r["probability"]) for k in tabs if k != "system"},
            "max_branch_loading_q999": float(np.quantile(load.max(axis=1), 0.999)),
            "imbalance_MW": {"std": float(D["imbalance_MW"].std()), "excess_kurtosis": float(stats.kurtosis(D["imbalance_MW"])),
                             "q005": float(np.quantile(D["imbalance_MW"], 0.005)), "q995": float(np.quantile(D["imbalance_MW"], 0.995)),
                             "min": float(D["imbalance_MW"].min()), "max": float(D["imbalance_MW"].max())},
        }
    res["violations"] = viol
    pd.DataFrame(graded_all).to_csv(out / f"graded_violation_probabilities_{name}.csv", index=False)
    res["graded_violations"] = graded_all
    res["elapsed_s"] = time.time() - t0

    # diagnostic figures
    fig, ax = plt.subplots(figsize=(6, 4))
    for key, c in (("Vm", "tab:blue"), ("Pf", "tab:orange"), ("Qf", "tab:green")):
        g = live[live.quantity == key]
        ax.scatter(g.renewable_variance_share_linear.clip(0, 1.2), g.tail_index_diff_z, s=8, c=c, label=key, alpha=0.7)
    ax.axhline(3, color="k", ls="--", lw=0.8)
    ax.axhline(-3, color="k", ls="--", lw=0.8)
    ax.set_xlabel("renewable share of variance (linear sensitivity)")
    ax.set_ylabel("tail index difference / sampling sd")
    ax.set_title(f"{name}: heavy vs light check sets")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / f"tail_difference_vs_renewable_share_{name}.png", dpi=130)
    plt.close(fig)

    for key in ("Vm", "Pf"):
        g = live[live.quantity == key].sort_values("tail_index_diff_z", ascending=False)
        if g.empty:
            continue
        obj = g.iloc[0]
        blocks = {k: (names, gh) for k, names, gh, _ in object_blocks(net, H)}
        blocks_l = {k: (names, gl) for k, names, gl, _ in object_blocks(net, L)}
        names, gh = blocks[key]
        col = names.index(obj.object)
        xh = gh()[:, col]
        xl = blocks_l[key][1]()[:, col]
        fig, ax = plt.subplots(figsize=(6, 4))
        for x, lab in ((xh, "t-type copula input (100000)"), (xl, "Gaussian copula input (20000)")):
            z = (x - np.median(x)) / (np.quantile(x, 0.75) - np.quantile(x, 0.25))
            for side, ls in ((1, "-"), (-1, "--")):
                v = np.sort(side * z)
                p = 1 - np.arange(1, len(v) + 1) / (len(v) + 1)
                ax.semilogy(v, p, ls, label=f"{lab}, {'upper' if side > 0 else 'lower'} tail")
        ax.set_xlim(0, None)
        ax.set_ylim(1e-5, 1)
        ax.set_xlabel("distance from median / IQR")
        ax.set_ylabel("exceedance probability")
        ax.set_title(f"{name}: {key} {obj.object}")
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(out / f"tail_exceedance_{key}_{name}.png", dpi=130)
        plt.close(fig)
    return res


def main():
    out = rr.new_result_dir(TASK, "stage_e_tail")
    ddir = find_stage_d()
    bdir = rr.latest_stage_dir(TASK, "stage_b_network")
    net = nw.from_matpower(mc.PGLIB118)
    summary = {"stage_d_result": rr.relpath(ddir)}
    for name in ("OP1", "OP2"):
        op = mcarlo.load_operating_point(bdir / f"operating_point_{name}.npz")
        summary[name] = analyse(net, op, name, ddir, out)
        print(name, "done", flush=True)
    rr.write_json(out / "stage_e_summary.json", summary)
    files = [mc.PGLIB118, bdir / "operating_point_OP1.npz", bdir / "operating_point_OP2.npz"]
    rr.write_manifest(
        out,
        command="python cloud/tests/run_stage_e_tail.py",
        data_files=files,
        parameters={"stage_d_result": rr.relpath(ddir), "voltage_violation_tolerance_pu": V_TOL,
                    "overload": "max(|S_from|, |S_to|) > rateA", "unit_limits": "units with Pmax > 0, Pg outside [Pmin, Pmax] by more than 1e-6 MW",
                    "tail_index": "(q99.5 - q0.5) / (q75 - q25)", "sampling_sd": f"sd of the index over {N_SUBSETS} disjoint 20000-scenario subsets of the heavy check set",
                    "fixed_object_std_tolerance": FIXED_TOL, "binomial_interval": "Clopper-Pearson 95%",
                    "sensitivity": "inverse Jacobian of the distributed-slack power-flow equations at the forecast operating point, most frequent unit mode of the heavy check set"},
        random_settings="not applicable (analysis of stored stage D scenarios)",
    )
    print(rr.relpath(out))


if __name__ == "__main__":
    main()
