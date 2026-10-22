"""Stage F: expansion-point batch power flow (01_IDEA/derivation.md) against point-wise Newton-Raphson on the
float64 modeling scenarios of stage D. Two runs per operating point: (a) the full procedure with reactive-limit
switching rounds and grouping by the current unit-mode pattern; (b) the final unit modes of the Newton-Raphson
solution given in advance, which isolates the expansion-point correction from the mode grouping."""

import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from tmm_cloud import batch_powerflow as bp
from tmm_cloud import matpower_case as mc
from tmm_cloud import montecarlo as mcarlo
from tmm_cloud import network as nw
from tmm_cloud import powerflow as pf
from tmm_cloud import run_record as rr

TASK = "phase1_rts118"
SET = "modeling"
V_ERR_TOL = 1e-7
NR_SUB_BATCH = 200


def find_stage_d():
    base = rr.CLOUD_ROOT / "results" / TASK
    cands = sorted(p for p in base.glob("*/stage_d_montecarlo") if all((p / o / "run_manifest.json").exists() for o in ("OP1", "OP2")))
    return cands[-1]


def known_mode_run(solver, Pd, Qd, Pren, mode):
    s = solver.pf
    P0, Q0 = s.injections(Pd, Qd, Pren)
    B = Pd.shape[0]
    isPV = np.zeros((B, s.n), bool)
    isPV[:, s.reg_bus] = mode == 0
    Qsp = Q0.copy()
    lim = np.where(mode == 1, s.Qmax_bus[s.reg_bus][None], s.Qmin_bus[s.reg_bus][None])
    rows, cols = np.nonzero(mode != 0)
    Qsp[rows, s.reg_bus[cols]] = lim[rows, cols] - Qd[rows, s.reg_bus[cols]]
    keys, inv = np.unique(mode, axis=0, return_inverse=True)
    inv = inv.ravel()
    x = np.full((B, 2 * (s.n - 1) + 1), np.nan)
    ok = np.zeros(B, bool)
    method = np.full(B, -1)
    sizes, n_exp = [], []
    for g in range(len(keys)):
        idx = np.flatnonzero(inv == g)
        xg, okg, mg, ne = solver.solve_group(P0[idx], Qsp[idx], isPV[idx])
        x[idx], ok[idx], method[idx] = xg, okg, mg
        sizes.append(len(idx))
        n_exp.append(ne)
    V, lam = solver._x_to_V(np.nan_to_num(x))
    return {"V": V, "ok": ok, "method": method, "group_sizes": np.array(sizes), "expansion_points": np.array(n_exp)}


def analyse(net, op, name, ddir, out):
    D = mcarlo.load_set(sorted((ddir / name).glob(f"{SET}_chunk*.npz")))
    ok_ref = D["converged"]
    k = D["load_factor"][ok_ref]
    Pd, Qd = op.Pd[None] * k, op.Qd[None] * k
    Pren = D["Pren_MW"][ok_ref] / net.baseMVA
    V_ref, mode_ref = D["V"][ok_ref], D["mode"][ok_ref].astype(int)
    B = len(k)

    nr = pf.BatchSolver(net, op)
    t0 = time.time()
    outs = [nr.solve(Pd[i:i + NR_SUB_BATCH], Qd[i:i + NR_SUB_BATCH], Pren[i:i + NR_SUB_BATCH]) for i in range(0, B, NR_SUB_BATCH)]
    t_nr = time.time() - t0
    V_nr = np.concatenate([o["V"] for o in outs])
    res = {"scenarios": int(B), "newton_raphson_time_s": t_nr,
           "newton_raphson_rerun_max_abs_V_diff_vs_stage_d": float(np.abs(V_nr - V_ref).max())}

    bs = bp.BatchExpansionSolver(net, op)
    t0 = time.time()
    r = bs.solve(Pd, Qd, Pren)
    t_b = time.time() - t0
    okb = r["fail"] == pf.FAIL_OK
    dV = np.abs(r["V"] - V_ref).max(axis=1)
    same_mode = (r["mode"] == mode_ref).all(axis=1)
    gs = np.array(bs.stats["group_sizes"])
    res["full_procedure"] = {
        "time_s": t_b, "speed_ratio_vs_newton_raphson": t_nr / t_b,
        "converged": int(okb.sum()), "identical_unit_modes": int((same_mode & okb).sum()),
        "max_abs_V_error_pu": float(dV[okb].max()), "scenarios_V_error_gt_tol": int((dV[okb] > V_ERR_TOL).sum()),
        "mode_groups_over_all_rounds": int(bs.stats["mode_groups"]),
        "groups_with_one_scenario": int((gs == 1).sum()), "median_group_size": float(np.median(gs)),
        "expansion_points": int(bs.stats["expansion_points"]),
        "newton_point_solutions": int(bs.stats["newton_point_solutions"]),
        "final_method_share": {m: float((r["method"][okb] == c).mean()) for c, m in ((0, "batch_correction"), (1, "expansion_point"), (2, "point_newton"))},
    }

    bs2 = bp.BatchExpansionSolver(net, op)
    t0 = time.time()
    r2 = known_mode_run(bs2, Pd, Qd, Pren, mode_ref)
    t_k = time.time() - t0
    dV2 = np.abs(r2["V"] - V_ref).max(axis=1)
    sizes = r2["group_sizes"]
    res["known_final_modes"] = {
        "time_s": t_k, "speed_ratio_vs_newton_raphson": t_nr / t_k,
        "converged": int(r2["ok"].sum()), "max_abs_V_error_pu": float(dV2[r2["ok"]].max()),
        "scenarios_V_error_gt_tol": int((dV2[r2["ok"]] > V_ERR_TOL).sum()),
        "distinct_final_modes": int(len(sizes)), "modes_with_one_scenario": int((sizes == 1).sum()),
        "largest_mode_share": float(sizes.max() / B), "scenarios_in_modes_with_ge_50": int(sizes[sizes >= 50].sum()),
        "expansion_points": int(r2["expansion_points"].sum()),
        "newton_point_solutions": int(bs2.stats["newton_point_solutions"]),
        "final_method_share": {m: float((r2["method"][r2["ok"]] == c).mean()) for c, m in ((0, "batch_correction"), (1, "expansion_point"), (2, "point_newton"))},
    }
    big = sizes >= 50
    if big.any():
        m_big = np.isin(np.unique(mode_ref, axis=0, return_inverse=True)[1].ravel(), np.flatnonzero(big))
        res["known_final_modes"]["batch_correction_share_in_modes_with_ge_50"] = float((r2["method"][m_big] == 0).mean())

    pd.DataFrame({"scenario": np.arange(B), "full_V_error_pu": dV, "full_identical_mode": same_mode,
                  "full_method": r["method"], "known_mode_V_error_pu": dV2, "known_mode_method": r2["method"]}
                 ).to_csv(out / f"batch_vs_newton_{name}.csv", index=False)

    fig, ax = plt.subplots(figsize=(6, 4))
    s_sorted = np.sort(sizes)[::-1]
    ax.loglog(np.arange(1, len(s_sorted) + 1), s_sorted, ".-")
    ax.set_xlabel("unit reactive-limit mode pattern, ranked")
    ax.set_ylabel("scenarios with this final pattern")
    ax.set_title(f"{name}: {len(sizes)} patterns in {B} modeling scenarios")
    fig.tight_layout()
    fig.savefig(out / f"mode_pattern_sizes_{name}.png", dpi=130)
    plt.close(fig)
    return res


def main():
    out = rr.new_result_dir(TASK, "stage_f_batch_powerflow")
    ddir = find_stage_d()
    bdir = rr.latest_stage_dir(TASK, "stage_b_network")
    net = nw.from_matpower(mc.PGLIB118)
    summary = {"stage_d_result": rr.relpath(ddir), "scenario_set": SET}
    for name in ("OP1", "OP2"):
        op = mcarlo.load_operating_point(bdir / f"operating_point_{name}.npz")
        summary[name] = analyse(net, op, name, ddir, out)
        print(name, summary[name], flush=True)
    rr.write_json(out / "stage_f_summary.json", summary)
    rr.write_manifest(
        out,
        command="python cloud/tests/run_stage_f_batch_pf.py",
        data_files=[mc.PGLIB118, bdir / "operating_point_OP1.npz", bdir / "operating_point_OP2.npz"],
        parameters={"stage_d_result": rr.relpath(ddir), "scenario_set": f"{SET} (float64 storage)",
                    "eps_y_pu": bp.EPS_Y, "eps_x_pu": bp.EPS_X, "max_corrections": bp.MAX_CORR,
                    "max_expansion_points_per_group": bp.MAX_EXPANSION, "step_lengths": bp.STEPS.tolist(),
                    "voltage_error_tolerance_pu": V_ERR_TOL, "newton_raphson_sub_batch": NR_SUB_BATCH,
                    "threads": "one BLAS thread, one process, for both methods"},
        random_settings="not applicable (stored stage D scenarios)",
    )
    print(rr.relpath(out))


if __name__ == "__main__":
    main()
