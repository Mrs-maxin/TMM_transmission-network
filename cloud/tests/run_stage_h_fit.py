"""Stage H: GMM, standard TMM and block TMM of the OP2 power-flow state (derivation.md Section 3).

State x = [e, f of the 117 non-reference buses, system imbalance lam in p.u.] (235). Estimation on the float64
modeling set (5256), component count by eq. (22) on the determination set (1752). Block TMM: one strength over
the active-power rows of the 10 renewable host buses (main), and a two-strength variant (wind rows, PV rows).

Usage: python cloud/tests/run_stage_h_fit.py --result-stamp YYYYmmddTHHMMSS [--workers 4]
"""

import os

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import argparse
import functools
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
from tmm_cloud import run_record as rr
from tmm_cloud import state_model as sm

TASK = "phase1b_tmm_op2"
OP = "OP2"
KS = [1, 2, 3, 4, 6, 8]
METHODS = ["GMM", "TMM", "block TMM", "two-strength block TMM"]
STAGE_D = rr.CLOUD_ROOT / "results" / "phase1_rts118" / "20261021T202200" / "stage_d_montecarlo" / OP
STAGE_B = rr.CLOUD_ROOT / "results" / "phase1_rts118" / "20261021T201240" / "stage_b_network"


def fit_all(net, op, mapping, X, Xv, modeling, workers, verbose=True):
    coords = mx.fit_coordinates(X)
    Z, Zv = coords.to_z(X), coords.to_z(Xv)
    dirs, dinfo = sm.strength_setup(net, op, mapping, modeling, coords)
    extra = {"two-strength block TMM": functools.partial(mx.fit_two_block_tmm, Q1=dirs["Q1"], Q2=dirs["Q2"], Qw=dirs["Qw2"])}
    t0 = time.time()
    best, selected, rows = mx.fit_orders(Z, Zv, KS, dirs["Qu"], dirs["Qw"], workers=workers, extra_block=extra, verbose=verbose)
    return coords, dirs, dinfo, best, selected, rows, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--result-stamp", required=True)
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    out = rr.stamped_stage_dir(TASK, a.result_stamp, "stage_h_fit")
    net = nw.from_matpower(mc.PGLIB118)
    op = mcarlo.load_operating_point(STAGE_B / f"operating_point_{OP}.npz")
    mapping = pd.read_csv(STAGE_B / "renewable_mapping.csv")
    mod_files = sorted(STAGE_D.glob("modeling_chunk*.npz"))
    det_files = sorted(STAGE_D.glob("determination_chunk*.npz"))
    X, modeling = sm.scenario_state(mod_files, net)
    Xv, _ = sm.scenario_state(det_files, net)

    coords, dirs, dinfo, best, selected, rows, secs = fit_all(net, op, mapping, X, Xv, modeling, a.workers)
    pd.DataFrame(rows).to_csv(out / f"fit_candidates_{OP}.csv", index=False)

    table = []
    for (method, K), m in sorted(best.items(), key=lambda t: (METHODS.index(t[0][0]), t[0][1])):
        au = m.audit
        row = {"method": method, "K": K, "selected": selected.get(method) == K,
               "modeling_mean_log_likelihood": au["mean_log_likelihood"],
               "determination_mean_log_likelihood": au["determination_mean_log_likelihood"],
               "iterations": au["iterations"], "seconds": au["seconds"], "restart_direction": au["restart_direction"],
               "initial_df": au.get("initial_df"), "weights": np.round(m.weights, 5).tolist(),
               "effective_counts": np.round(au["effective_counts"], 1).tolist(), "df": np.round(m.df, 3).tolist()}
        if hasattr(m, "df2"):
            row["df_pv"] = np.round(m.df2, 3).tolist()
        table.append(row)
    tab = pd.DataFrame(table)
    tab.to_csv(out / f"fit_by_order_{OP}.csv", index=False)

    np.savez_compressed(out / f"coordinates_{OP}.npz", center=coords.center, basis=coords.basis, scale=coords.scale,
                        eigenvalues=coords.eigenvalues)
    np.savez_compressed(out / f"strength_directions_{OP}.npz", **{k: dirs[k] for k in ("Qu", "Qw", "Q1", "Q2", "Qw2", "H", "mode")})
    for method, K in selected.items():
        tag = method.replace(" ", "_").replace("-", "_")
        np.savez_compressed(out / f"model_{tag}_K{K}_{OP}.npz", **mx.to_npz_dict(best[(method, K)]))

    fig, ax = plt.subplots(figsize=(6, 4))
    for method, mk in zip(METHODS, "osd^"):
        g = tab[tab.method == method]
        ax.plot(g.K, g.determination_mean_log_likelihood, marker=mk, label=method)
    ax.set_xlabel("number of components K")
    ax.set_ylabel("determination-set mean log-likelihood")
    ax.set_title(f"{OP}: d = {coords.d}, r = {dinfo['r']}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / f"order_selection_{OP}.png", dpi=130)
    plt.close(fig)

    summary = {"operating_point": OP, "modeling_scenarios": int(len(X)), "determination_scenarios": int(len(Xv)),
               "state_dimension": int(X.shape[1]), "d": int(coords.d),
               "dropped_variance_share": float(np.clip(coords.eigenvalues[coords.d:], 0, None).sum() / np.clip(coords.eigenvalues, 0, None).sum()),
               "smallest_kept_eigenvalue": float(coords.eigenvalues[coords.d - 1]),
               "largest_dropped_eigenvalue": float(coords.eigenvalues[coords.d]) if coords.d < len(coords.eigenvalues) else None,
               "directions": dinfo, "selected_K": selected, "fit_wall_time_s": secs,
               "candidate_counts": pd.DataFrame(rows).groupby(["method", "status"]).size().reset_index().values.tolist() if rows else []}
    rr.write_json(out / f"stage_h_summary_{OP}.json", summary)
    rr.write_manifest(
        out,
        command=f"python cloud/tests/run_stage_h_fit.py --result-stamp {a.result_stamp} --workers {a.workers}",
        data_files=[mc.PGLIB118, STAGE_B / f"operating_point_{OP}.npz", STAGE_B / "renewable_mapping.csv", *mod_files, *det_files],
        parameters={"state": "x = [e(117 non-reference buses), f(117), lam (p.u.)]", "standardization": "eq. (11), dropped variance share <= 1e-10, eigenvalues > p*eps*lambda_1",
                    "candidate_K": KS, "covariance_floor": mx.COV_FLOOR, "df_bounds": mx.DF_BOUNDS,
                    "stopping": f"absolute change of the mean log-likelihood < {mx.LL_TOL}, at most {mx.MAX_ITER} iterations",
                    "decrease_tolerance": mx.DECREASE_TOL, "initial_df": mx.INITIAL_DF, "initial_splits": "principal directions 1, 2, 3, K equal-count groups",
                    "acceptance": "converged, no decrease, N_k > d + 2", "per_K_choice": "highest modeling mean log-likelihood among accepted fits",
                    "K_selection": "eq. (22), determination mean log-likelihood",
                    "block_L_W": "active-power rows of the 10 renewable host buses (stage-B renewable_mapping.csv)",
                    "jacobian_point": "most frequent unit mode of the modeling set, mean modeling-set load factors and renewable injections",
                    "two_strength_variant": "wind rows and PV rows each with an independent inverse-gamma strength within a component; PV directions orthogonalized against wind directions",
                    "workers": a.workers},
        random_settings="not applicable (deterministic initial splits; stored stage D scenarios)",
    )
    print(rr.relpath(out), selected, f"{secs:.0f} s")


if __name__ == "__main__":
    main()
