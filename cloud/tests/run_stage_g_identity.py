"""Stage G: confirm that the regenerated OP2 scenario sets of phase-1 stage D are identical to the original run.

The OP2 chunks are regenerated beforehand with
    python cloud/tests/run_stage_d_montecarlo.py OP2 --result-stamp 20261021T202200 --workers 4
This entry recomputes, from the regenerated chunks, the per-set counts that the committed stage_d_summary_OP2.json
records (scenarios, converged, failures, mean and max iterations, rounds distribution), recounts two phase-1
stage E events on check_heavy, compares per-object statistics with phase-1 stage E, and records SHA256 of the
regenerated chunks.

Usage: python cloud/tests/run_stage_g_identity.py --result-stamp YYYYmmddTHHMMSS [--op OP1]
(OP1 is stage K and writes to stage_k_OP1/stage_g_identity.)
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from tmm_cloud import matpower_case as mc
from tmm_cloud import montecarlo as mcarlo
from tmm_cloud import network as nw
from tmm_cloud import run_record as rr
from tmm_cloud import tail_stats as ts

TASK = "phase1b_tmm_op2"
STAGE_D_ROOT = rr.CLOUD_ROOT / "results" / "phase1_rts118" / "20261021T202200" / "stage_d_montecarlo"
# phase-1 stage E graded_violation_probabilities_<OP>.csv: scenarios with ANY branch above the loading level
EXPECTED = {"OP2": {"any_branch_loading_gt_1.05": 153, "imbalance_below_-325.75MW": 228},
            "OP1": {"any_branch_loading_gt_1.00": 27745, "any_branch_loading_gt_1.30": 30}}
ENTROPY = {"OP1": 2026102111, "OP2": 2026102112}
STAGE_E = rr.CLOUD_ROOT / "results" / "phase1_rts118" / "20261021T203828" / "stage_e_tail"


def committed_summary(path):
    rel = rr.relpath(path)
    return json.loads(subprocess.check_output(["git", "show", f"HEAD:{rel}"], cwd=rr.PROJECT_ROOT, text=True))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--result-stamp", required=True)
    ap.add_argument("--op", default="OP2", choices=["OP1", "OP2"])
    a = ap.parse_args()
    op_name = a.op
    STAGE_D = STAGE_D_ROOT / op_name
    expected = EXPECTED[op_name]
    out = rr.stamped_stage_dir(TASK, a.result_stamp, "stage_g_identity" if op_name == "OP2" else f"stage_k_{op_name}/stage_g_identity")
    old = committed_summary(STAGE_D / f"stage_d_summary_{op_name}.json")
    rows, identical = [], True
    for set_name, info in old["sets"].items():
        d = mcarlo.load_set(sorted(STAGE_D.glob(f"{set_name}_chunk*.npz")), ["converged", "fail_code", "iterations", "rounds"])
        rk, rc = np.unique(d["rounds"], return_counts=True)
        new = {"scenarios": int(len(d["converged"])), "converged": int(d["converged"].sum()),
               "failed_scenario_index": np.flatnonzero(~d["converged"]).tolist(),
               "iterations_mean": float(d["iterations"].mean()), "iterations_max": int(d["iterations"].max()),
               "rounds_distribution": {str(int(a)): int(b) for a, b in zip(rk, rc)}}
        for key, v in new.items():
            same = info[key] == v if key != "iterations_mean" else abs(info[key] - v) < 1e-12
            identical &= bool(same)
            rows.append({"set": set_name, "quantity": key, "committed": info[key], "regenerated": v, "identical": bool(same)})

    net = nw.from_matpower(mc.PGLIB118)
    H = mcarlo.load_set(sorted(STAGE_D.glob("check_heavy_chunk*.npz")), ["V", "Pf_MW", "Qf_Mvar", "Pt_MW", "Qt_Mvar", "imbalance_MW", "converged"])
    k = int(np.flatnonzero((net.bus_id[net.f] == 100) & (net.bus_id[net.t] == 103))[0])
    load = np.maximum(np.hypot(H["Pf_MW"], H["Qf_Mvar"]), np.hypot(H["Pt_MW"], H["Qt_Mvar"])) / net.rateA[None]
    level = 1.05 if op_name == "OP2" else 1.30
    over = load > level
    counts = {"branch_100_103_loading_gt_1.05": int((load[:, k] > 1.05).sum()),
              "any_branch_loading_gt_1.00": int((load > 1.00).any(axis=1).sum()),
              "any_branch_loading_gt_1.05": int((load > 1.05).any(axis=1).sum()),
              "any_branch_loading_gt_1.30": int((load > 1.30).any(axis=1).sum()),
              "imbalance_below_-325.75MW": int((H["imbalance_MW"] < -325.75).sum())}
    counts[f"per_branch_loading_gt_{level:.2f}"] = {f"{net.bus_id[net.f[j]]}-{net.bus_id[net.t[j]]}": int(c) for j, c in enumerate(over.sum(axis=0)) if c}
    for key, v in expected.items():
        identical &= counts[key] == v

    # per-object statistics against phase-1 stage E, which was computed from the original chunks
    tail = pd.read_csv(STAGE_E / f"tail_statistics_{op_name}.csv")
    nr = np.setdiff1d(np.arange(net.n), [net.ref])
    blocks = {"Vm": np.abs(H["V"][:, nr]), "Pf": H["Pf_MW"], "Qf": H["Qf_Mvar"], "Pt": H["Pt_MW"], "Qt": H["Qt_Mvar"]}
    stat_diff = {}
    for key, X in blocks.items():
        s = ts.column_statistics(X.astype(np.float64))
        ref = tail[tail.quantity == key]
        stat_diff[key] = {f: float(np.nanmax(np.abs(s[f] - ref[f"{f}_heavy"].values))) for f in ("mean", "std", "excess_kurtosis", "q005", "q995")}
    worst = max(v for d in stat_diff.values() for v in d.values())
    identical &= worst < 1e-9

    chunks = sorted(STAGE_D.glob("*_chunk*.npz"))
    result = {"identical": bool(identical), "comparison": rows, "event_counts": counts, "expected_event_counts": expected, "max_abs_statistic_difference_vs_stage_e": stat_diff,
              "branch_100_103_index": k, "regenerated_chunks": [{"file": rr.relpath(p), "sha256": rr.sha256_file(p), "bytes": p.stat().st_size} for p in chunks],
              "committed_total_wall_time_s": old["total_wall_time_s"]}
    rr.write_json(out / "stage_g_identity.json", result)
    rr.write_manifest(
        out,
        command=f"python cloud/tests/run_stage_g_identity.py --result-stamp {a.result_stamp}" + (f" --op {op_name}" if op_name != "OP2" else ""),
        data_files=[mc.PGLIB118, STAGE_E / f"tail_statistics_{op_name}.csv"],
        parameters={"regeneration_command": f"python cloud/tests/run_stage_d_montecarlo.py {op_name} --result-stamp 20261021T202200 --workers 4",
                    "compared_summary_fields": ["scenarios", "converged", "failed_scenario_index", "iterations_mean", "iterations_max", "rounds_distribution"], "expected_event_counts": expected,
                    "loading": "max(|S_from|, |S_to|) / rateA from float32 check_heavy powers"},
        random_settings={"generator": "numpy default_rng (PCG64)", "seed_sequence_entropy": ENTROPY[op_name], "spawn": "as in phase-1 stage D"},
    )
    print("identical" if identical else "DIFFERENT", counts, worst)


if __name__ == "__main__":
    main()
