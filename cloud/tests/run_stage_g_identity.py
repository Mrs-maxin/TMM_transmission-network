"""Stage G: confirm that the regenerated OP2 scenario sets of phase-1 stage D are identical to the original run.

The OP2 chunks are regenerated beforehand with
    python cloud/tests/run_stage_d_montecarlo.py OP2 --result-stamp 20261021T202200 --workers 4
which overwrites stage_d_summary_OP2.json in the working tree. This entry compares that regenerated summary
with the committed one (git HEAD), recounts two phase-1 stage E events on check_heavy, records SHA256 of the
regenerated chunks, and copies the regenerated summary into its own result folder.

Usage: python cloud/tests/run_stage_g_identity.py --result-stamp YYYYmmddTHHMMSS
"""

import argparse
import json
import shutil
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
STAGE_D = rr.CLOUD_ROOT / "results" / "phase1_rts118" / "20261021T202200" / "stage_d_montecarlo" / "OP2"
# phase-1 stage E counted 153 scenarios with ANY branch above 1.05 rateA (graded_violation_probabilities_OP2.csv)
EXPECTED = {"any_branch_loading_gt_1.05": 153, "imbalance_below_-325.75MW": 228}
STAGE_E = rr.CLOUD_ROOT / "results" / "phase1_rts118" / "20261021T203828" / "stage_e_tail"
COMPARED = ["scenarios", "converged", "failures", "failed_scenario_index", "iterations_mean", "iterations_max",
            "rounds_distribution", "storage_precision", "input_model"]


def committed_summary(path):
    rel = rr.relpath(path)
    return json.loads(subprocess.check_output(["git", "show", f"HEAD:{rel}"], cwd=rr.PROJECT_ROOT, text=True))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--result-stamp", required=True)
    a = ap.parse_args()
    out = rr.stamped_stage_dir(TASK, a.result_stamp, "stage_g_identity")
    spath = STAGE_D / "stage_d_summary_OP2.json"
    new = json.loads(spath.read_text(encoding="utf-8"))
    old = committed_summary(spath)
    shutil.copy(spath, out / "regenerated_stage_d_summary_OP2.json")

    rows, identical = [], True
    for set_name, info in old["sets"].items():
        for key in COMPARED:
            same = info[key] == new["sets"][set_name][key]
            identical &= same
            rows.append({"set": set_name, "quantity": key, "committed": info[key], "regenerated": new["sets"][set_name][key], "identical": same})

    net = nw.from_matpower(mc.PGLIB118)
    H = mcarlo.load_set(sorted(STAGE_D.glob("check_heavy_chunk*.npz")), ["V", "Pf_MW", "Qf_Mvar", "Pt_MW", "Qt_Mvar", "imbalance_MW", "converged"])
    k = int(np.flatnonzero((net.bus_id[net.f] == 100) & (net.bus_id[net.t] == 103))[0])
    load = np.maximum(np.hypot(H["Pf_MW"], H["Qf_Mvar"]), np.hypot(H["Pt_MW"], H["Qt_Mvar"])) / net.rateA[None]
    over = load > 1.05
    counts = {"branch_100_103_loading_gt_1.05": int(over[:, k].sum()),
              "any_branch_loading_gt_1.05": int(over.any(axis=1).sum()),
              "imbalance_below_-325.75MW": int((H["imbalance_MW"] < -325.75).sum())}
    counts["per_branch_loading_gt_1.05"] = {f"{net.bus_id[net.f[j]]}-{net.bus_id[net.t[j]]}": int(c) for j, c in enumerate(over.sum(axis=0)) if c}
    identical &= counts["any_branch_loading_gt_1.05"] == EXPECTED["any_branch_loading_gt_1.05"]
    identical &= counts["imbalance_below_-325.75MW"] == EXPECTED["imbalance_below_-325.75MW"]

    # per-object statistics against phase-1 stage E, which was computed from the original chunks
    tail = pd.read_csv(STAGE_E / "tail_statistics_OP2.csv")
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
    result = {"identical": bool(identical), "comparison": rows, "event_counts": counts, "expected_event_counts": EXPECTED, "max_abs_statistic_difference_vs_stage_e": stat_diff,
              "branch_100_103_index": k, "regenerated_chunks": [{"file": rr.relpath(p), "sha256": rr.sha256_file(p), "bytes": p.stat().st_size} for p in chunks],
              "regenerated_total_wall_time_s": new["total_wall_time_s"], "committed_total_wall_time_s": old["total_wall_time_s"]}
    rr.write_json(out / "stage_g_identity.json", result)
    rr.write_manifest(
        out,
        command=f"python cloud/tests/run_stage_g_identity.py --result-stamp {a.result_stamp}",
        data_files=[mc.PGLIB118, STAGE_E / "tail_statistics_OP2.csv"],
        parameters={"regeneration_command": "python cloud/tests/run_stage_d_montecarlo.py OP2 --result-stamp 20261021T202200 --workers 4",
                    "compared_summary_fields": COMPARED, "expected_event_counts": EXPECTED,
                    "loading": "max(|S_from|, |S_to|) / rateA from float32 check_heavy powers"},
        random_settings={"generator": "numpy default_rng (PCG64)", "seed_sequence_entropy": 2026102112, "spawn": "as in phase-1 stage D"},
    )
    print("identical" if identical else "DIFFERENT", counts, worst)


if __name__ == "__main__":
    main()
