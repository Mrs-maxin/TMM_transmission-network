"""Stage D: scenario sets and AC Monte Carlo for one operating point.

Usage: python cloud/tests/run_stage_d_montecarlo.py OP1 [--result-stamp YYYYmmddTHHMMSS] [--workers 4]
"""

import argparse
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from tmm_cloud import CLOUD_ROOT
from tmm_cloud import error_model as em
from tmm_cloud import matpower_case as mc
from tmm_cloud import montecarlo as mcarlo
from tmm_cloud import run_record as rr

TASK = "phase1_rts118"
SETS = [("modeling", 5256, "heavy"), ("determination", 1752, "heavy"), ("check_heavy", 100000, "heavy"), ("check_light", 20000, "gauss")]
RNG_BASE = {"OP1": 2026102111, "OP2": 2026102112}
LOAD_CV = 0.01
CHUNK = 5000


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("op")
    ap.add_argument("--result-stamp")
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    name = a.op
    if a.result_stamp:
        out = CLOUD_ROOT / "results" / TASK / a.result_stamp / "stage_d_montecarlo" / name
        out.mkdir(parents=True, exist_ok=True)
    else:
        out = rr.new_result_dir(TASK, f"stage_d_montecarlo/{name}")
    bdir = rr.latest_stage_dir(TASK, "stage_b_network")
    cdir = rr.latest_stage_dir(TASK, "stage_c_error_model")
    op_path = bdir / f"operating_point_{name}.npz"
    model_path = cdir / f"error_model_{name}.npz"
    model = em.from_npz(np.load(model_path))
    op = mcarlo.load_operating_point(op_path)
    base = 100.0
    f_op = op.Pren * base
    n_bus = len(op.Pd)

    ss = np.random.SeedSequence(RNG_BASE[name])
    children = ss.spawn(len(SETS))
    t0 = time.time()
    summary = {"operating_point": name, "heavy_tailed_model": model.heavy_kind, "nu_groups": model.nu_groups.tolist(), "sets": {}}
    for (set_name, n, kind), child in zip(SETS, children):
        rng = np.random.default_rng(child)
        ts = time.time()
        e, w = em.sample_errors(model, n, rng, kind)
        Pren = (f_op[None, :] + e) / base
        kload = 1.0 + LOAD_CV * rng.standard_normal((n, n_bus))
        t_sample = time.time() - ts
        info = mcarlo.run_set(op_path, Pren, kload, a.workers, CHUNK, out, set_name,
                              {"strength_w": w.astype(np.float32)}, single=set_name.startswith("check"))
        info["storage_precision"] = "float32" if set_name.startswith("check") else "float64"
        info["sampling_time_s"] = t_sample
        info["input_model"] = model.heavy_kind if kind == "heavy" else "gauss"
        info["files"] = [rr.relpath(p) for p in info["files"]]
        info["disk_MB"] = float(sum((rr.PROJECT_ROOT / p).stat().st_size for p in info["files"]) / 2 ** 20)
        summary["sets"][set_name] = info
        print(set_name, {k: v for k, v in info.items() if k not in ("files", "failed_scenario_index")}, flush=True)
    summary["total_wall_time_s"] = time.time() - t0
    rr.write_json(out / f"stage_d_summary_{name}.json", summary)
    rr.write_manifest(
        out,
        command=f"python cloud/tests/run_stage_d_montecarlo.py {name}" + (f" --result-stamp {a.result_stamp}" if a.result_stamp else "") + f" --workers {a.workers}",
        data_files=[mc.PGLIB118, op_path, model_path],
        parameters={
            "sets": [{"name": s, "scenarios": n, "input_model": k} for s, n, k in SETS],
            "renewable_injection": "day-ahead forecast + sampled error, within [0, capacity], unity power factor",
            "load": "bus P and Q multiplied by 1 + 0.01 xi, xi independent standard normal",
            "solver": "tmm_cloud.powerflow.BatchSolver: rectangular Newton-Raphson, distributed imbalance by Pmax, Q limits (all violators switch together, no return), tolerance 1e-10 p.u., max 30 iterations per round, max 30 rounds",
            "storage": "chunks of 5000 scenarios, numpy savez_compressed; check sets: V complex64, unit and branch powers and load factors float32; modeling and determination sets: float64; imbalance and renewable injections float64; NaN for non-converged scenarios",
            "stage_b_result": rr.relpath(bdir), "stage_c_result": rr.relpath(cdir),
        },
        random_settings={"generator": "numpy default_rng (PCG64)", "seed_sequence_entropy": RNG_BASE[name],
                         "spawn": "SeedSequence(entropy).spawn(4) in order modeling, determination, check_heavy, check_light; per set: copula Gaussian vector, chi-square strengths, then load normals"},
    )
    print(rr.relpath(out), f"{summary['total_wall_time_s']:.0f} s")


if __name__ == "__main__":
    main()
