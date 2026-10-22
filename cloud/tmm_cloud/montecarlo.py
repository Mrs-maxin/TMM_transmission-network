"""Parallel AC Monte Carlo with the batched Newton-Raphson solver and chunked compressed storage."""

import os
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from . import matpower_case as mc
from . import network as nw
from . import powerflow as pf

_STATE = {}
SUB_BATCH = 200


def load_operating_point(path):
    z = np.load(path)
    op = pf.OperatingPoint(**{k: z[k] for k in ["Pg0", "Vset", "alpha", "Pd", "Qd", "ren_bus", "Pren", "ren_cap"]})
    return op


def _init(op_path):
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    net = nw.from_matpower(mc.PGLIB118)
    op = load_operating_point(op_path)
    _STATE["net"] = net
    _STATE["op"] = op
    _STATE["solver"] = pf.BatchSolver(net, op)


def _solve(args):
    Pren, kload = args
    net, op, solver = _STATE["net"], _STATE["op"], _STATE["solver"]
    out = []
    t0 = time.time()
    for i in range(0, len(Pren), SUB_BATCH):
        pr = Pren[i:i + SUB_BATCH]
        k = kload[i:i + SUB_BATCH]
        out.append(solver.solve(op.Pd[None, :] * k, op.Qd[None, :] * k, pr))
    res = {key: np.concatenate([o[key] for o in out]) for key in out[0]}
    res["elapsed_s"] = time.time() - t0
    return res


def pack(res, single):
    """Convert solver output to storage arrays; single=True stores voltages and powers in float32."""
    ok = res["fail"] == pf.FAIL_OK
    fl, cx = (np.float32, np.complex64) if single else (np.float64, np.complex128)
    Sf = res["Sf_MVA"]
    St = res["St_MVA"]
    d = {
        "V": res["V"].astype(cx),
        "Pg_MW": res["Pg_MW"].astype(fl),
        "Qg_Mvar": res["Qg_Mvar"].astype(fl),
        "mode": res["mode"].astype(np.int8),
        "Pf_MW": Sf.real.astype(fl), "Qf_Mvar": Sf.imag.astype(fl),
        "Pt_MW": St.real.astype(fl), "Qt_Mvar": St.imag.astype(fl),
        "imbalance_MW": res["lam_MW"].astype(np.float64),
        "converged": ok,
        "fail_code": res["fail"].astype(np.int8),
        "iterations": res["iterations"].astype(np.int16),
        "rounds": res["rounds"].astype(np.int8),
        "max_P_residual_pu": res["max_P_residual_pu"].astype(np.float32),
    }
    for key in ("V", "Pg_MW", "Qg_Mvar", "Pf_MW", "Qf_Mvar", "Pt_MW", "Qt_Mvar"):
        d[key][~ok] = np.nan
    return d


def run_set(op_path, Pren, kload, n_workers, chunk_size, out_dir, prefix, extra_inputs, single=True):
    """Solve all scenarios of one set and write chunk files <prefix>_chunkNNN.npz. Returns timing and failures."""
    net = nw.from_matpower(mc.PGLIB118)
    op = load_operating_point(op_path)
    solver = pf.BatchSolver(net, op)
    n = len(Pren)
    starts = list(range(0, n, chunk_size))
    per_worker = max(SUB_BATCH, int(np.ceil(chunk_size / n_workers / SUB_BATCH)) * SUB_BATCH)
    t0 = time.time()
    files, fail_codes, iters, rounds, worker_time = [], [], [], [], 0.0
    with ProcessPoolExecutor(max_workers=n_workers, initializer=_init, initargs=(str(op_path),)) as ex:
        for c, s in enumerate(starts):
            e = min(n, s + chunk_size)
            parts = [(Pren[a:min(e, a + per_worker)], kload[a:min(e, a + per_worker)]) for a in range(s, e, per_worker)]
            results = list(ex.map(_solve, parts))
            worker_time += sum(r.pop("elapsed_s") for r in results)
            res = {key: np.concatenate([r[key] for r in results]) for key in results[0]}
            d = pack(res, single)
            d["scenario_index"] = np.arange(s, e, dtype=np.int32)
            d["Pren_MW"] = (Pren[s:e] * net.baseMVA).astype(np.float64)
            d["load_factor"] = kload[s:e].astype(np.float32 if single else np.float64)
            for key, arr in extra_inputs.items():
                d[key] = arr[s:e]
            path = out_dir / f"{prefix}_chunk{c:03d}.npz"
            np.savez_compressed(path, **d)
            files.append(path)
            fail_codes.append(d["fail_code"])
            iters.append(d["iterations"])
            rounds.append(d["rounds"])
            print(f"{prefix}: {e}/{n} scenarios, {time.time() - t0:.1f} s", flush=True)
    fc = np.concatenate(fail_codes)
    return {
        "scenarios": int(n),
        "converged": int((fc == pf.FAIL_OK).sum()),
        "failures": {pf.FAIL_REASONS[int(k)]: int((fc == k).sum()) for k in np.unique(fc) if k != pf.FAIL_OK},
        "failed_scenario_index": np.flatnonzero(fc != pf.FAIL_OK).tolist(),
        "wall_time_s": time.time() - t0,
        "summed_worker_time_s": worker_time,
        "workers": n_workers,
        "iterations_mean": float(np.concatenate(iters).mean()),
        "iterations_max": int(np.concatenate(iters).max()),
        "rounds_distribution": {int(k): int(v) for k, v in zip(*np.unique(np.concatenate(rounds), return_counts=True))},
        "files": files,
    }


def load_set(paths, keys=None):
    parts = [np.load(p) for p in sorted(paths)]
    keys = keys or parts[0].files
    return {k: np.concatenate([p[k] for p in parts]) for k in keys}
