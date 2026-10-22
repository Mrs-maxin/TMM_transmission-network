"""Stage J: violation probabilities of OP2 from the AC Monte Carlo check sets, the modeling scenarios and 2^20
stratified Sobol points of each selected state model (part 'table'), and their repeatability over R new
independent modeling (5256) + determination (1752) pairs refitted with their own K selection (part 'repeat').

Usage: python cloud/tests/run_stage_j_violations.py --result-stamp YYYYmmddTHHMMSS --part table [--op OP1]
       python cloud/tests/run_stage_j_violations.py --result-stamp YYYYmmddTHHMMSS --part repeat --replications 10
"""

import os

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from tmm_cloud import error_model as em
from tmm_cloud import matpower_case as mc
from tmm_cloud import montecarlo as mcarlo
from tmm_cloud import network as nw
from tmm_cloud import run_record as rr
from tmm_cloud import state_model as sm
from tmm_cloud import violation_events as ve

TASK = "phase1b_tmm_op2"
OP = "OP2"
N_POINTS = 1 << 20
SOBOL_ENTROPY = 2026102201
REPEAT_ENTROPY = 2026102200
REPEAT_SOBOL_ENTROPY = 2026102300
REPEAT_EVENTS = ["a_branch_100_103_above_1.05", "b_any_unit_below_Pmin", "d_branch_77_82_above_1.00",
                 "f_any_branch_above_1.05", "f_any_branch_above_1.10"]
METHODS = ["GMM", "TMM", "block TMM", "two-strength block TMM"]
LOAD_CV = 0.01
CHUNK = 5000
STAGE_D_ROOT = rr.CLOUD_ROOT / "results" / "phase1_rts118" / "20261021T202200" / "stage_d_montecarlo"
STAGE_D = STAGE_D_ROOT / OP
ANY_LEVELS = {"OP2": ve.ANY_BRANCH_LEVELS, "OP1": (1.05, 1.10, 1.20, 1.30)}
STAGE_B = rr.CLOUD_ROOT / "results" / "phase1_rts118" / "20261021T201240" / "stage_b_network"
STAGE_C = rr.CLOUD_ROOT / "results" / "phase1_rts118" / "20261021T201850" / "stage_c_error_model"


def load_models(hdir):
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import run_stage_i_marginals as si
    si.OP = OP
    return si.load_stage_h(hdir)


def counted_rows(label, q, thr):
    k, kc, n = ve.counted(q, thr, ANY_LEVELS[OP])
    rows = []
    for e, ke in zip(ve.event_names(ANY_LEVELS[OP]), k):
        lo, hi = ve.clopper_pearson(int(ke), n)
        rows.append({"source": label, "event": e, "scenarios": n, "count": int(ke), "probability": ke / n, "ci95_low": lo, "ci95_high": hi})
    return rows, kc / n


def part_table(a, net, op, outmap, thr, hdir):
    events = ve.event_names(ANY_LEVELS[OP])
    out = rr.stamped_stage_dir(TASK, a.result_stamp, "stage_j_violations" if OP == "OP2" else f"stage_k_{OP}/stage_j_violations")
    t0 = time.time()
    coords, models = load_models(hdir)
    check_files = sorted(STAGE_D.glob("check_heavy_chunk*.npz"))
    light_files = sorted(STAGE_D.glob("check_light_chunk*.npz"))
    mod_files = sorted(STAGE_D.glob("modeling_chunk*.npz"))
    dc = mcarlo.load_set(check_files)
    check = outmap.from_stage_d(dc)
    rows, curves = [], {}
    r, curves["AC Monte Carlo check (100000)"] = counted_rows("check_heavy", check, thr)
    rows += r
    ref = {x["event"]: x for x in r}

    # equivalence of the unit-limit events with the imbalance thresholds on the check set
    below, above = ve.unit_limit_check(check, thr)
    ev, _ = ve.indicators(check, thr, ANY_LEVELS[OP])
    eq = {"below_Pmin_direct": int(below.sum()), "below_Pmin_by_lam": int(ev[:, 1].sum()),
          "below_Pmin_disagreements": int((below != ev[:, 1]).sum()),
          "above_Pmax_direct": int(above.sum()), "above_Pmax_by_lam": int(ev[:, 2].sum()),
          "above_Pmax_disagreements": int((above != ev[:, 2]).sum()),
          "max_abs_Pg_minus_Pg0_plus_alpha_lam_MW": float(np.nanmax(np.abs(
              check["Pg"] - (op.Pg0[net.Pmax > 0][None] * net.baseMVA + check["lam"] * op.alpha[net.Pmax > 0][None]))))}
    del dc, check

    dl = mcarlo.load_set(light_files)
    light = outmap.from_stage_d(dl)
    r, curves["light-tailed check (20000)"] = counted_rows("check_light", light, thr)
    rows += r
    del dl, light
    _, dmod = sm.scenario_state(mod_files, net)
    r, curves["modeling scenarios (5256)"] = counted_rows("modeling", outmap.from_stage_d(dmod), thr)
    rows += r

    timing, comp = {}, {}
    for i, method in enumerate(METHODS):
        if method not in models:
            continue
        ts = time.time()
        res = ve.model_probabilities(models[method][0], coords, outmap, thr, N_POINTS, SOBOL_ENTROPY + i, any_levels=ANY_LEVELS[OP])
        timing[method] = time.time() - ts
        curves[method] = res["curve"]
        comp[method] = {"allocation": res["allocation"].tolist(), "per_component": res["per_component"].tolist()}
        for e, p, h in zip(events, res["events"], res["point_hits"]):
            rows.append({"source": method, "event": e, "scenarios": N_POINTS, "count": int(h), "probability": float(p),
                         "ci95_low": np.nan, "ci95_high": np.nan})
        print(method, f"{timing[method]:.0f} s", flush=True)

    df = pd.DataFrame(rows)
    df["check_probability"] = df.event.map({e: x["probability"] for e, x in ref.items()})
    df["relative_error_to_check"] = df.probability / df.check_probability - 1
    df["inside_check_ci95"] = (df.probability >= df.event.map({e: x["ci95_low"] for e, x in ref.items()})) & \
                              (df.probability <= df.event.map({e: x["ci95_high"] for e, x in ref.items()}))
    df.to_csv(out / f"violation_probabilities_{OP}.csv", index=False)
    wide = df.pivot(index="event", columns="source", values="probability").reindex(events)
    wide.to_csv(out / f"violation_probabilities_wide_{OP}.csv")
    cdf = pd.DataFrame({"loading_threshold": ve.CURVE_LEVELS, **{k: v for k, v in curves.items()}})
    cdf.to_csv(out / f"branch_100_103_exceedance_curve_{OP}.csv", index=False)

    kc = np.rint(curves["AC Monte Carlo check (100000)"] * 100000).astype(int)
    fig, ax = plt.subplots(figsize=(6, 4))
    lo = [ve.clopper_pearson(int(k), 100000)[0] for k in kc]
    hi = [ve.clopper_pearson(int(k), 100000)[1] for k in kc]
    ax.fill_between(ve.CURVE_LEVELS, np.clip(lo, 1e-7, None), hi, color="0.85", label="check 95% Clopper-Pearson")
    styles = {"AC Monte Carlo check (100000)": ("k", "-"), "light-tailed check (20000)": ("0.4", ":"),
              "modeling scenarios (5256)": ("0.55", "--"), "GMM": ("tab:blue", "--"), "TMM": ("tab:red", "-"),
              "block TMM": ("tab:green", "-."), "two-strength block TMM": ("tab:purple", ":")}
    for k, v in curves.items():
        ax.semilogy(ve.CURVE_LEVELS, np.where(v > 0, v, np.nan), color=styles[k][0], ls=styles[k][1], label=k)
    ax.set_xlabel("loading threshold of branch 100-103 (fraction of rateA)")
    ax.set_ylabel("P(loading > threshold)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out / f"branch_100_103_exceedance_curve_{OP}.png", dpi=130)
    plt.close(fig)

    thr_out = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in thr.items()}
    rr.write_json(out / f"stage_j_summary_{OP}.json", {"operating_point": OP, "thresholds": thr_out, "unit_limit_equivalence": eq,
                                                     "model_components": comp, "propagation_seconds": timing, "wall_time_s": time.time() - t0})
    rr.write_manifest(
        out,
        command=f"python cloud/tests/run_stage_j_violations.py --result-stamp {a.result_stamp} --part table" + (f" --op {OP}" if OP != "OP2" else ""),
        data_files=[mc.PGLIB118, STAGE_B / f"operating_point_{OP}.npz", hdir / f"coordinates_{OP}.npz", *[p for _, p in models.values()], *mod_files],
        parameters={"events": events, "loading": "max(|S_from|, |S_to|)/rateA", "curve_levels": ve.CURVE_LEVELS.tolist(),
                    "unit_limits": "units with alpha > 0; below Pmin iff lam < max_i (Pmin_i - Pg0_i)/alpha_i, above Pmax iff lam > min_i (Pmax_i - Pg0_i)/alpha_i",
                    "voltage_events": "non-reference buses (reference bus voltage fixed at its setpoint)",
                    "check_sets": f"check_heavy 100000 and check_light 20000 of phase-1 stage D {OP}, chunk SHA256 in stage_g_identity.json",
                    "interval": "Clopper-Pearson 95%", "points_per_model": N_POINTS},
        random_settings={"sobol": "same points as stage I", "entropy": {m: SOBOL_ENTROPY + i for i, m in enumerate(METHODS)}},
    )
    print(rr.relpath(out), f"{time.time() - t0:.0f} s")


def part_repeat(a, net, op, outmap, thr):
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from run_stage_h_fit import fit_all
    out = rr.stamped_stage_dir(TASK, a.result_stamp, "stage_j_repeatability")
    t0 = time.time()
    op_path = STAGE_B / f"operating_point_{OP}.npz"
    model_path = STAGE_C / f"error_model_{OP}.npz"
    emodel = em.from_npz(np.load(model_path))
    mapping = pd.read_csv(STAGE_B / "renewable_mapping.csv")
    idx = [ve.EVENTS.index(e) for e in REPEAT_EVENTS]
    rows_path = out / f"replication_probabilities_{OP}.csv"
    rows = pd.read_csv(rows_path).to_dict("records") if rows_path.exists() else []
    done = {int(r["replication"]) for r in rows}
    info_path = out / f"replication_info_{OP}.json"
    info = json.loads(info_path.read_text()) if info_path.exists() else {}
    for rep in range(1, a.replications + 1):
        if rep in done:
            continue
        tr = time.time()
        rdir = out / f"rep{rep:02d}"
        rdir.mkdir(exist_ok=True)
        children = np.random.SeedSequence(REPEAT_ENTROPY + rep).spawn(4)
        pf_info = {}
        for (set_name, n), child in zip([("modeling", 5256), ("determination", 1752)], children[:2]):
            rng = np.random.default_rng(child)
            e, _ = em.sample_errors(emodel, n, rng, "heavy")
            Pren = (op.Pren[None, :] * net.baseMVA + e) / net.baseMVA
            kload = 1.0 + LOAD_CV * rng.standard_normal((n, len(op.Pd)))
            r = mcarlo.run_set(op_path, Pren, kload, a.workers, CHUNK, rdir, set_name, {}, single=False)
            pf_info[set_name] = {k: v for k, v in r.items() if k != "files"}
        X, dmod = sm.scenario_state(sorted(rdir.glob("modeling_chunk*.npz")), net, drop_failed=True)
        Xv, ddet = sm.scenario_state(sorted(rdir.glob("determination_chunk*.npz")), net, drop_failed=True)
        tf = time.time()
        coords, dirs, dinfo, best, selected, cand, secs = fit_all(net, op, mapping, X, Xv, dmod, a.workers, verbose=False)
        fit_s = time.time() - tf
        status = pd.DataFrame(cand).groupby(["method", "status"]).size().reset_index().values.tolist() if cand else []
        qm = outmap.from_stage_d(dmod)
        qd = outmap.from_stage_d(ddet)
        for label, qs in [("counting modeling (5256)", [qm]), ("counting modeling + determination (7008)", [qm, qd])]:
            k = sum(ve.counted(q, thr)[0] for q in qs)
            n = sum(len(q["lam"]) for q in qs)
            for e, j in zip(REPEAT_EVENTS, idx):
                rows.append({"replication": rep, "source": label, "K": np.nan, "event": e, "probability": k[j] / n, "count": int(k[j]), "scenarios": n})
        prop_s = {}
        for i, method in enumerate(METHODS):
            if method not in selected:
                continue
            ts = time.time()
            res = ve.model_probabilities(best[(method, selected[method])], coords, outmap, thr, N_POINTS,
                                         REPEAT_SOBOL_ENTROPY + 10 * rep + i)
            prop_s[method] = time.time() - ts
            for e, j in zip(REPEAT_EVENTS, idx):
                rows.append({"replication": rep, "source": method, "K": selected[method], "event": e,
                             "probability": float(res["events"][j]), "count": int(res["point_hits"][j]), "scenarios": N_POINTS})
        info[str(rep)] = {"entropy": REPEAT_ENTROPY + rep, "power_flow": pf_info, "d": int(coords.d), "r": dinfo["r"],
                          "selected_K": selected, "candidate_status_counts": status,
                          "df": {m: np.round(best[(m, K)].df, 3).tolist() for m, K in selected.items() if m != "GMM"},
                          "fit_seconds": fit_s, "propagation_seconds": prop_s, "replication_seconds": time.time() - tr}
        pd.DataFrame(rows).to_csv(rows_path, index=False)
        rr.write_json(info_path, info)
        print(f"replication {rep}: K* {selected}, {time.time() - tr:.0f} s", flush=True)

    df = pd.DataFrame(rows)
    df = df[df.replication <= a.replications]
    check = mcarlo.load_set(sorted(STAGE_D.glob("check_heavy_chunk*.npz")))
    kc, _, nc = ve.counted(outmap.from_stage_d(check), thr)
    pref = {e: kc[j] / nc for e, j in zip(REPEAT_EVENTS, idx)}
    summ = []
    for (src, e), g in df.groupby(["source", "event"], sort=False):
        p = g.probability.to_numpy()
        rel = p / pref[e] - 1
        summ.append({"source": src, "event": e, "replications": len(p), "check_probability": pref[e],
                     "median": float(np.median(p)), "p10": float(np.quantile(p, 0.1)), "p90": float(np.quantile(p, 0.9)),
                     "min": float(p.min()), "max": float(p.max()),
                     "median_relative_error": float(np.median(rel)), "median_abs_relative_error": float(np.median(np.abs(rel))),
                     "rms_relative_error": float(np.sqrt(np.mean(rel ** 2))), "zero_count_replications": int((p == 0).sum())})
    sdf = pd.DataFrame(summ)
    sdf.to_csv(out / f"repeatability_summary_{OP}.csv", index=False)

    srcs = ["counting modeling (5256)", "counting modeling + determination (7008)", *METHODS]
    for e in REPEAT_EVENTS:
        fig, ax = plt.subplots(figsize=(6.5, 4))
        g = df[df.event == e]
        for i, s in enumerate(srcs):
            v = g[g.source == s].probability.to_numpy()
            ax.scatter(np.full(len(v), i) + np.linspace(-0.15, 0.15, len(v)), v, s=14)
        ax.axhline(pref[e], color="k", lw=1, label="check (100000)")
        ax.set_xticks(range(len(srcs)))
        ax.set_xticklabels(["count 5256", "count 7008", "GMM", "TMM", "block TMM", "two-strength"], fontsize=8)
        ax.set_ylabel("probability")
        ax.set_title(e, fontsize=9)
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(out / f"repeatability_{e}_{OP}.png", dpi=130)
        plt.close(fig)

    rr.write_manifest(
        out,
        command=f"python cloud/tests/run_stage_j_violations.py --result-stamp {a.result_stamp} --part repeat --replications {a.replications} --workers {a.workers}",
        data_files=[mc.PGLIB118, op_path, model_path, STAGE_B / "renewable_mapping.csv"],
        parameters={"replications": a.replications, "sets_per_replication": {"modeling": 5256, "determination": 1752},
                    "solver": "phase-1 tmm_cloud.powerflow.BatchSolver through montecarlo.run_set, float64 storage",
                    "fitting": "stage H settings, every method with its own K selection by eq. (22)", "events": REPEAT_EVENTS,
                    "points_per_model": N_POINTS, "reference": "check_heavy (100000) of phase-1 stage D OP2",
                    "counting": "empirical share of the 5256 modeling scenarios, and of the 7008 modeling + determination scenarios (equal power-flow count)"},
        random_settings={"scenario_entropy": f"{REPEAT_ENTROPY} + r, SeedSequence spawn(4), children 0 and 1 for modeling and determination as in stage D",
                         "sobol_entropy": f"{REPEAT_SOBOL_ENTROPY} + 10 r + method index"},
    )
    print(rr.relpath(out), f"{time.time() - t0:.0f} s")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--result-stamp", required=True)
    ap.add_argument("--part", choices=["table", "repeat"], required=True)
    ap.add_argument("--replications", type=int, default=10)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--op", default="OP2", choices=["OP1", "OP2"])
    a = ap.parse_args()
    global OP, STAGE_D
    OP, STAGE_D = a.op, STAGE_D_ROOT / a.op
    if a.part == "repeat" and OP != "OP2":
        raise SystemExit("the repeatability study is defined for OP2 only")
    net = nw.from_matpower(mc.PGLIB118)
    op = mcarlo.load_operating_point(STAGE_B / f"operating_point_{OP}.npz")
    outmap = sm.OutputMap(net, op)
    thr = ve.thresholds(net, op)
    hdir = rr.CLOUD_ROOT / "results" / TASK / a.result_stamp / ("stage_h_fit" if OP == "OP2" else f"stage_k_{OP}/stage_h_fit")
    if a.part == "table":
        part_table(a, net, op, outmap, thr, hdir)
    else:
        part_repeat(a, net, op, outmap, thr)


if __name__ == "__main__":
    main()
