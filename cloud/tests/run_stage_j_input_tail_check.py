"""Stage J diagnostic: upper tail of the total renewable injection error of the OP2 heavy check set (100000)
against 10^6 fresh samples of the same input model. The system imbalance is almost a linear function of this
total (correlation reported), so the comparison shows how far the check-set reference of the unit-limit event
is from the input model by sampling variation alone, without power flows.

Usage: python cloud/tests/run_stage_j_input_tail_check.py --result-stamp YYYYmmddTHHMMSS
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from tmm_cloud import error_model as em
from tmm_cloud import montecarlo as mcarlo
from tmm_cloud import run_record as rr
from tmm_cloud import violation_events as ve

TASK = "phase1b_tmm_op2"
OP = "OP2"
FRESH_ENTROPY = 2026102400
FRESH_BLOCKS, BLOCK = 10, 100000
LEVELS = (0.99, 0.995, 0.9977, 0.999)
STAGE_D = rr.CLOUD_ROOT / "results" / "phase1_rts118" / "20261021T202200" / "stage_d_montecarlo" / OP
STAGE_B = rr.CLOUD_ROOT / "results" / "phase1_rts118" / "20261021T201240" / "stage_b_network"
STAGE_C = rr.CLOUD_ROOT / "results" / "phase1_rts118" / "20261021T201850" / "stage_c_error_model"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--result-stamp", required=True)
    a = ap.parse_args()
    out = rr.stamped_stage_dir(TASK, a.result_stamp, "stage_j_repeatability/input_tail_check")
    model = em.from_npz(np.load(STAGE_C / f"error_model_{OP}.npz"))
    op = mcarlo.load_operating_point(STAGE_B / f"operating_point_{OP}.npz")
    files = sorted(STAGE_D.glob("check_heavy_chunk*.npz"))
    chk = mcarlo.load_set(files, ["Pren_MW", "imbalance_MW"])
    forecast = op.Pren.sum() * 100.0
    E = chk["Pren_MW"].sum(axis=1) - forecast
    corr = float(np.corrcoef(chk["imbalance_MW"], E)[0, 1])
    fresh = []
    for b in range(FRESH_BLOCKS):
        e, _ = em.sample_errors(model, BLOCK, np.random.default_rng(np.random.SeedSequence(FRESH_ENTROPY + b)), "heavy")
        fresh.append(e.sum(axis=1))
    Ef = np.concatenate(fresh)
    rows = []
    for p in LEVELS:
        q = float(np.quantile(E, p))
        kc, kf = int((E > q).sum()), int((Ef > q).sum())
        pc, pf = kc / len(E), kf / len(Ef)
        se = np.sqrt(pc * (1 - pc) / len(E) + pf * (1 - pf) / len(Ef))
        rows.append({"check_quantile_level": p, "threshold_MW": q, "check_count": kc, "check_probability": pc,
                     "check_ci95": ve.clopper_pearson(kc, len(E)), "fresh_count": kf, "fresh_probability": pf,
                     "fresh_ci95": ve.clopper_pearson(kf, len(Ef)), "difference_in_standard_errors": float((pc - pf) / se)})
    pd.DataFrame(rows).to_csv(out / f"input_tail_check_{OP}.csv", index=False)
    rr.write_json(out / f"input_tail_check_{OP}.json", {"correlation_imbalance_vs_total_error": corr, "rows": rows,
                                                       "check_imbalance_below_-325.75MW": int((chk["imbalance_MW"] < -325.75).sum())})
    rr.write_manifest(
        out,
        command=f"python cloud/tests/run_stage_j_input_tail_check.py --result-stamp {a.result_stamp}",
        data_files=[STAGE_C / f"error_model_{OP}.npz", STAGE_B / f"operating_point_{OP}.npz"],
        parameters={"quantity": "total renewable injection minus total forecast (MW), within [0, capacity] per station",
                    "thresholds": "check-set quantiles at the listed levels", "levels": LEVELS,
                    "fresh_samples": FRESH_BLOCKS * BLOCK, "check_set": "phase-1 stage D OP2 check_heavy (chunk SHA256 in stage_g_identity.json)"},
        random_settings={"fresh": f"numpy default_rng(SeedSequence({FRESH_ENTROPY} + b)), b = 0..{FRESH_BLOCKS - 1}, {BLOCK} heavy-tailed error vectors each"},
    )
    print(pd.DataFrame(rows).to_string(), corr)


if __name__ == "__main__":
    main()
