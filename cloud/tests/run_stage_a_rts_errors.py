"""Stage A: RTS-GMLC 2020 wind, UPV and regional load forecast errors (real-time hourly mean minus day-ahead)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from tmm_cloud import rts_gmlc as rts
from tmm_cloud import run_record as rr

TASK = "phase1_rts118"


def offdiag_mean(c):
    a = c.values
    m = ~np.eye(len(a), dtype=bool) & np.isfinite(a)
    return float(a[m].mean())


def main():
    out = rr.new_result_dir(TASK, "stage_a_rts_errors")
    s = rts.load_series()
    cap = s["capacity"]

    wind_err = s["wind_rt"] - s["wind_da"]
    day = rts.daylight_mask(s["pv_da"], s["pv_rt"])
    pv_err = (s["pv_rt"] - s["pv_da"]).where(day)
    load_err = s["load_rt"] - s["load_da"]

    wind_stats = rts.site_statistics(wind_err, s["wind_da"], cap)
    pv_stats = rts.site_statistics(pv_err, s["pv_da"], cap)
    wind_stats.to_csv(out / "wind_error_site_statistics.csv", index=False)
    pv_stats.to_csv(out / "pv_error_site_statistics.csv", index=False)

    wc, wc2 = rts.pairwise_corr(wind_err)
    pc, pc2 = rts.pairwise_corr(pv_err)
    wc.to_csv(out / "wind_error_corr.csv")
    wc2.to_csv(out / "wind_sq_std_error_corr.csv")
    pc.to_csv(out / "pv_error_corr.csv")
    pc2.to_csv(out / "pv_sq_std_error_corr.csv")

    wz = rts.conditional_standardise(wind_err, s["wind_da"], cap)
    pz = rts.conditional_standardise(pv_err, s["pv_da"], cap)
    wcz2 = (wz ** 2).corr(min_periods=200)
    pcz2 = (pz ** 2).corr(min_periods=200)
    wcz2.to_csv(out / "wind_sq_condstd_error_corr.csv")
    pcz2.to_csv(out / "pv_sq_condstd_error_corr.csv")
    wind_cond_kurt = {c: rts.excess_kurtosis(wz[c].dropna()) for c in wz.columns}
    pv_cond_kurt = {c: rts.excess_kurtosis(pz[c].dropna()) for c in pz.columns}

    agg_wind = wind_err.sum(axis=1)
    both = pd.concat([agg_wind.rename("wind_total"), pv_err.sum(axis=1, min_count=1).rename("pv_total")], axis=1).dropna()
    wind_pv_corr = float(both.corr().iloc[0, 1])

    spread = pd.concat([
        rts.spread_by_forecast_level(wind_err, s["wind_da"], cap).assign(kind="wind"),
        rts.spread_by_forecast_level(pv_err, s["pv_da"], cap).assign(kind="pv"),
    ])
    spread.to_csv(out / "error_spread_by_forecast_level.csv", index=False)

    load_rows = []
    for r in load_err.columns:
        e = load_err[r]
        rel = e / s["load_da"][r]
        load_rows.append({"region": r, "mean_load_MW": s["load_da"][r].mean(), "mean_error_MW": e.mean(), "std_MW": e.std(),
                          "std_pct_mean_load": 100 * e.std() / s["load_da"][r].mean(),
                          "relative_error_mean": rel.mean(), "relative_error_std": rel.std(),
                          "excess_kurtosis": rts.excess_kurtosis(e), "q01_MW": e.quantile(0.01), "q99_MW": e.quantile(0.99)})
    load_stats = pd.DataFrame(load_rows)
    load_stats.to_csv(out / "load_error_region_statistics.csv", index=False)

    agg_stats = {
        "wind_total_capacity_MW": float(cap[rts.WIND_SITES].sum()),
        "wind_total_error_std_MW": float(agg_wind.std()),
        "wind_total_error_excess_kurtosis": float(rts.excess_kurtosis(agg_wind)),
        "pv_total_capacity_MW": float(cap[s["pv_sites"]].sum()),
        "pv_total_error_std_MW_daylight": float(both.pv_total.std()),
        "pv_total_error_excess_kurtosis_daylight": float(rts.excess_kurtosis(both.pv_total)),
        "corr_wind_total_pv_total_daylight": wind_pv_corr,
        "wind_mean_pairwise_error_corr": offdiag_mean(wc),
        "wind_mean_pairwise_sq_std_error_corr": offdiag_mean(wc2),
        "wind_mean_pairwise_sq_condstd_error_corr": offdiag_mean(wcz2),
        "pv_mean_pairwise_error_corr": offdiag_mean(pc),
        "pv_mean_pairwise_sq_std_error_corr": offdiag_mean(pc2),
        "pv_mean_pairwise_sq_condstd_error_corr": offdiag_mean(pcz2),
        "wind_condstd_excess_kurtosis": wind_cond_kurt,
        "pv_condstd_excess_kurtosis_median": float(np.median(list(pv_cond_kurt.values()))),
        "hours_total": int(len(wind_err)),
        "pv_daylight_hours_per_site_median": int(day.sum().median()),
    }
    rr.write_json(out / "aggregate_statistics.json", agg_stats)

    fig, ax = plt.subplots(figsize=(6, 4))
    for site in rts.WIND_SITES:
        d = spread[(spread.site == site)]
        ax.plot((d.bin_low + d.bin_high) / 2, d.std_pct_cap, marker="o", label=site)
    ax.set_xlabel("day-ahead forecast / capacity")
    ax.set_ylabel("error std (% of capacity)")
    ax.set_title("RTS-GMLC wind: error spread vs forecast level")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "wind_error_spread_vs_level.png", dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4))
    for site in s["pv_sites"]:
        d = spread[(spread.site == site)]
        ax.plot((d.bin_low + d.bin_high) / 2, d.std_pct_cap, color="tab:orange", alpha=0.4)
    ax.set_xlabel("day-ahead forecast / capacity")
    ax.set_ylabel("error std (% of capacity)")
    ax.set_title("RTS-GMLC UPV (25 plants, daylight): error spread vs forecast level")
    fig.tight_layout()
    fig.savefig(out / "pv_error_spread_vs_level.png", dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4))
    z = (agg_wind - agg_wind.mean()) / agg_wind.std()
    bins = np.linspace(-6, 6, 121)
    ax.hist(z, bins=bins, density=True, alpha=0.6, label="aggregate wind error (standardised)")
    ax.plot(bins, stats.norm.pdf(bins), "k-", label="standard normal")
    ax.set_yscale("log")
    ax.set_ylim(1e-5, 1)
    ax.set_xlabel("standardised error")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "wind_total_error_hist_log.png", dpi=130)
    plt.close(fig)

    data_files = [rts.FILES[k] for k in rts.FILES]
    rr.write_manifest(
        out,
        command="python cloud/tests/run_stage_a_rts_errors.py",
        data_files=data_files,
        parameters={
            "error_definition": "real-time 5-min values averaged to hourly minus day-ahead hourly forecast",
            "pv_daylight_rule": "hour kept for a UPV plant when day-ahead forecast > 0 or real-time hourly mean > 0",
            "squared_error_corr": "Pearson correlation of squared errors standardised by site mean and std; condstd variant standardises within 10 forecast-level bins",
            "forecast_level_bins": 10,
        },
        random_settings="not applicable (deterministic statistics)",
    )
    print(agg_stats)
    print(wind_stats.round(3).to_string())
    print(pv_stats.describe().round(3).to_string())
    print(load_stats.round(3).to_string())
    print(out)


if __name__ == "__main__":
    main()
