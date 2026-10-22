"""Stage C: renewable forecast-error model with common fluctuation strength at each operating point."""

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
from tmm_cloud import operating_points as opm
from tmm_cloud import rts_gmlc as rts
from tmm_cloud import run_record as rr

TASK = "phase1_rts118"
K_NEIGHBOURS = 750
N_MODEL_CHECK = 200000
RNG_SEED = 2026102103
LOAD_CV = 0.01


def excess_kurtosis(x):
    z = (x - x.mean()) / x.std()
    return float((z ** 4).mean() - 3)


def check_table(model, E_hist, samples):
    """E_hist: historical neighbour errors clipped to the physical bounds at the operating point."""
    names = model.sites[model.active]
    rows = []
    src_all = {"historical": E_hist} | {k: v[:, model.active] for k, v in samples.items()}
    stats_ = {k: em.dependence_statistics(v) for k, v in src_all.items()}
    for j, nm in enumerate(names):
        row = {"site": nm}
        for k in stats_:
            src = src_all[k]
            row[f"mean_MW_{k}"] = float(src[:, j].mean())
            row[f"std_MW_{k}"] = float(src[:, j].std())
            row[f"excess_kurtosis_{k}"] = float(stats_[k][0][j])
        rows.append(row)
    per_site = pd.DataFrame(rows)
    pairs = []
    m = len(names)
    for i in range(m):
        for j in range(i + 1, m):
            row = {"site_i": names[i], "site_j": names[j]}
            for k in stats_:
                row[f"corr_{k}"] = float(stats_[k][1][i, j])
                row[f"sq_corr_{k}"] = float(stats_[k][2][i, j])
            pairs.append(row)
    pairs = pd.DataFrame(pairs)
    agg = []
    for k in stats_:
        tot = src_all[k].sum(axis=1)
        q = np.quantile(tot, [0.005, 0.01, 0.99, 0.995])
        agg.append({"source": k, "n": len(tot), "total_std_MW": float(tot.std()), "total_excess_kurtosis": excess_kurtosis(tot),
                    "total_q005_MW": q[0], "total_q01_MW": q[1], "total_q99_MW": q[2], "total_q995_MW": q[3]})
    return per_site, pairs, pd.DataFrame(agg)


def pair_type_means(pairs, mapping, cols):
    kind = dict(zip(mapping.site, mapping.kind))
    t = pairs.site_i.map(kind) + "-" + pairs.site_j.map(kind)
    t = t.replace({"pv-wind": "wind-pv"})
    return {c: pairs.groupby(t)[f"sq_corr_{c}"].mean().round(4).to_dict() for c in cols}


def main():
    out = rr.new_result_dir(TASK, "stage_c_error_model")
    t0 = time.time()
    bdir = rr.latest_stage_dir(TASK, "stage_b_network")
    mapping = pd.read_csv(bdir / "renewable_mapping.csv")
    pv_tab = pd.read_csv(bdir / "upv_plant_groups.csv")
    site_table, da, rt, s = opm.site_series(pv_tab)
    if site_table.site.tolist() != mapping.site.tolist():
        raise ValueError("site order differs from stage B mapping")
    sites = mapping.site.values
    rts_cap = mapping.rts_capacity_MW.values
    cap = mapping.capacity_MW.values
    scale = cap / rts_cap
    f_norm = (da.values / rts_cap[None, :])
    E_all = (rt.values - da.values) * scale[None, :]
    rng = np.random.default_rng(RNG_SEED)
    summary = {}
    for name in ("OP1", "OP2"):
        z = np.load(bdir / f"operating_point_{name}.npz")
        base = 100.0
        f_op = z["Pren"] * base
        active = f_op > 0
        group = (mapping.kind.values[active] == "pv").astype(int)
        model, info = em.fit(sites, f_norm, E_all, f_op, cap, active, K_NEIGHBOURS, group=group)
        np.savez(out / f"error_model_{name}.npz", **em.to_npz_dict(model))
        hours = da.index[model.neighbours]
        pd.DataFrame({"hour_start": hours, "distance": model.distance}).to_csv(out / f"neighbour_hours_{name}.csv", index=False)
        samples = {}
        bounds = {}
        kinds = ("t", "grouped_t", "gauss") if len(np.unique(group)) > 1 else ("t", "gauss")
        for kind in kinds:
            e, w = em.sample_errors(model, N_MODEL_CHECK, rng, kind)
            samples[kind] = e
            lo, hi = model.lower(), model.upper()
            ea = e[:, active]
            bounds[kind] = {"fraction_at_lower_bound": dict(zip(sites[active], (ea <= lo + 1e-9).mean(axis=0).round(5))),
                            "fraction_at_upper_bound": dict(zip(sites[active], (ea >= hi - 1e-9).mean(axis=0).round(5)))}
        E_hist_op = np.clip(model.errors, model.lower(), model.upper())
        per_site, pairs, agg = check_table(model, E_hist_op, samples)
        raw_tot = model.errors.sum(axis=1)
        per_site.to_csv(out / f"model_check_sites_{name}.csv", index=False)
        pairs.to_csv(out / f"model_check_pairs_{name}.csv", index=False)
        agg.to_csv(out / f"model_check_total_{name}.csv", index=False)
        pd.DataFrame(model.R_t, index=sites[active], columns=sites[active]).to_csv(out / f"t_copula_corr_{name}.csv")
        pd.DataFrame(model.R_g, index=sites[active], columns=sites[active]).to_csv(out / f"gauss_copula_corr_{name}.csv")
        hist_hours_in_window = int(np.isin(hours.hour, opm.NIGHT_HOURS).sum()) if name == "OP1" else int(np.isin(hours.hour, opm.MIDDAY_HOURS).sum())
        summary[name] = {
            "active_sites": sites[active].tolist(),
            "k": K_NEIGHBOURS,
            "distance_rule": "Euclidean distance between day-ahead forecast / capacity vectors of all 10 mapped sites",
            "distance_median": float(np.median(model.distance)), "distance_max": float(model.distance.max()),
            "neighbour_hour_of_day_counts": pd.Series(hours.hour).value_counts().sort_index().to_dict(),
            "neighbour_hours_in_selection_window": hist_hours_in_window,
            "neighbour_months": pd.Series(hours.month).value_counts().sort_index().to_dict(),
            "kde_bandwidth_MW": dict(zip(sites[active], model.bandwidth.round(4))),
            "t_copula_nu": model.nu, "t_copula_nu_ci95_profile_likelihood": model.nu_ci95,
            "heavy_tailed_model": model.heavy_kind, "heavy_tailed_nu_by_group": model.nu_groups.tolist(),
            "groups": {"0": "wind", "1": "photovoltaic"},
            **info,
            "raw_neighbour_total_error_std_MW": float(raw_tot.std()),
            "raw_neighbour_total_error_excess_kurtosis": excess_kurtosis(raw_tot),
            "historical_reference_rule": "neighbour errors clipped to [-forecast, capacity - forecast] of the operating point",
            "bounds": bounds,
            "mean_pairwise_corr": {c: float(pairs[f"corr_{c}"].mean()) for c in ("historical",) + kinds},
            "mean_pairwise_sq_corr": {c: float(pairs[f"sq_corr_{c}"].mean()) for c in ("historical",) + kinds},
            "mean_pairwise_sq_corr_by_pair_type": pair_type_means(pairs, mapping, ("historical",) + kinds),
            "site_excess_kurtosis": {c: per_site[f"excess_kurtosis_{c}"].round(3).tolist() for c in ("historical",) + kinds},
            "total_error": agg.set_index("source").round(4).to_dict(orient="index"),
        }

        fig, ax = plt.subplots(figsize=(5, 3.5))
        ax.semilogx(model.nu_profile[:, 0], model.nu_profile[:, 1] - model.nu_profile[:, 1].max(), "-")
        ax.axvline(model.nu, color="k", ls="--", lw=0.8)
        ax.set_ylim(-30, 1)
        ax.set_xlabel("t copula degrees of freedom")
        ax.set_ylabel("log-likelihood - max")
        ax.set_title(f"{name}: profile likelihood, nu = {model.nu:.2f}")
        fig.tight_layout()
        fig.savefig(out / f"nu_profile_{name}.png", dpi=130)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(6, 4))
        tot_h = E_hist_op.sum(axis=1)
        hk = "grouped_t" if model.heavy_kind == "grouped_t" else "t"
        bins = np.linspace(min(tot_h.min(), samples[hk][:, active].sum(1).min()), max(tot_h.max(), samples[hk][:, active].sum(1).max()), 80)
        ax.hist(tot_h, bins=bins, density=True, histtype="step", lw=1.5, label=f"historical neighbours at OP bounds (k={K_NEIGHBOURS})")
        ax.hist(samples[hk][:, active].sum(1), bins=bins, density=True, histtype="step", lw=1.2, label=f"{hk} copula model")
        ax.hist(samples["gauss"][:, active].sum(1), bins=bins, density=True, histtype="step", lw=1.2, label="Gaussian copula model")
        ax.set_yscale("log")
        ax.set_xlabel("total renewable forecast error (MW)")
        ax.set_ylabel("density")
        ax.set_title(f"{name}: total renewable error")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(out / f"total_error_density_{name}.png", dpi=130)
        plt.close(fig)

    summary["load_fluctuation"] = {"rule": "independent Gaussian per bus, P and Q multiplied by the same factor 1 + 0.01 xi (constant power factor)", "cv": LOAD_CV}
    summary["elapsed_s"] = time.time() - t0
    rr.write_json(out / "stage_c_summary.json", summary)
    rr.write_manifest(
        out,
        command="python cloud/tests/run_stage_c_error_model.py",
        data_files=[rts.FILES[k] for k in rts.FILES] + [mc.PGLIB118, bdir / "operating_point_OP1.npz", bdir / "operating_point_OP2.npz",
                                                         bdir / "renewable_mapping.csv", bdir / "upv_plant_groups.csv"],
        parameters={"k_neighbours": K_NEIGHBOURS, "marginal": "Gaussian kernel, Silverman bandwidth, values outside [0, capacity] set to the bound",
                    "t_copula": "R from Kendall tau (sin(pi tau / 2), nearest correlation matrix), nu by bounded maximum likelihood on log nu in [2.05, 500]",
                    "gauss_copula": "R = correlation of normal scores of pseudo-observations",
                    "model_check_samples": N_MODEL_CHECK, "load_cv": LOAD_CV, "stage_b_result": rr.relpath(bdir)},
        random_settings={"generator": "numpy default_rng (PCG64)", "seed_value": RNG_SEED, "use": "model check samples only"},
    )
    import json
    print(json.dumps({k: {kk: v[kk] for kk in ("t_copula_nu", "t_copula_nu_ci95_profile_likelihood", "lr_statistic_t_vs_gauss", "mean_pairwise_corr", "mean_pairwise_sq_corr", "site_excess_kurtosis", "total_error", "distance_max", "neighbour_hour_of_day_counts")} for k, v in summary.items() if k.startswith("OP")}, indent=1, default=str))
    print(rr.relpath(out))


if __name__ == "__main__":
    main()
