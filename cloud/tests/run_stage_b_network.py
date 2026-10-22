"""Stage B: PGLib IEEE118 network, three regions, RTS-GMLC renewable placement, two operating points,
AC OPF base dispatch, distributed active-power balancing, and cross-check of the batched Newton-Raphson
solver against pandapower on random scenarios."""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from tmm_cloud import matpower_case as mc
from tmm_cloud import network as nw
from tmm_cloud import operating_points as opm
from tmm_cloud import pandapower_bridge as ppb
from tmm_cloud import powerflow as pf
from tmm_cloud import regions as rg
from tmm_cloud import rts_gmlc as rts
from tmm_cloud import run_record as rr

TASK = "phase1_rts118"
RENEWABLE_SHARE = 0.25
RENEWABLE_SHARE_SECOND = 0.20
PV_GROUPS = {1: 2, 2: 1, 3: 3}
N_CROSSCHECK = 200
LOAD_CV = 0.01
RNG_SEED = 2026102101
OPF_HEADROOM = 0.05
V_TOL = 1e-4
LOADING_TOL = 1e-3


def place_sites(net, regions, site_table, scaled_cap):
    region_of_area = {1: 1, 2: 2, 3: 3}
    cand = rg.host_candidates(net, regions, net.Pd * net.baseMVA)
    used = set()
    host = {}
    order = site_table.assign(cap=site_table.site.map(scaled_cap)).sort_values(["kind", "cap"], ascending=[False, False])
    for _, row in order.iterrows():
        reg = region_of_area[row.rts_area]
        c = cand[(cand.region == reg) & (~cand.bus.isin(used)) & (cand.connected_rateA_MVA >= 2 * scaled_cap[row.site])]
        if c.empty:
            raise RuntimeError(f"no host bus for {row.site}")
        b = int(c.iloc[0].bus)
        used.add(b)
        host[row.site] = b
    return host, cand


def imbalance_range(net, op):
    """Largest negative and positive imbalance for which every sharing unit stays within [Pmin, Pmax]."""
    a = op.alpha > 0
    up = ((net.Pmax - op.Pg0)[a] / op.alpha[a]).min() * net.baseMVA
    dn = ((net.Pmin - op.Pg0)[a] / op.alpha[a]).max() * net.baseMVA
    return [float(dn), float(up)]


def base_case_report(net, res, op, regions, name, out):
    V = res["V"][0]
    Vm = np.abs(V)
    bus = pd.DataFrame({"bus": net.bus_id, "region": regions, "Vm_pu": Vm, "Va_deg": np.rad2deg(np.angle(V)),
                        "Vmin": net.Vmin, "Vmax": net.Vmax, "Pd_MW": op.Pd * net.baseMVA, "Qd_Mvar": op.Qd * net.baseMVA})
    bus.to_csv(out / f"base_bus_{name}.csv", index=False)
    Sf, St = np.abs(res["Sf_MVA"][0]), np.abs(res["St_MVA"][0])
    br = pd.DataFrame({"from_bus": net.bus_id[net.f], "to_bus": net.bus_id[net.t], "rateA_MVA": net.rateA,
                       "S_from_MVA": Sf, "S_to_MVA": St, "loading": np.maximum(Sf, St) / net.rateA})
    br.to_csv(out / f"base_branch_{name}.csv", index=False)
    gen = pd.DataFrame({"bus": net.bus_id[net.gen_bus], "Pmax_MW": net.Pmax * net.baseMVA, "Pmin_MW": net.Pmin * net.baseMVA,
                        "Qmax_Mvar": net.Qmax * net.baseMVA, "Qmin_Mvar": net.Qmin * net.baseMVA,
                        "P0_MW": op.Pg0 * net.baseMVA, "Vset_pu": op.Vset, "participation": op.alpha,
                        "P_base_MW": res["Pg_MW"][0], "Q_base_Mvar": res["Qg_Mvar"][0]})
    mode = np.zeros(net.n_gen, int)
    reg_gen = np.flatnonzero(net.gen_bus != net.ref)
    mode[reg_gen] = res["mode"][0]
    gen["mode_base"] = mode
    gen["is_reference"] = net.gen_bus == net.ref
    gen.to_csv(out / f"base_dispatch_{name}.csv", index=False)
    return {
        "Vm_min_pu": float(Vm.min()), "Vm_min_bus": int(net.bus_id[Vm.argmin()]),
        "Vm_max_pu": float(Vm.max()), "Vm_max_bus": int(net.bus_id[Vm.argmax()]),
        "buses_below_Vmin": int((Vm < net.Vmin - V_TOL).sum()), "buses_above_Vmax": int((Vm > net.Vmax + V_TOL).sum()),
        "min_margin_to_Vmin_pu": float((Vm - net.Vmin).min()), "min_margin_to_Vmax_pu": float((net.Vmax - Vm).min()),
        "max_branch_loading": float(br.loading.max()),
        "max_loading_branch": f"{int(br.from_bus[br.loading.idxmax()])}-{int(br.to_bus[br.loading.idxmax()])}",
        "branches_over_rateA": int((br.loading > 1 + LOADING_TOL).sum()),
        "branches_above_90pct": int((br.loading > 0.9).sum()),
        "units_at_Q_limit_base": int((mode != 0).sum()),
        "units_within_1Mvar_of_Q_limit_base": [int(net.bus_id[net.gen_bus[g]]) for g in reg_gen
                                              if min(abs(res["Qg_Mvar"][0][g] - net.Qmax[g] * net.baseMVA),
                                                     abs(res["Qg_Mvar"][0][g] - net.Qmin[g] * net.baseMVA)) < 1.0],
        "units_P_outside_limits_base": int(((res["Pg_MW"][0] > net.Pmax * net.baseMVA + 1e-6) | (res["Pg_MW"][0] < net.Pmin * net.baseMVA - 1e-6)).sum()),
        "imbalance_base_MW": float(res["lam_MW"][0]),
        "total_load_MW": float(op.Pd.sum() * net.baseMVA),
        "total_renewable_forecast_MW": float(op.Pren.sum() * net.baseMVA),
        "total_conventional_MW": float(op.Pg0.sum() * net.baseMVA),
        "losses_MW": float(op.Pg0.sum() * net.baseMVA + op.Pren.sum() * net.baseMVA - op.Pd.sum() * net.baseMVA + res["lam_MW"][0]),
        "base_pf_fail_code": int(res["fail"][0]),
    }


def random_scenarios(rng, n, op, site_kind, errors_hist, is_night, net):
    """Scenarios for the solver cross-check: historical hourly errors of the mapped sites drawn at random
    (same hour for all sites), clipped to [0, capacity]; independent Gaussian load fluctuation."""
    E = errors_hist.values
    rows = rng.integers(0, len(E), n)
    Pren = op.Pren[None, :] + E[rows] / net.baseMVA
    Pren = np.clip(Pren, 0, op.ren_cap[None, :])
    if is_night:
        Pren[:, site_kind == "pv"] = 0.0
    k = 1 + LOAD_CV * rng.standard_normal((n, net.n))
    return op.Pd[None, :] * k, op.Qd[None, :] * k, Pren


def build_operating_point(net, pm, regions, host, site_table, scaled_cap, da_scaled, hour, rts_load_da):
    fac = opm.regional_load_factor(rts_load_da, hour)
    fvec = np.array([fac[str(r)] for r in regions])
    Pd = net.Pd * fvec
    Qd = net.Qd * fvec
    sites = site_table.site.tolist()
    ren_bus = net.idx([host[s] for s in sites])
    Pren_MW = da_scaled.loc[hour, sites].values
    cap_MW = np.array([scaled_cap[s] for s in sites])
    pm.set_injections(Pd * net.baseMVA, Qd * net.baseMVA, Pren_MW)
    opf0 = pm.run_opf(headroom=0.0)
    opf = pm.run_opf(headroom=OPF_HEADROOM)
    Pmax_MW = net.Pmax * net.baseMVA
    act = Pmax_MW > 0
    opf["first_dispatch_without_headroom"] = {
        "Pg_MW": opf0["Pg_MW"], "cost": opf0["cost"],
        "units_at_Pmin": int((act & (opf0["Pg_MW"] < net.Pmin * net.baseMVA + 0.01)).sum()),
        "units_at_Pmax": int((act & (opf0["Pg_MW"] > Pmax_MW - 0.01)).sum()),
    }
    alpha = net.Pmax / net.Pmax.sum()
    op = pf.OperatingPoint(Pg0=opf["Pg_MW"] / net.baseMVA, Vset=opf["Vset"], alpha=alpha, Pd=Pd, Qd=Qd,
                           ren_bus=ren_bus, Pren=Pren_MW / net.baseMVA, ren_cap=cap_MW / net.baseMVA)
    return op, opf, fac


def run_share(share, net, case, regions, site_table, da, rt, s, hours, out, tag):
    total_load = net.Pd.sum() * net.baseMVA
    scale = share * total_load / site_table.rts_capacity_MW.sum()
    scaled_cap = {r.site: r.rts_capacity_MW * scale for r in site_table.itertuples()}
    host, cand = place_sites(net, regions, site_table, scaled_cap)
    pm = ppb.PandapowerModel(case, net)
    sites = site_table.site.tolist()
    pm.add_renewables(net.idx([host[x] for x in sites]))
    da_scaled = da * scale
    ops, reports = {}, {}
    for name, hour in hours.items():
        op, opf, fac = build_operating_point(net, pm, regions, host, site_table, scaled_cap, da_scaled, hour, s["load_da"])
        solver = pf.BatchSolver(net, op)
        res = solver.solve(op.Pd[None, :], op.Qd[None, :], op.Pren[None, :])
        rep = base_case_report(net, res, op, regions, f"{name}{tag}", out)
        rep.update({"hour_start": str(hour), "regional_load_factor": fac, "opf_cost": opf["cost"],
                    "opf_headroom_fraction_of_Pmax": OPF_HEADROOM,
                    "first_dispatch_without_headroom": {k: v for k, v in opf["first_dispatch_without_headroom"].items() if k != "Pg_MW"},
                    "upward_reserve_MW": float(((net.Pmax - op.Pg0) * net.baseMVA)[net.Pmax > 0].sum()),
                    "downward_reserve_MW": float(((op.Pg0 - net.Pmin) * net.baseMVA)[net.Pmax > 0].sum()),
                    "imbalance_range_without_any_unit_leaving_limits_MW": imbalance_range(net, op)})
        pd.DataFrame({"bus": net.bus_id[net.gen_bus], "P_opf_headroom_MW": op.Pg0 * net.baseMVA,
                      "P_opf_without_headroom_MW": opf["first_dispatch_without_headroom"]["Pg_MW"]}
                     ).to_csv(out / f"opf_dispatch_comparison_{name}{tag}.csv", index=False)
        ops[name] = (op, opf)
        reports[name] = rep
    feasible = all(r["base_pf_fail_code"] == 0 and r["branches_over_rateA"] == 0 and r["buses_below_Vmin"] == 0
                   and r["buses_above_Vmax"] == 0 for r in reports.values())
    return {"scale": scale, "scaled_cap": scaled_cap, "host": host, "cand": cand, "pm": pm, "ops": ops,
            "reports": reports, "feasible": feasible, "da_scaled": da_scaled}


def plot_base_cases(out, names):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for name in names:
        bus = pd.read_csv(out / f"base_bus_{name}.csv")
        br = pd.read_csv(out / f"base_branch_{name}.csv")
        fig, ax = plt.subplots(figsize=(8, 3.5))
        for r, c in zip((1, 2, 3), ("tab:blue", "tab:orange", "tab:green")):
            d = bus[bus.region == r]
            ax.plot(d.bus, d.Vm_pu, "o", ms=3, color=c, label=f"region {r}")
        ax.axhline(bus.Vmax.iloc[0], color="k", ls="--", lw=0.8)
        ax.axhline(bus.Vmin.iloc[0], color="k", ls="--", lw=0.8)
        ax.set_xlabel("bus")
        ax.set_ylabel("|V| (p.u.)")
        ax.set_title(f"{name} base case bus voltage magnitude")
        ax.legend(ncol=3, fontsize=8)
        fig.tight_layout()
        fig.savefig(out / f"base_voltage_{name}.png", dpi=130)
        plt.close(fig)
        fig, ax = plt.subplots(figsize=(8, 3.5))
        ax.plot(np.sort(br.loading.values)[::-1] * 100, ".-", ms=3)
        ax.axhline(100, color="k", ls="--", lw=0.8)
        ax.set_xlabel("branch rank")
        ax.set_ylabel("max(|S_from|, |S_to|) / rateA (%)")
        ax.set_title(f"{name} base case branch loading")
        fig.tight_layout()
        fig.savefig(out / f"base_branch_loading_{name}.png", dpi=130)
        plt.close(fig)


def main():
    out = rr.new_result_dir(TASK, "stage_b_network")
    t0 = time.time()
    case = mc.load_case()
    net = nw.from_matpower(case)

    pm0 = ppb.PandapowerModel(case, net)
    check = {
        "file": {"bus": int(case["bus"].shape[0]), "gen": int(case["gen"].shape[0]), "branch": int(case["branch"].shape[0]),
                 "Pd_MW": float(case["bus"][:, 2].sum()), "Qd_Mvar": float(case["bus"][:, 3].sum()),
                 "transformer_branches": int((case["branch"][:, 8] != 0).sum()), "reference_bus": int(net.bus_id[net.ref]),
                 "units_Pmax_positive": int((net.Pmax > 0).sum()), "total_Pmax_MW": float(net.Pmax.sum() * net.baseMVA)},
        "parsed_network": {"bus": net.n, "gen": net.n_gen, "branch": len(net.f), "Pd_MW": float(net.Pd.sum() * net.baseMVA),
                           "Qd_Mvar": float(net.Qd.sum() * net.baseMVA)},
        "pandapower": pm0.counts,
        "Ybus_max_abs_difference_pu": float(np.abs(pm0.pnet._ppc["internal"]["Ybus"].toarray() - net.Y).max())
        if pm0.pnet.get("_ppc") is not None else None,
    }
    import pandapower as pp
    pp.runpp(pm0.pnet)
    check["Ybus_max_abs_difference_pu"] = float(np.abs(pm0.pnet._ppc["internal"]["Ybus"].toarray() - net.Y).max())
    rr.write_json(out / "network_check.json", check)

    regions, eig = rg.spectral_regions(net)
    conn = rg.region_connectivity(net, regions)
    ties = rg.tie_lines(net, regions)
    pd.DataFrame({"bus": net.bus_id, "region": regions}).to_csv(out / "regions.csv", index=False)
    pd.DataFrame({"from_bus": net.bus_id[net.f[ties]], "to_bus": net.bus_id[net.t[ties]],
                  "from_region": regions[net.f[ties]], "to_region": regions[net.t[ties]], "rateA_MVA": net.rateA[ties]}
                 ).to_csv(out / "tie_lines.csv", index=False)
    region_summary = {int(r): {"buses": int((regions == r).sum()), "load_MW": float(net.Pd[regions == r].sum() * net.baseMVA),
                               "units": int(np.isin(net.gen_bus, np.flatnonzero(regions == r)).sum()),
                               "Pmax_MW": float(net.Pmax[np.isin(net.gen_bus, np.flatnonzero(regions == r))].sum() * net.baseMVA),
                               "connected_components": conn[int(r)],
                               "bus_list": net.bus_id[regions == r].tolist()} for r in (1, 2, 3)}

    pv_tab = rg.pv_groups(PV_GROUPS)
    pv_tab.to_csv(out / "upv_plant_groups.csv", index=False)
    site_table, da, rt, s = opm.site_series(pv_tab)
    kinds = dict(zip(site_table.site, site_table.kind))
    hours, sel_info = opm.choose_operating_points(da, kinds, s["load_da"])

    run = run_share(RENEWABLE_SHARE, net, case, regions, site_table, da, rt, s, hours, out, "")
    share_used = RENEWABLE_SHARE
    second = None
    if not run["feasible"]:
        second = run_share(RENEWABLE_SHARE_SECOND, net, case, regions, site_table, da, rt, s, hours, out, "_share20")
        if second["feasible"]:
            run, share_used = second, RENEWABLE_SHARE_SECOND

    mapping = site_table.copy()
    mapping["ieee118_region"] = mapping.rts_area
    mapping["host_bus"] = mapping.site.map(run["host"])
    mapping["capacity_MW"] = mapping.site.map(run["scaled_cap"])
    mapping["host_connected_rateA_MVA"] = mapping.host_bus.map(run["cand"].set_index("bus").connected_rateA_MVA)
    mapping["host_base_load_MW"] = mapping.host_bus.map(run["cand"].set_index("bus").load_MW)
    mapping["upv_plants"] = mapping.site.map(pv_tab.groupby("group")["GEN UID"].apply(lambda x: " ".join(x)))
    mapping.to_csv(out / "renewable_mapping.csv", index=False)
    run["cand"].to_csv(out / "host_bus_candidates.csv", index=False)

    op_info = {}
    for name, hour in hours.items():
        op, opf = run["ops"][name]
        np.savez(out / f"operating_point_{name}.npz", Pg0=op.Pg0, Vset=op.Vset, alpha=op.alpha, Pd=op.Pd, Qd=op.Qd,
                 ren_bus=op.ren_bus, Pren=op.Pren, ren_cap=op.ren_cap, sites=np.array(site_table.site),
                 site_kind=np.array(site_table.site.map(kinds)), regions=regions, hour=str(hour))
        op_info[name] = {"hour_start": str(hour),
                         "forecast_MW": dict(zip(site_table.site, (op.Pren * net.baseMVA).round(3))),
                         "forecast_wind_MW": float(op.Pren[site_table.kind.values == "wind"].sum() * net.baseMVA),
                         "forecast_pv_MW": float(op.Pren[site_table.kind.values == "pv"].sum() * net.baseMVA),
                         **run["reports"][name]}

    rng = np.random.default_rng(RNG_SEED)
    cross = {}
    err_hist = (rt - da) * run["scale"]
    for name, hour in hours.items():
        op, opf = run["ops"][name]
        night = hour.hour in opm.NIGHT_HOURS
        hist = err_hist.copy()
        if not night:
            hist = hist[np.isin(hist.index.hour, opm.MIDDAY_HOURS)]
        Pd, Qd, Pren = random_scenarios(rng, N_CROSSCHECK, op, site_table.kind.values, hist, night, net)
        solver = pf.BatchSolver(net, op)
        ts = time.time()
        res = solver.solve(Pd, Qd, Pren)
        t_own = time.time() - ts
        pm = run["pm"]
        pm.set_dispatch(op.Pg0 * net.baseMVA, op.Vset, op.alpha)
        rows = []
        ts = time.time()
        for k in range(N_CROSSCHECK):
            pm.set_injections(Pd[k] * net.baseMVA, Qd[k] * net.baseMVA, Pren[k] * net.baseMVA)
            try:
                r = pm.run_pf()
                ok = r["converged"]
            except Exception as exc:
                ok = False
                r = None
                print("pandapower failure", name, k, exc)
            own_mode = np.zeros(net.n_gen, int)
            own_mode[solver.reg_gen] = res["mode"][k]
            row = {"scenario": k, "own_fail_code": int(res["fail"][k]), "pandapower_converged": ok,
                   "own_iterations": int(res["iterations"][k]), "own_rounds": int(res["rounds"][k]),
                   "imbalance_MW": float(res["lam_MW"][k]), "units_at_limit_own": int((own_mode != 0).sum())}
            if ok:
                ppm = r["mode"]
                row.update({
                    "max_abs_V_diff_pu": float(np.abs(r["V"] - res["V"][k]).max()),
                    "max_abs_Pg_diff_MW": float(np.abs(r["Pg_MW"] - res["Pg_MW"][k]).max()),
                    "max_abs_Qg_diff_Mvar": float(np.abs(r["Qg_Mvar"] - res["Qg_Mvar"][k]).max()),
                    "units_at_limit_pandapower": int((ppm != 0).sum()),
                    "modes_identical": bool(np.array_equal(ppm, own_mode)),
                })
            rows.append(row)
        t_pp = time.time() - ts
        df = pd.DataFrame(rows)
        df.to_csv(out / f"crosscheck_{name}.csv", index=False)
        both = df[(df.own_fail_code == 0) & df.pandapower_converged]
        cross[name] = {
            "scenarios": N_CROSSCHECK, "own_converged": int((df.own_fail_code == 0).sum()),
            "pandapower_converged": int(df.pandapower_converged.sum()),
            "compared": int(len(both)),
            "max_abs_V_diff_pu": float(both.max_abs_V_diff_pu.max()),
            "max_abs_Pg_diff_MW": float(both.max_abs_Pg_diff_MW.max()),
            "max_abs_Qg_diff_Mvar": float(both.max_abs_Qg_diff_Mvar.max()),
            "modes_identical": int(both.modes_identical.sum()),
            "scenarios_with_units_at_limit": int((both.units_at_limit_own > 0).sum()),
            "max_units_at_limit": int(both.units_at_limit_own.max()),
            "imbalance_MW_range": [float(df.imbalance_MW.min()), float(df.imbalance_MW.max())],
            "own_time_s": t_own, "pandapower_time_s": t_pp,
            "own_iterations_mean": float(df.own_iterations.mean()),
        }

    summary = {
        "network_check": check,
        "region_rule": "normalised spectral clustering (Ng-Jordan-Weiss, k=3) on branch |series admittance| graph, k-means 50 restarts",
        "laplacian_smallest_eigenvalues": eig,
        "regions": region_summary,
        "tie_lines": int(len(ties)),
        "rts_area_to_region": {"1": 1, "2": 2, "3": 3},
        "renewable_share_of_base_load": share_used,
        "share_25pct_feasible": run["feasible"] if share_used == RENEWABLE_SHARE else False,
        "share_20pct_tried": second is not None,
        "capacity_scale_factor": run["scale"],
        "total_renewable_capacity_MW": float(sum(run["scaled_cap"].values())),
        "operating_point_selection": sel_info,
        "operating_points": op_info,
        "participation_rule": "alpha_g = Pmax_g / sum(Pmax), units with Pmax = 0 (synchronous condensers) do not share",
        "crosscheck": cross,
        "elapsed_s": time.time() - t0,
    }
    rr.write_json(out / "stage_b_summary.json", summary)
    plot_base_cases(out, list(hours))
    rr.write_manifest(
        out,
        command="python cloud/tests/run_stage_b_network.py",
        data_files=[mc.PGLIB118] + [rts.FILES[k] for k in rts.FILES],
        parameters={
            "renewable_share": RENEWABLE_SHARE, "second_share": RENEWABLE_SHARE_SECOND, "pv_groups_per_rts_area": PV_GROUPS,
            "host_rule": "PQ bus with positive base load, degree >= 2, connected rateA sum >= 2 x site capacity; sites placed largest first (wind before UPV) on the remaining candidate with largest connected rateA in the mapped region",
            "load_rule": "bus load x (RTS regional day-ahead load at the hour / RTS regional annual mean), constant power factor",
            "opf_headroom_fraction_of_Pmax": OPF_HEADROOM,
            "dispatch": "pandapower runopp (PIPS AC OPF, PGLib quadratic costs, rateA, bus voltage limits, unit P/Q limits); OPF unit P and bus voltage at unit buses used as base schedule and setpoints",
            "balancing": "all units with Pmax > 0 share the imbalance by alpha = Pmax / sum Pmax, reference bus fixes angle and magnitude",
            "q_limits": "all violating units switch together, no return to regulation, eps_q = 5e-8 p.u. (5e-6 Mvar), reference unit unlimited",
            "newton_tolerance_pu": 1e-10, "pandapower_tolerance_mva": 1e-8,
            "crosscheck_scenarios": N_CROSSCHECK, "load_cv": LOAD_CV,
        },
        random_settings={"generator": "numpy default_rng (PCG64)", "seed_value": RNG_SEED,
                         "use": "cross-check scenarios: historical hourly errors of a random hour (midday hours for OP2), load factor 1 + 0.01 N(0,1) per bus"},
    )
    print(rr.relpath(out))
    import json
    print(json.dumps({k: summary[k] for k in ["renewable_share_of_base_load", "capacity_scale_factor", "crosscheck"]}, indent=1, default=str))
    for n_, v in op_info.items():
        print(n_, {k: v[k] for k in v if k != "forecast_MW"})


if __name__ == "__main__":
    main()
