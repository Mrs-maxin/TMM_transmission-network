"""Mapped renewable sites on IEEE118, their RTS-GMLC series, and the two representative operating hours."""

import numpy as np
import pandas as pd

from . import rts_gmlc as rts

NIGHT_HOURS = [22, 23, 0, 1, 2, 3, 4]      # hour_start values counted as night
MIDDAY_HOURS = [11, 12, 13]                 # hour_start values counted as midday


def site_series(pv_group_table):
    """Day-ahead and real-time hourly series (MW, RTS capacity) for the 4 wind farms and the UPV groups."""
    s = rts.load_series()
    cap = s["capacity"]
    sites, da, rt, kind, rts_area = [], {}, {}, {}, {}
    for w in rts.WIND_SITES:
        sites.append(w)
        da[w] = s["wind_da"][w]
        rt[w] = s["wind_rt"][w]
        kind[w] = "wind"
        rts_area[w] = int(w[0])
    groups = pv_group_table.groupby("group")
    for gname, g in groups:
        plants = g["GEN UID"].tolist()
        sites.append(gname)
        da[gname] = s["pv_da"][plants].sum(axis=1)
        rt[gname] = s["pv_rt"][plants].sum(axis=1)
        kind[gname] = "pv"
        rts_area[gname] = int(g["area"].iloc[0])
    capacity = {w: float(cap[w]) for w in rts.WIND_SITES}
    capacity.update({gname: float(g["PMax MW"].sum()) for gname, g in groups})
    table = pd.DataFrame({"site": sites, "kind": [kind[x] for x in sites], "rts_area": [rts_area[x] for x in sites],
                          "rts_capacity_MW": [capacity[x] for x in sites]})
    return table, pd.DataFrame(da)[sites], pd.DataFrame(rt)[sites], s


def choose_operating_points(da_norm, kinds, load_da):
    """OP1: night hour with aggregate wind forecast in the annual 75th-90th percentile band, closest to the
    band centre (82.5th percentile). OP2: midday hour with aggregate UPV forecast at or above the 90th
    percentile of midday hours and aggregate wind forecast in the annual 40th-60th percentile band; the
    hour with wind closest to the annual median is taken."""
    wind_cols = [c for c in da_norm.columns if kinds[c] == "wind"]
    pv_cols = [c for c in da_norm.columns if kinds[c] == "pv"]
    wind = da_norm[wind_cols].sum(axis=1)
    pv = da_norm[pv_cols].sum(axis=1)
    hour = da_norm.index.hour
    q75, q825, q90 = np.quantile(wind, [0.75, 0.825, 0.90])
    night = np.isin(hour, NIGHT_HOURS)
    c1 = night & (wind >= q75) & (wind <= q90)
    i1 = (wind[c1] - q825).abs().idxmin()
    mid = np.isin(hour, MIDDAY_HOURS)
    pv90 = np.quantile(pv[mid], 0.90)
    q40, q50, q60 = np.quantile(wind, [0.40, 0.50, 0.60])
    c2 = mid & (pv >= pv90) & (wind >= q40) & (wind <= q60)
    i2 = (wind[c2] - q50).abs().idxmin()
    info = {
        "wind_aggregate_quantiles_MWnorm": {"q40": q40, "q50": q50, "q60": q60, "q75": q75, "q825": q825, "q90": q90},
        "pv_midday_q90_MWnorm": pv90,
        "op1_candidates": int(c1.sum()),
        "op2_candidates": int(c2.sum()),
    }
    return {"OP1": i1, "OP2": i2}, info


def regional_load_factor(load_da, hour):
    return (load_da.loc[hour] / load_da.mean()).to_dict()
