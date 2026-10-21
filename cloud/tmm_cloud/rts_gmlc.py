"""RTS-GMLC 2020 renewable and regional load time series and their hourly forecast errors.

Day-ahead files hold 8784 hourly periods (Period 1..24 per day). Real-time files hold
105408 five-minute periods (Period 1..288 per day). Real-time hour h is the mean of
periods 12(h-1)+1 .. 12h of the same day.
"""

import numpy as np
import pandas as pd

from . import DATA_ROOT

RTS_DIR = DATA_ROOT / "rts_gmlc"
WIND_SITES = ["309_WIND_1", "317_WIND_1", "303_WIND_1", "122_WIND_1"]
FILES = {
    "wind_da": RTS_DIR / "DAY_AHEAD_wind.csv",
    "wind_rt": RTS_DIR / "REAL_TIME_wind.csv",
    "pv_da": RTS_DIR / "DAY_AHEAD_pv.csv",
    "pv_rt": RTS_DIR / "REAL_TIME_pv.csv",
    "load_da": RTS_DIR / "DAY_AHEAD_regional_Load.csv",
    "load_rt": RTS_DIR / "REAL_TIME_regional_Load.csv",
    "gen": RTS_DIR / "gen.csv",
    "bus": RTS_DIR / "bus.csv",
}


def _hourly_index(df):
    return pd.to_datetime(dict(year=df.Year, month=df.Month, day=df.Day)) + pd.to_timedelta(df.Period - 1, unit="h")


def read_day_ahead(path):
    df = pd.read_csv(path)
    idx = _hourly_index(df)
    out = df.drop(columns=["Year", "Month", "Day", "Period"])
    out.index = pd.DatetimeIndex(idx, name="hour_start")
    return out


def read_real_time_hourly(path):
    df = pd.read_csv(path)
    hour = (df.Period - 1) // 12 + 1
    key = pd.to_datetime(dict(year=df.Year, month=df.Month, day=df.Day)) + pd.to_timedelta(hour - 1, unit="h")
    out = df.drop(columns=["Year", "Month", "Day", "Period"]).groupby(key.values).mean()
    out.index = pd.DatetimeIndex(out.index, name="hour_start")
    return out


def generator_table():
    return pd.read_csv(FILES["gen"])


def bus_table():
    return pd.read_csv(FILES["bus"])


def load_series():
    """Return dict with day-ahead and hourly real-time series for wind, UPV and regional load."""
    gen = generator_table()
    pv_sites = gen.loc[gen["Unit Type"] == "PV", "GEN UID"].tolist()
    cap = gen.set_index("GEN UID")["PMax MW"]
    s = {
        "wind_da": read_day_ahead(FILES["wind_da"])[WIND_SITES],
        "wind_rt": read_real_time_hourly(FILES["wind_rt"])[WIND_SITES],
        "pv_da": read_day_ahead(FILES["pv_da"])[pv_sites],
        "pv_rt": read_real_time_hourly(FILES["pv_rt"])[pv_sites],
        "load_da": read_day_ahead(FILES["load_da"]),
        "load_rt": read_real_time_hourly(FILES["load_rt"]),
        "capacity": cap,
        "pv_sites": pv_sites,
        "gen": gen,
    }
    for a, b in [("wind_da", "wind_rt"), ("pv_da", "pv_rt"), ("load_da", "load_rt")]:
        if not s[a].index.equals(s[b].index):
            raise ValueError(f"time index mismatch between {a} and {b}")
    return s


def daylight_mask(pv_da, pv_rt):
    """Hour counts as daylight for a PV site when either forecast or realised hourly output is positive."""
    return (pv_da > 0) | (pv_rt > 0)


def excess_kurtosis(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    m = x.mean()
    v = ((x - m) ** 2).mean()
    return ((x - m) ** 4).mean() / v ** 2 - 3.0


def pairwise_corr(errors):
    """Pearson correlation of errors and of squared standardised errors, using pairwise-complete hours."""
    z = (errors - errors.mean()) / errors.std()
    return errors.corr(min_periods=200), (z ** 2).corr(min_periods=200)


def site_statistics(err, fc, cap):
    rows = []
    for site in err.columns:
        e = err[site].dropna()
        rows.append({
            "site": site,
            "capacity_MW": cap[site],
            "hours": len(e),
            "mean_MW": e.mean(),
            "std_MW": e.std(),
            "std_pct_cap": 100 * e.std() / cap[site],
            "excess_kurtosis": excess_kurtosis(e),
            "q01_MW": e.quantile(0.01),
            "q99_MW": e.quantile(0.99),
            "forecast_mean_MW": fc[site].loc[e.index].mean(),
        })
    return pd.DataFrame(rows)


def spread_by_forecast_level(err, fc, cap, n_bins=10):
    """Standard deviation of error (percent of capacity) in bins of capacity-normalised forecast."""
    edges = np.linspace(0, 1, n_bins + 1)
    rows = []
    for site in err.columns:
        e = err[site].dropna()
        level = fc[site].loc[e.index] / cap[site]
        b = np.clip(np.digitize(level, edges) - 1, 0, n_bins - 1)
        for k in range(n_bins):
            sel = b == k
            rows.append({
                "site": site,
                "bin_low": edges[k],
                "bin_high": edges[k + 1],
                "hours": int(sel.sum()),
                "std_pct_cap": 100 * e[sel].std() / cap[site] if sel.sum() > 20 else np.nan,
            })
    return pd.DataFrame(rows)


def conditional_standardise(err, fc, cap, n_bins=10):
    """Errors divided by their forecast-level-bin standard deviation, removing heteroscedasticity."""
    edges = np.linspace(0, 1, n_bins + 1)
    out = err.copy() * np.nan
    for site in err.columns:
        e = err[site].dropna()
        level = fc[site].loc[e.index] / cap[site]
        b = np.clip(np.digitize(level, edges) - 1, 0, n_bins - 1)
        z = pd.Series(np.nan, index=e.index)
        for k in range(n_bins):
            sel = b == k
            if sel.sum() > 20:
                ek = e[sel]
                z[sel] = (ek - ek.mean()) / ek.std()
        out.loc[e.index, site] = z
    return out
