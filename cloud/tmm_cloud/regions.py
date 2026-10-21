"""Three-region partition of a transmission network and placement of RTS-GMLC renewable sites."""

import numpy as np
import pandas as pd
from scipy.cluster.vq import kmeans2
from scipy.sparse.csgraph import connected_components

from . import rts_gmlc as rts


def spectral_regions(net, k=3, restarts=50, rng_seed=20261021):
    """Normalised spectral clustering (Ng-Jordan-Weiss) on the branch series-admittance magnitude graph.

    Returns region labels ordered so that region 1 holds the lowest bus number, region 2 the next.
    """
    A = net.adjacency()
    d = A.sum(axis=1)
    Dm = 1.0 / np.sqrt(d)
    L = np.eye(net.n) - Dm[:, None] * A * Dm[None, :]
    w, U = np.linalg.eigh(L)
    X = U[:, :k]
    X = X / np.linalg.norm(X, axis=1, keepdims=True)
    rng = np.random.default_rng(rng_seed)
    best = None
    for _ in range(restarts):
        cent, lab = kmeans2(X, k, minit="++", seed=rng)
        inertia = ((X - cent[lab]) ** 2).sum()
        if best is None or inertia < best[0] - 1e-12:
            best = (inertia, lab)
    lab = best[1]
    order = np.argsort([net.bus_id[lab == c].min() for c in range(k)])
    remap = np.empty(k, int)
    remap[order] = np.arange(k)
    regions = remap[lab] + 1
    return regions, w[:k + 1]


def region_connectivity(net, regions):
    out = {}
    A = net.adjacency() > 0
    for r in np.unique(regions):
        m = regions == r
        nc, _ = connected_components(A[np.ix_(m, m)], directed=False)
        out[int(r)] = int(nc)
    return out


def tie_lines(net, regions):
    rf, rt = regions[net.f], regions[net.t]
    return np.flatnonzero(rf != rt)


def pv_groups(n_groups_by_area):
    """Group RTS-GMLC UPV plants of each RTS area by k-means on plant bus latitude/longitude."""
    gen = rts.generator_table()
    bus = rts.bus_table().set_index("Bus ID")
    pv = gen[gen["Unit Type"] == "PV"][["GEN UID", "Bus ID", "PMax MW"]].copy()
    pv["area"] = pv["Bus ID"] // 100
    pv["lat"] = bus.loc[pv["Bus ID"], "lat"].values
    pv["lng"] = bus.loc[pv["Bus ID"], "lng"].values
    pv["group"] = ""
    rng = np.random.default_rng(7)
    for area, ng in n_groups_by_area.items():
        sel = pv.area == area
        pts = pv.loc[sel, ["lat", "lng"]].values
        if ng == 1 or len(np.unique(pts, axis=0)) <= 1:
            lab = np.zeros(sel.sum(), int)
        else:
            best = None
            for _ in range(30):
                cent, lab_ = kmeans2(pts, ng, minit="++", seed=rng)
                if len(np.unique(lab_)) < ng:
                    continue
                inertia = ((pts - cent[lab_]) ** 2).sum()
                if best is None or inertia < best[0] - 1e-12:
                    best = (inertia, lab_)
            lab = best[1]
            order = np.argsort([-pv.loc[sel, "PMax MW"].values[lab == c].sum() for c in range(ng)])
            remap = np.empty(ng, int)
            remap[order] = np.arange(ng)
            lab = remap[lab]
        pv.loc[sel, "group"] = [f"UPV_A{area}_G{g + 1}" for g in lab]
    return pv


def host_candidates(net, regions, base_load_MW):
    """PQ buses with positive load, at least two branches and their connected rateA sum."""
    is_gen = np.zeros(net.n, bool)
    is_gen[net.gen_bus] = True
    deg = net.degree()
    cap = np.zeros(net.n)
    np.add.at(cap, net.f, net.rateA)
    np.add.at(cap, net.t, net.rateA)
    df = pd.DataFrame({
        "bus": net.bus_id, "idx": np.arange(net.n), "region": regions, "degree": deg,
        "is_gen_bus": is_gen, "load_MW": base_load_MW, "connected_rateA_MVA": cap,
    })
    return df[(~df.is_gen_bus) & (df.load_MW > 0) & (df.degree >= 2)].sort_values(
        ["region", "connected_rateA_MVA"], ascending=[True, False])
