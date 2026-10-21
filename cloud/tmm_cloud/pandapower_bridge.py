"""pandapower model of a MATPOWER case: AC OPF for base dispatch and distributed-slack power flow
used as an independent reference for the batched Newton-Raphson solver."""

import warnings

import numpy as np
import pandapower as pp
from pandapower.converter.pypower import from_ppc

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)


class PandapowerModel:
    def __init__(self, case, net):
        ppc = {k: (v.copy() if hasattr(v, "copy") else v) for k, v in case.items() if k in ("bus", "gen", "branch", "gencost", "baseMVA")}
        ppc["version"] = "2"
        self.pnet = from_ppc(ppc, f_hz=60, validate_conversion=False)
        self.net = net
        p = self.pnet
        if not np.array_equal(np.sort(p.bus.index.values), np.sort(net.bus_id)):
            raise ValueError("unexpected pandapower bus labels")
        self.label = net.bus_id
        pos = {int(b): i for i, b in enumerate(net.bus_id)}
        self.gen_row = -np.ones(net.n_gen, int)
        bus_to_unit = {int(net.bus_id[b]): g for g, b in enumerate(net.gen_bus)}
        for row, b in zip(p.gen.index, p.gen.bus.values):
            self.gen_row[bus_to_unit[int(b)]] = row
        self.ext_unit = bus_to_unit[int(p.ext_grid.bus.iloc[0])]
        self.counts = {"bus": len(p.bus), "gen": len(p.gen), "ext_grid": len(p.ext_grid), "line": len(p.line),
                       "trafo": len(p.trafo), "impedance": len(p.impedance), "load_P_MW": float(p.load.p_mw.sum()),
                       "load_Q_Mvar": float(p.load.q_mvar.sum()), "shunt": len(p.shunt)}
        p.load["p_mw"] = 0.0
        p.load["q_mvar"] = 0.0
        self.load_idx = {pos[int(b)]: i for i, b in zip(p.load.index, p.load.bus.values)}
        for b in range(net.n):
            if b not in self.load_idx:
                self.load_idx[b] = pp.create_load(p, bus=int(net.bus_id[b]), p_mw=0.0, q_mvar=0.0, controllable=False)
        self.load_rows = np.array([self.load_idx[b] for b in range(net.n)])
        p.load["controllable"] = False
        self.sgen_rows = None
        # pandapower gives the ext_grid a zero active schedule under distributed slack; the reference
        # unit schedule is carried by this fixed injection so that its output is Pg0_ref + alpha_ref * lam.
        self.ref_sched = pp.create_sgen(p, bus=int(p.ext_grid.bus.iloc[0]), p_mw=0.0, q_mvar=0.0, controllable=False)

    def add_renewables(self, ren_bus):
        p = self.pnet
        self.sgen_rows = np.array([pp.create_sgen(p, bus=int(self.label[b]), p_mw=0.0, q_mvar=0.0, controllable=False) for b in ren_bus])

    def set_injections(self, Pd_MW, Qd_MVAr, Pren_MW):
        p = self.pnet
        p.load.loc[self.load_rows, "p_mw"] = Pd_MW
        p.load.loc[self.load_rows, "q_mvar"] = Qd_MVAr
        p.sgen.loc[self.sgen_rows, "p_mw"] = Pren_MW
        p.sgen.loc[self.sgen_rows, "q_mvar"] = 0.0

    def run_opf(self, headroom=0.0):
        """AC OPF with branch apparent power kept within rateA. pandapower limits branch current, so the
        loading limit is scaled by 1/max(Vmax) to keep |S| <= rateA at any admissible voltage.
        headroom: fraction of Pmax kept free above Pmin and below Pmax for units with Pmax > 0, as
        balancing reserve for the participation-factor sharing of the imbalance."""
        p = self.pnet
        net = self.net
        vmax = float(net.Vmax.max())
        p.line["max_loading_percent"] = 100.0 / vmax
        p.trafo["max_loading_percent"] = 100.0 / vmax
        for g in range(net.n_gen):
            pmax, pmin = net.Pmax[g] * net.baseMVA, net.Pmin[g] * net.baseMVA
            lo, hi = (pmin + headroom * pmax, pmax - headroom * pmax) if pmax > 0 else (pmin, pmax)
            if g == self.ext_unit:
                p.ext_grid.loc[p.ext_grid.index[0], ["min_p_mw", "max_p_mw"]] = [lo, hi]
            else:
                p.gen.loc[self.gen_row[g], ["min_p_mw", "max_p_mw"]] = [lo, hi]
        p.sgen.loc[self.ref_sched, "p_mw"] = 0.0
        p.gen["controllable"] = True
        p.ext_grid["controllable"] = True
        pp.runopp(p, init="pf", calculate_voltage_angles=True, delta=1e-10)
        net = self.net
        Pg = np.zeros(net.n_gen)
        Vg = np.zeros(net.n_gen)
        for g in range(net.n_gen):
            if g == self.ext_unit:
                Pg[g] = p.res_ext_grid.p_mw.iloc[0]
            else:
                Pg[g] = p.res_gen.p_mw.loc[self.gen_row[g]]
            Vg[g] = p.res_bus.vm_pu.loc[self.label[net.gen_bus[g]]]
        return {"Pg_MW": Pg, "Vset": Vg, "cost": float(p.res_cost), "vm": p.res_bus.vm_pu.loc[self.label].values.copy(),
                "va_deg": p.res_bus.va_degree.loc[self.label].values.copy()}

    def set_dispatch(self, Pg_MW, Vset, alpha):
        p = self.pnet
        for g in range(self.net.n_gen):
            if g == self.ext_unit:
                p.ext_grid.loc[p.ext_grid.index[0], "vm_pu"] = Vset[g]
                p.ext_grid.loc[p.ext_grid.index[0], "slack_weight"] = alpha[g]
            else:
                r = self.gen_row[g]
                p.gen.loc[r, "p_mw"] = Pg_MW[g]
                p.gen.loc[r, "vm_pu"] = Vset[g]
                p.gen.loc[r, "slack_weight"] = alpha[g]
        p.sgen.loc[self.ref_sched, "p_mw"] = Pg_MW[self.ext_unit]

    def run_pf(self, tol_mva=1e-8):
        """Distributed-slack Newton-Raphson power flow with unit reactive limits enforced."""
        p = self.pnet
        pp.runpp(p, algorithm="nr", distributed_slack=True, enforce_q_lims=True, tolerance_mva=tol_mva,
                 init="flat", max_iteration=50, calculate_voltage_angles=True)
        net = self.net
        rb = p.res_bus.loc[self.label]
        V = rb.vm_pu.values * np.exp(1j * np.deg2rad(rb.va_degree.values))
        Pg = np.zeros(net.n_gen)
        Qg = np.zeros(net.n_gen)
        for g in range(net.n_gen):
            if g == self.ext_unit:
                Pg[g] = p.res_ext_grid.p_mw.iloc[0] + p.sgen.p_mw.loc[self.ref_sched]
                Qg[g] = p.res_ext_grid.q_mvar.iloc[0]
            else:
                Pg[g] = p.res_gen.p_mw.loc[self.gen_row[g]]
                Qg[g] = p.res_gen.q_mvar.loc[self.gen_row[g]]
        # pandapower restores bus types after the run; a unit fixed at a limit reports exactly that limit.
        mode = np.zeros(net.n_gen, int)
        for g in range(net.n_gen):
            if g == self.ext_unit:
                continue
            q = Qg[g] / net.baseMVA
            if abs(q - net.Qmax[g]) < 1e-9:
                mode[g] = 1
            elif abs(q - net.Qmin[g]) < 1e-9:
                mode[g] = -1
        return {"V": V, "Pg_MW": Pg, "Qg_Mvar": Qg, "mode": mode, "converged": bool(p.converged)}
