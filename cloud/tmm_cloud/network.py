"""Bus-branch network model built from a MATPOWER case (MATPOWER makeYbus conventions)."""

from dataclasses import dataclass, field

import numpy as np

from .matpower_case import BUS_COLS, GEN_COLS, BRANCH_COLS, load_case


def _c(name, cols):
    return cols.index(name)


@dataclass
class Network:
    baseMVA: float
    bus_id: np.ndarray          # original bus numbers
    Pd: np.ndarray              # p.u.
    Qd: np.ndarray
    Vmax: np.ndarray
    Vmin: np.ndarray
    ref: int                    # internal index of reference bus
    gen_bus: np.ndarray         # internal bus index per unit
    Pg: np.ndarray              # p.u., file dispatch
    Qmax: np.ndarray
    Qmin: np.ndarray
    Pmax: np.ndarray
    Pmin: np.ndarray
    Vg: np.ndarray
    f: np.ndarray
    t: np.ndarray
    Yff: np.ndarray
    Yft: np.ndarray
    Ytf: np.ndarray
    Ytt: np.ndarray
    Ysh: np.ndarray
    rateA: np.ndarray           # MVA
    gencost: np.ndarray = None
    Y: np.ndarray = field(default=None, repr=False)

    @property
    def n(self):
        return len(self.bus_id)

    @property
    def n_gen(self):
        return len(self.gen_bus)

    def idx(self, bus_numbers):
        lookup = {b: i for i, b in enumerate(self.bus_id)}
        return np.array([lookup[int(b)] for b in np.atleast_1d(bus_numbers)])

    def build_ybus(self):
        n = self.n
        Y = np.zeros((n, n), complex)
        np.add.at(Y, (self.f, self.f), self.Yff)
        np.add.at(Y, (self.f, self.t), self.Yft)
        np.add.at(Y, (self.t, self.f), self.Ytf)
        np.add.at(Y, (self.t, self.t), self.Ytt)
        Y[np.arange(n), np.arange(n)] += self.Ysh
        self.Y = Y
        return Y

    def branch_flows(self, V):
        """Branch complex power at from and to ends in MVA; V has shape (..., n)."""
        Vf = V[..., self.f]
        Vt = V[..., self.t]
        Sf = Vf * np.conj(self.Yff * Vf + self.Yft * Vt) * self.baseMVA
        St = Vt * np.conj(self.Ytf * Vf + self.Ytt * Vt) * self.baseMVA
        return Sf, St

    def adjacency(self):
        A = np.zeros((self.n, self.n))
        w = np.abs(self.Yft)
        np.add.at(A, (self.f, self.t), w)
        np.add.at(A, (self.t, self.f), w)
        return A

    def degree(self):
        d = np.zeros(self.n, int)
        np.add.at(d, self.f, 1)
        np.add.at(d, self.t, 1)
        return d


def from_matpower(case):
    if isinstance(case, (str, bytes)) or hasattr(case, "suffix"):
        case = load_case(case)
    base = case["baseMVA"]
    bus, gen, br = case["bus"], case["gen"], case["branch"]
    bc = lambda k: bus[:, _c(k, BUS_COLS)]
    gc = lambda k: gen[:, _c(k, GEN_COLS)]
    rc = lambda k: br[:, _c(k, BRANCH_COLS)]

    bus_id = bc("bus_i").astype(int)
    lookup = {b: i for i, b in enumerate(bus_id)}
    on_gen = gc("status") > 0
    on_br = rc("status") > 0
    gen = gen[on_gen]
    br = br[on_br]

    ref = np.flatnonzero(bc("type") == 3)
    if len(ref) != 1:
        raise ValueError("expected exactly one reference bus")

    r, x, b = rc("r"), rc("x"), rc("b")
    ratio = rc("ratio")
    tap = np.where(ratio == 0, 1.0, ratio) * np.exp(1j * np.deg2rad(rc("angle")))
    Ys = 1.0 / (r + 1j * x)
    Ytt = Ys + 1j * b / 2
    Yff = Ytt / (tap * np.conj(tap))
    Yft = -Ys / np.conj(tap)
    Ytf = -Ys / tap

    net = Network(
        baseMVA=base,
        bus_id=bus_id,
        Pd=bc("Pd") / base,
        Qd=bc("Qd") / base,
        Vmax=bc("Vmax"),
        Vmin=bc("Vmin"),
        ref=int(ref[0]),
        gen_bus=np.array([lookup[int(v)] for v in gen[:, _c("bus", GEN_COLS)]]),
        Pg=gen[:, _c("Pg", GEN_COLS)] / base,
        Qmax=gen[:, _c("Qmax", GEN_COLS)] / base,
        Qmin=gen[:, _c("Qmin", GEN_COLS)] / base,
        Pmax=gen[:, _c("Pmax", GEN_COLS)] / base,
        Pmin=gen[:, _c("Pmin", GEN_COLS)] / base,
        Vg=gen[:, _c("Vg", GEN_COLS)],
        f=np.array([lookup[int(v)] for v in br[:, _c("fbus", BRANCH_COLS)]]),
        t=np.array([lookup[int(v)] for v in br[:, _c("tbus", BRANCH_COLS)]]),
        Yff=Yff, Yft=Yft, Ytf=Ytf, Ytt=Ytt,
        Ysh=(bc("Gs") + 1j * bc("Bs")) / base,
        rateA=br[:, _c("rateA", BRANCH_COLS)],
        gencost=case.get("gencost")[on_gen] if case.get("gencost") is not None else None,
    )
    if len(np.unique(net.gen_bus)) != net.n_gen:
        raise ValueError("solver assumes at most one unit per bus")
    net.build_ybus()
    return net
