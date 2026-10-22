"""Read MATPOWER version-2 case files (.m) into numpy arrays."""

import re
from pathlib import Path

import numpy as np

from . import DATA_ROOT

PGLIB118 = DATA_ROOT / "networks" / "pglib118" / "pglib_opf_case118_ieee.m"

BUS_COLS = ["bus_i", "type", "Pd", "Qd", "Gs", "Bs", "area", "Vm", "Va", "baseKV", "zone", "Vmax", "Vmin"]
GEN_COLS = ["bus", "Pg", "Qg", "Qmax", "Qmin", "Vg", "mBase", "status", "Pmax", "Pmin"]
BRANCH_COLS = ["fbus", "tbus", "r", "x", "b", "rateA", "rateB", "rateC", "ratio", "angle", "status", "angmin", "angmax"]


def _matrix(text, name):
    m = re.search(rf"mpc\.{name}\s*=\s*\[(.*?)\];", text, re.S)
    if m is None:
        raise ValueError(f"mpc.{name} not found")
    rows = []
    for line in m.group(1).splitlines():
        line = line.split("%", 1)[0].strip().rstrip(";").strip()
        if line:
            rows.append([float(v) for v in line.replace(",", " ").split()])
    width = max(len(r) for r in rows)
    if any(len(r) != width for r in rows):
        raise ValueError(f"ragged rows in mpc.{name}")
    return np.array(rows)


def load_case(path=PGLIB118):
    text = Path(path).read_text()
    base = float(re.search(r"mpc\.baseMVA\s*=\s*([0-9.eE+-]+)", text).group(1))
    case = {
        "path": Path(path),
        "baseMVA": base,
        "bus": _matrix(text, "bus"),
        "gen": _matrix(text, "gen"),
        "branch": _matrix(text, "branch"),
    }
    if "mpc.gencost" in text:
        case["gencost"] = _matrix(text, "gencost")
    return case


def col(name, cols):
    return cols.index(name)
