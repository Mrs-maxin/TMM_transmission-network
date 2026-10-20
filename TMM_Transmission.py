"""New-study material checks and environment recording, without AC computation."""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def load_json(relative_path):
    with (ROOT / relative_path).open("r", encoding="utf-8") as handle:
        return json.load(handle)


def file_sha256(relative_path):
    digest = hashlib.sha256()
    with (ROOT / relative_path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def connected_bus_ids(network, excluded_branch_row=None):
    active = {int(row[0]) for row in network["bus"] if int(row[1]) != 4}
    graph = {bus: set() for bus in active}
    for index, row in enumerate(network["branch"]):
        if row[10] == 0 or index == excluded_branch_row:
            continue
        left, right = int(row[0]), int(row[1])
        if left in active and right in active:
            graph[left].add(right)
            graph[right].add(left)
    seen, pending = set(), [min(active)] if active else []
    while pending:
        bus = pending.pop()
        if bus in seen:
            continue
        seen.add(bus)
        pending.extend(graph[bus] - seen)
    return seen == active and bool(active)


def matching_branches(network, pair):
    requested = {int(pair[0]), int(pair[1])}
    return [
        index for index, row in enumerate(network["branch"])
        if {int(row[0]), int(row[1])} == requested and row[10] != 0
    ]


def check_materials():
    study = load_json("data/protocol/study.json")
    registry = load_json("data/protocol/case_registry.json")
    sources = load_json("data/networks/source_registry.json")
    errors, checks, networks = [], [], {}
    if study["study_id"] != "ieee14-baseline-v1":
        errors.append("Unexpected study id")
    if study["TMM_family"] != "elliptical-ordinary-multivariate-t":
        errors.append("Ordinary multivariate t family is required")
    if study["evidence"]["old_results_allowed"]:
        errors.append("Old numerical results must not be used as new evidence")
    for source in sources["networks"]:
        network_id = source["network_id"]
        network = load_json(source["json_path"])
        networks[network_id] = network
        bus_ids = [int(row[0]) for row in network["bus"]]
        ids = set(bus_ids)
        if len(ids) != len(bus_ids):
            errors.append(network_id + ": duplicate external bus ids")
        if file_sha256(source["json_path"]) != source["json_sha256"]:
            errors.append(network_id + ": input JSON hash differs")
        raw_path = source["local_source_path"]
        if file_sha256(raw_path) != source["source_sha256"]:
            errors.append(network_id + ": original network source hash differs")
        if any(len(row) < 13 for row in network["bus"]):
            errors.append(network_id + ": short bus row")
        if any(len(row) < 10 for row in network["gen"]):
            errors.append(network_id + ": fewer than 10 required gen columns")
        if any(len(row) < 13 for row in network["branch"]):
            errors.append(network_id + ": short branch row")
        online_generators = [row for row in network["gen"] if row[7] != 0]
        generator_buses = [int(row[0]) for row in online_generators]
        if len(set(generator_buses)) != len(generator_buses):
            errors.append(network_id + ": multiple online generators at a bus need a tested allocation rule")
        if any(bus not in ids for bus in generator_buses):
            errors.append(network_id + ": generator refers to an unknown bus")
        if any(int(row[0]) not in ids or int(row[1]) not in ids for row in network["branch"]):
            errors.append(network_id + ": branch refers to an unknown bus")
        if not connected_bus_ids(network):
            errors.append(network_id + ": source network is disconnected")
        checks.append({
            "network": network_id,
            "bus_count": len(network["bus"]),
            "generator_count": len(network["gen"]),
            "branch_count": len(network["branch"]),
            "source_gen_columns": len(network["gen"][0]),
            "AC_evaluated": False,
        })
    case_ids = [case["id"] for case in registry["probabilistic_cases"]]
    if len(set(case_ids)) != len(case_ids):
        errors.append("Duplicate case ids")
    for case in registry["probabilistic_cases"]:
        network = networks[case["network"]]
        if case["status"] != "planned-not-run":
            errors.append(case["id"] + ": expected an unrun case in initial materials")
        for monitor in case["monitor_pairs"]:
            rows = matching_branches(network, monitor["branch"])
            if len(rows) != 1 or monitor["end_bus"] not in monitor["branch"]:
                errors.append(case["id"] + ": invalid or ambiguous branch-end monitor")
        if case["outage"]:
            rows = matching_branches(network, case["outage"])
            if len(rows) != 1:
                errors.append(case["id"] + ": outage branch is absent or ambiguous")
            elif not connected_bus_ids(network, rows[0]):
                errors.append(case["id"] + ": registered outage would disconnect the network")
    return {
        "kind": "material-check-not-power-flow-test",
        "study_id": study["study_id"],
        "UTC_time": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if not errors else "failed",
        "errors": errors,
        "network_checks": checks,
        "case_count": len(case_ids),
        "AC_solves": 0,
        "TMM_fits": 0,
    }


def record_environment():
    packages = {}
    for name in ("numpy", "scipy", "matplotlib", "scikit-learn", "PYPOWER", "pytest"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    disk = shutil.disk_usage(ROOT)
    affinity = sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None
    return {
        "kind": "environment-probe-not-numerical-validation",
        "study_id": "ieee14-baseline-v1",
        "provider": "unverified-by-this-script",
        "tool_execution_location_requires_external_confirmation": True,
        "UTC_time": datetime.now(timezone.utc).isoformat(),
        "OS": platform.platform(),
        "machine": platform.machine(),
        "cwd": str(Path.cwd()),
        "project_root": str(ROOT),
        "python": sys.version,
        "python_executable": sys.executable,
        "CPU_threads": os.cpu_count(),
        "CPU_affinity": affinity,
        "disk_total_bytes": disk.total,
        "available_space_bytes": disk.free,
        "packages": packages,
        "research_execution_policy": "OpenAI-hosted only; no research computation on user Mac",
        "AC_solves": 0,
        "TMM_fits": 0,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--check-materials", action="store_true")
    modes.add_argument("--preflight", action="store_true")
    args = parser.parse_args(argv)
    report = check_materials() if args.check_materials else record_environment()
    folder = ROOT / "result" / "startup"
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / ("materials_check.json" if args.check_materials else "environment_probe.json")
    path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"report": str(path), **report}, ensure_ascii=False, indent=2))
    return 1 if report.get("status") == "failed" else 0


if __name__ == "__main__":
    raise SystemExit(main())
