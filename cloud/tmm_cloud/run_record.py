"""Result folders and run_manifest.json records for cloud computations."""

import datetime as _dt
import hashlib
import json
import platform
import subprocess
import sys
from importlib import metadata
from pathlib import Path

import numpy as np

from . import PROJECT_ROOT, CLOUD_ROOT

PACKAGES = ["numpy", "scipy", "pandas", "matplotlib", "pandapower", "numba"]


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def git_state():
    def run(*args):
        try:
            return subprocess.check_output(["git", *args], cwd=PROJECT_ROOT, text=True).strip()
        except Exception:
            return None

    status = run("status", "--porcelain", "--untracked-files=no")
    return {
        "branch": run("rev-parse", "--abbrev-ref", "HEAD"),
        "commit": run("rev-parse", "HEAD"),
        "dirty": bool(status),
    }


def environment():
    versions = {}
    for name in PACKAGES:
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = None
    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "packages": versions,
    }


def new_result_dir(task, stage):
    stamp = _dt.datetime.now(_dt.timezone.utc).strftime("%Y%m%dT%H%M%S")
    path = CLOUD_ROOT / "results" / task / stamp / stage
    path.mkdir(parents=True, exist_ok=True)
    return path


def latest_stage_dir(task, stage):
    base = CLOUD_ROOT / "results" / task
    candidates = sorted(p for p in base.glob(f"*/{stage}") if (p / "run_manifest.json").exists())
    if not candidates:
        raise FileNotFoundError(f"no finished result for {task}/{stage}")
    return candidates[-1]


def relpath(path):
    return str(Path(path).resolve().relative_to(PROJECT_ROOT))


def _default(obj):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, Path):
        return relpath(obj)
    raise TypeError(f"not serialisable: {type(obj)}")


def write_json(path, obj):
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, ensure_ascii=False, indent=2, default=_default)


def write_manifest(result_dir, command, data_files, parameters, random_settings, outputs=None, extra=None):
    manifest = {
        "created_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "command": command,
        "git": git_state(),
        "environment": environment(),
        "data_sources": [{"path": relpath(p), "sha256": sha256_file(p)} for p in data_files],
        "parameters": parameters,
        "random_sampling": random_settings,
        "outputs": sorted(relpath(p) for p in Path(result_dir).rglob("*") if p.is_file() and p.name != "run_manifest.json")
        if outputs is None else outputs,
    }
    if extra:
        manifest.update(extra)
    write_json(Path(result_dir) / "run_manifest.json", manifest)
    return manifest


def refresh_outputs(result_dir):
    """Re-list output files in an existing manifest, e.g. after README.md is written."""
    import json

    path = Path(result_dir) / "run_manifest.json"
    manifest = json.loads(path.read_text(encoding="utf-8"))
    manifest["outputs"] = sorted(relpath(p) for p in Path(result_dir).rglob("*") if p.is_file() and p.name != "run_manifest.json")
    write_json(path, manifest)
