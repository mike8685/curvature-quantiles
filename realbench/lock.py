"""Freeze (pre-registration) and verification. See PROTOCOL.md section 4."""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import subprocess
from pathlib import Path

from . import config as C
from .datasets import DATASETS, ROOT, data_root, raw_files

HERE = Path(__file__).resolve().parent
REGISTRY = HERE / "registry.json"
RUNLOG = ROOT / "results" / "realbench" / "runlog.jsonl"


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def code_files() -> list[Path]:
    files = [HERE / "PROTOCOL.md"] + sorted(HERE.glob("*.py")) + sorted((ROOT / "src" / "curvquant").glob("*.py"))
    return files


def code_hash() -> str:
    h = hashlib.sha256()
    for p in code_files():
        h.update(str(p.relative_to(ROOT)).encode())
        h.update(_sha(p).encode())
    return h.hexdigest()


def data_manifest(dataset: str) -> dict:
    root = data_root()
    files = {str(Path(f).relative_to(root)): _sha(Path(f)) for f in raw_files(dataset)}
    h = hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest()
    return {"files": files, "hash": h}


def git_head() -> str | None:
    try:
        return subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        return None


def load_registry() -> dict:
    return json.loads(REGISTRY.read_text()) if REGISTRY.exists() else {}


def tested(dataset: str) -> bool:
    if not RUNLOG.exists():
        return False
    return any(json.loads(l).get("dataset") == dataset and json.loads(l).get("split") == "test"
               for l in RUNLOG.read_text().splitlines() if l.strip())


def freeze(dataset: str, tuned: dict, confirmatory_flag: bool, amend: str | None = None) -> dict:
    reg = load_registry()
    cat = DATASETS[dataset].category
    if dataset in reg and tested(dataset) and not amend:
        raise SystemExit(f"{dataset} already has a test run under its frozen entry; re-freezing needs --amend REASON "
                         "and both versions will be reported (protocol 4).")
    if dataset in C.FIXED_CONFIRMATORY:
        confirmatory = True
    elif confirmatory_flag:
        if cat not in C.LATE_CONFIRMATORY_CATEGORIES:
            raise SystemExit(f"category '{cat}' has no late confirmatory slot")
        taken = [d for d, e in reg.items() if e.get("confirmatory") and DATASETS[d].category == cat and d != dataset]
        if taken:
            raise SystemExit(f"the confirmatory slot for '{cat}' is already taken by {taken[0]}")
        confirmatory = True
    else:
        confirmatory = False
    man = data_manifest(dataset)
    if not man["files"]:
        raise SystemExit(f"no raw files for {dataset} under {data_root()}")
    entry = {"dataset": dataset, "category": cat, "confirmatory": confirmatory,
             "frozen_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
             "protocol_version": C.PROTOCOL_VERSION, "code_hash": code_hash(), "git_head": git_head(),
             "data_hash": man["hash"], "data_files": man["files"], "tuned": tuned,
             "history": (reg.get(dataset, {}).get("history", []) +
                         ([{"replaced": reg[dataset]["frozen_at"], "reason": amend}] if dataset in reg else []))}
    reg[dataset] = entry
    REGISTRY.write_text(json.dumps(reg, indent=2))
    return entry


def verify(dataset: str) -> dict:
    reg = load_registry()
    if dataset not in reg:
        raise SystemExit(f"{dataset} is not frozen; run dev, then freeze, then push the freeze commit (protocol 4)")
    e = reg[dataset]
    problems = []
    if e["code_hash"] != code_hash():
        problems.append("protocol/harness/curvquant source changed since freeze")
    if e["data_hash"] != data_manifest(dataset)["hash"]:
        problems.append("raw data files changed since freeze")
    if problems:
        raise SystemExit("REFUSING test run: " + "; ".join(problems))
    return e


def public_check() -> list[str]:
    """Problems that would make the freeze not a public pre-registration."""
    issues = []
    try:
        dirty = subprocess.check_output(["git", "-C", str(ROOT), "status", "--porcelain", "--",
                                         "realbench", "src/curvquant"], text=True).strip()
        if dirty:
            issues.append("uncommitted changes in realbench/ or src/curvquant/")
        head = git_head()
        remote = subprocess.check_output(["git", "-C", str(ROOT), "branch", "-r", "--contains", head],
                                         text=True, stderr=subprocess.DEVNULL).strip()
        if not remote:
            issues.append("HEAD (with the frozen registry) is not on any pushed remote branch")
    except Exception as ex:
        issues.append(f"git check failed: {ex}")
    return issues


def log_run(record: dict) -> None:
    RUNLOG.parent.mkdir(parents=True, exist_ok=True)
    record = {"at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"), **record}
    with open(RUNLOG, "a") as f:
        f.write(json.dumps(record, default=str) + "\n")
