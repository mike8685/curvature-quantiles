"""Freeze (pre-registration) and verification. See PROTOCOL.md section 4."""

from __future__ import annotations

import datetime as dt
import hashlib
import importlib
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
    import curvquant
    src = Path(curvquant.__file__).resolve().parent
    if src != (ROOT / "src" / "curvquant").resolve():
        raise SystemExit(f"curvquant is imported from {src}, not from this repository's src/; "
                         "install with `pip install -e .` so the frozen code is the code that runs")
    return [HERE / "PROTOCOL.md"] + sorted(HERE.glob("*.py")) + sorted(src.glob("*.py"))


def code_hash() -> str:
    h = hashlib.sha256()
    for p in code_files():
        h.update(str(p.resolve().relative_to(ROOT)).encode())
        h.update(_sha(p).encode())
    return h.hexdigest()


def env_versions() -> dict:
    out = {}
    for m in C.TRACKED_PACKAGES:
        try:
            out[m] = importlib.import_module(m).__version__
        except Exception:
            out[m] = None
    return out


def data_manifest(dataset: str) -> dict:
    root = data_root()
    files = {str(Path(f).relative_to(root)): _sha(Path(f)) for f in raw_files(dataset)}
    h = hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest()
    return {"files": files, "hash": h}


def _git(*args, check=True) -> str:
    return subprocess.check_output(["git", "-C", str(ROOT), *args], text=True,
                                   stderr=subprocess.DEVNULL if not check else None).strip()


def git_head() -> str | None:
    try:
        return _git("rev-parse", "HEAD")
    except Exception:
        return None


def load_registry() -> dict:
    return json.loads(REGISTRY.read_text()) if REGISTRY.exists() else {}


def runlog() -> list[dict]:
    if not RUNLOG.exists():
        return []
    return [json.loads(l) for l in RUNLOG.read_text().splitlines() if l.strip()]


def tested(dataset: str) -> bool:
    """True if any test run was started (run log) or any test output exists (results dir)."""
    started = any(r.get("dataset") == dataset and r.get("split") == "test" for r in runlog())
    outputs = (ROOT / "results" / "realbench" / "test" / dataset).exists()
    return started or outputs


def freeze(dataset: str, amend: str | None = None) -> dict:
    """Pre-register one dataset. Confirmatory status comes from config (hashed), never a flag."""
    from .runner import RESULTS
    reg = load_registry()
    if dataset in reg and tested(dataset) and not amend:
        raise SystemExit(f"{dataset} already has a test run; re-freezing needs --amend REASON and both "
                         "versions will be reported (protocol 4).")
    tf = RESULTS / "dev" / dataset / "tuned.json"
    rows = RESULTS / "dev" / dataset / "rows.csv"
    if not tf.exists() or not rows.exists():
        raise SystemExit(f"run `python -m realbench dev {dataset}` first")
    tuned = json.loads(tf.read_text())
    if tuned.get("rows_sha256") != _sha(rows):
        raise SystemExit("dev rows.csv does not match tuned.json; rerun the dev split")
    if tuned.get("code_hash") != code_hash():
        raise SystemExit("the dev split was run with different code; rerun `python -m realbench dev "
                         f"{dataset}` with the current code before freezing")
    man = data_manifest(dataset)
    if not man["files"]:
        raise SystemExit(f"no raw files for {dataset} under {data_root()}")
    entry = {"dataset": dataset, "category": DATASETS[dataset].category,
             "confirmatory": dataset in C.CONFIRMATORY_DATASETS,
             "frozen_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
             "protocol_version": C.PROTOCOL_VERSION, "code_hash": code_hash(), "git_head_at_freeze": git_head(),
             "versions": env_versions(), "data_root": str(data_root()), "data_hash": man["hash"],
             "data_files": man["files"], "tuned": tuned["chosen"], "dev_rows_sha256": tuned["rows_sha256"],
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
    if e.get("versions") != env_versions():
        problems.append(f"package versions differ from freeze ({e.get('versions')} vs {env_versions()})")
    if e["confirmatory"] != (dataset in C.CONFIRMATORY_DATASETS):
        problems.append("confirmatory status differs from config")
    if problems:
        raise SystemExit("REFUSING test run: " + "; ".join(problems))
    return e


def public_check() -> list[str]:
    """Problems that would make the freeze not a public pre-registration.

    Checks the real remote (git ls-remote), not local tracking refs. Whether the remote
    repository is public cannot be checked offline; its URL is recorded in the run log.
    """
    issues = []
    try:
        if _git("status", "--porcelain", "--", "realbench", "src/curvquant"):
            issues.append("uncommitted changes in realbench/ or src/curvquant/")
        head = git_head()
        url = _git("remote", "get-url", "origin")
        if "github.com" not in url:
            issues.append(f"origin is not a GitHub remote ({url})")
        tips = {l.split()[0] for l in _git("ls-remote", "origin").splitlines() if l.strip()}
        if head not in tips:
            issues.append("HEAD (with the frozen registry) is not the tip of any branch on origin")
    except Exception as ex:
        issues.append(f"git check failed: {ex}")
    return issues


def log_run(record: dict) -> None:
    RUNLOG.parent.mkdir(parents=True, exist_ok=True)
    try:
        record = {"remote": _git("remote", "get-url", "origin"), **record}
    except Exception:
        pass
    record = {"at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"), **record}
    with open(RUNLOG, "a") as f:
        f.write(json.dumps(record, default=str) + "\n")


def publish_runlog(message: str) -> None:
    """Commit and push the run log so the start of a test run is public before any computation."""
    try:
        _git("add", str(RUNLOG.relative_to(ROOT)))
        _git("commit", "-q", "-m", message)
        _git("push", "-q", "origin", "HEAD")
    except Exception as ex:
        raise SystemExit(f"could not commit and push the run-log start record ({ex}); nothing was computed")
