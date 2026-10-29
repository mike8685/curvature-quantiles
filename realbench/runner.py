"""Run methods on curves and record every error (one row per curve x K x method x smoothing)."""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import os
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

from curvquant.core import turning_angles

from . import config as C
from .datasets import ROOT, load_split
from .geometry import arclength, diameter, gaussian_smooth, index_to_s, locate, points_at_index
from .methods import ABLATIONS, METHODS, VERTEX_SUBSET, vertex_indices
from .metrics import errors

RESULTS = ROOT / "results" / "realbench"


def budgets(N: int) -> list[tuple[float, int]]:
    out, seen = [], set()
    for r in C.BUDGET_RATIOS:
        K = max(C.K_MIN, int(round(r * N)))
        if K > C.K_MAX_FRACTION * N or K in seen:
            continue
        seen.add(K)
        out.append((r, K))
    return out


def curve_features(pts: np.ndarray) -> dict:
    lengths, w = turning_angles(pts)
    dual = 0.5 * (lengths[:-1] + lengths[1:])
    return {"total_turning": float(w.sum()),
            "straight_frac": float(dual[w < np.radians(0.5)].sum() / max(dual.sum(), 1e-300))}


def _uniform_frac(dec: np.ndarray, K: int) -> np.ndarray:
    """Fallback (protocol 5.4): K points at equal arc length, as fractional vertex indices."""
    s = arclength(dec)
    return np.interp(np.linspace(0.0, s[-1], K), s, np.arange(len(dec), dtype=np.float64))


def run_one(pts: np.ndarray, closed: bool, method: str, K: int, c: int, smoothed: np.ndarray,
            vertices_only: bool = False) -> dict:
    """Run one method at one smoothing level; return metrics and bookkeeping.

    The method decides on ``smoothed`` (== pts when c = 0); its output is converted to
    fractional vertex indices and placed on the ORIGINAL polyline. Any exception, in the
    method or in recovering positions, is a recorded failure and falls back to uniform.
    """
    dec = smoothed
    failed, reason = False, ""
    t0 = time.perf_counter()
    try:
        out_dec = METHODS[method](dec, K, closed)
        ms = 1e3 * (time.perf_counter() - t0)
        if method in VERTEX_SUBSET:
            frac = vertex_indices(dec, out_dec)
        else:
            frac = locate(dec, out_dec, C.ONCURVE_TOL * diameter(dec))
        out = out_dec if c == 0 else points_at_index(pts, frac)
    except Exception as ex:
        ms = 1e3 * (time.perf_counter() - t0)
        failed, reason = True, f"{type(ex).__name__}: {ex}"[:200]
        frac = _uniform_frac(dec, K)
        out = points_at_index(pts, frac)
        if closed:
            out[-1] = out[0]
    out_s = index_to_s(pts, frac)
    e = errors(pts, out, out_s, vertices_only=vertices_only)
    pin_active = np.nan
    if method in ("cq", "cq_snapped"):
        _, w = turning_angles(dec)
        ncorner = int((w > C.CQ_CORNER_THRESHOLD).sum())
        pin_active = float(K > 4 and 0 < ncorner <= K - 2)
    return {"failed": failed, "fail_reason": reason, "n_out": len(out), "ms": ms, "pin_active": pin_active,
            **e, "_out": out}


def plan(split: str, tuned: dict | None) -> list[tuple[str, int]]:
    if split == "dev":
        mc = [(m, c) for m in C.TUNED_METHODS for c in C.SMOOTHING_GRID]
    else:
        mc = [(m, 0) for m in C.TUNED_METHODS]
        mc += [(m, int(tuned[m])) for m in C.TUNED_METHODS if int(tuned[m]) != 0]
    mc += [(m, 0) for m in ABLATIONS] + [("oracle", 0)]
    return mc


def _work(args):
    cv, mc = args
    pts = cv.pts
    N = len(pts)
    feats = curve_features(pts)
    diam = diameter(pts)
    cache = {}
    rows = []
    for r, K in budgets(N):
        for m, c in mc:
            if m == "oracle" and N > C.ORACLE_MAX_N:
                continue
            if c not in cache:
                cache[c] = gaussian_smooth(pts, c, closed=cv.closed)
            res = run_one(pts, cv.closed, m, K, c, cache[c])
            res.pop("_out")
            rows.append({"dataset": cv.dataset, "curve_id": cv.curve_id, "cluster": cv.cluster, "N": N,
                         "closed": cv.closed, "diam": diam, **feats, "r": r, "K": K, "method": m, "c": c, **res})
    return rows


def sha_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _compute(dataset, split, mc, procs, cap):
    stats: dict = {}
    curves = load_split(dataset, split, cap=cap, stats=stats)
    stats.pop("_cluster_ids", None)
    procs = procs or min(os.cpu_count() or 1, 4)
    rows = []
    if procs > 1 and len(curves) > 1:
        with Pool(procs) as pool:
            for rr in pool.imap_unordered(_work, [(cv, mc) for cv in curves], chunksize=1):
                rows.extend(rr)
    else:
        for cv in curves:
            rows.extend(_work((cv, mc)))
    df = pd.DataFrame(rows)
    if len(df):
        df = df.sort_values(["curve_id", "K", "method", "c"]).reset_index(drop=True)
    return df, curves, stats


def run_dev(dataset: str, procs: int | None = None, cap: int | None = None) -> tuple[pd.DataFrame, dict]:
    from .lock import code_hash, env_versions
    t0 = time.time()
    df, curves, stats = _compute(dataset, "dev", plan("dev", None), procs, cap)
    out_dir = RESULTS / "dev" / dataset
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / "rows.csv", index=False)
    tuned = tune(df)
    tuned["rows_sha256"] = sha_file(out_dir / "rows.csv")
    tuned["code_hash"] = code_hash()
    (out_dir / "tuned.json").write_text(json.dumps(tuned, indent=2))
    meta = {"dataset": dataset, "split": "dev", "n_curves": len(curves), "load_stats": stats,
            "plan": plan("dev", None), "tuned": tuned["chosen"], "code_hash": tuned["code_hash"],
            "versions": env_versions(), "seconds": round(time.time() - t0, 1)}
    (out_dir / "meta.json").write_text(json.dumps(meta, indent=2, default=str))
    return df, meta


def run_test(dataset: str, procs: int | None = None, allow_unpushed: bool = False) -> tuple[pd.DataFrame, dict]:
    """The only way to compute on the test split: verify the freeze, log the start
    (committed and pushed unless allow_unpushed), compute into a fresh run directory,
    log the result hash. Earlier runs are never overwritten."""
    from .lock import log_run, publish_runlog, public_check, verify
    entry = verify(dataset)
    issues = public_check()
    if issues and not allow_unpushed:
        raise SystemExit("REFUSING test run: " + "; ".join(issues) + " (protocol 4, step 3)")
    run_id = dt.datetime.now(dt.timezone.utc).strftime("run-%Y%m%dT%H%M%SZ")
    log_run({"event": "start", "dataset": dataset, "split": "test", "run_id": run_id,
             "code_hash": entry["code_hash"], "data_hash": entry["data_hash"], "tuned": entry["tuned"],
             "confirmatory": entry["confirmatory"], "public_issues": issues, "allow_unpushed": allow_unpushed})
    if not allow_unpushed:
        publish_runlog(f"realbench: start test run {dataset} {run_id}")
    t0 = time.time()
    mc = plan("test", entry["tuned"])
    df, curves, stats = _compute(dataset, "test", mc, procs, None)
    out_dir = RESULTS / "test" / dataset / run_id
    out_dir.mkdir(parents=True, exist_ok=False)
    df.to_csv(out_dir / "rows.csv", index=False)
    meta = {"dataset": dataset, "split": "test", "run_id": run_id, "n_curves": len(curves), "load_stats": stats,
            "plan": mc, "tuned": entry["tuned"], "confirmatory": entry["confirmatory"],
            "code_hash": entry["code_hash"], "seconds": round(time.time() - t0, 1)}
    (out_dir / "meta.json").write_text(json.dumps(meta, indent=2, default=str))
    log_run({"event": "done", "dataset": dataset, "split": "test", "run_id": run_id, "n_curves": len(curves),
             "rows": len(df), "rows_sha256": sha_file(out_dir / "rows.csv"), "seconds": meta["seconds"]})
    return df, meta


def tune(df: pd.DataFrame) -> dict:
    """Pick the smoothing level per method on dev (protocol 2.1): lowest mean cluster log-ratio vs c=0."""
    from .stats import cluster_scores, log_ratios
    chosen, table = {}, {}
    for m in C.TUNED_METHODS:
        scores = {}
        for c in C.SMOOTHING_GRID:
            s = cluster_scores(log_ratios(df, m, m, c, 0, C.PRIMARY_METRIC))
            scores[c] = float(s.mean()) if len(s) else 0.0
        chosen[m] = int(min(C.SMOOTHING_GRID, key=lambda c: (round(scores[c], 9), c)))
        table[m] = scores
    return {"chosen": chosen, "mean_cluster_log_ratio_vs_c0": table}
