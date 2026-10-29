"""Run methods on curves and record every error (one row per curve x K x method x smoothing)."""

from __future__ import annotations

import json
import os
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

from curvquant import uniform_arclength
from curvquant.core import turning_angles

from . import config as C
from .datasets import ROOT, Curve, load_split
from .geometry import diameter, gaussian_smooth, index_to_s, locate, points_at_index
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


def run_one(pts: np.ndarray, closed: bool, method: str, K: int, c: int, smoothed: np.ndarray,
            vertices_only: bool = False) -> dict:
    """Run one method at one smoothing level; return metrics and bookkeeping."""
    tol = C.ONCURVE_TOL * diameter(pts)
    dec = smoothed
    failed = False
    t0 = time.perf_counter()
    try:
        out_dec = METHODS[method](dec, K, closed)
    except Exception:
        failed = True
        out_dec = uniform_arclength(dec, K, closed)
    ms = 1e3 * (time.perf_counter() - t0)
    if method in VERTEX_SUBSET and not failed:
        frac = vertex_indices(dec, out_dec)
    else:
        frac = locate(dec, out_dec, C.ONCURVE_TOL * diameter(dec))
    out = out_dec if c == 0 else points_at_index(pts, frac)
    out_s = index_to_s(pts, frac)
    e = errors(pts, out, out_s, vertices_only=vertices_only)
    return {"failed": failed, "n_out": len(out), "ms": ms, **e, "_out": out, "_tol": tol}


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
                cache[c] = gaussian_smooth(pts, c)
            res = run_one(pts, cv.closed, m, K, c, cache[c])
            res.pop("_out"); res.pop("_tol")
            rows.append({"dataset": cv.dataset, "curve_id": cv.curve_id, "cluster": cv.cluster, "N": N,
                         "closed": cv.closed, "diam": diam, **feats, "r": r, "K": K, "method": m, "c": c, **res})
    return rows


def run_split(dataset: str, split: str, tuned: dict | None = None, procs: int | None = None,
              cap: int | None = None, out_dir: Path | None = None) -> tuple[pd.DataFrame, dict]:
    stats: dict = {}
    curves = load_split(dataset, split, cap=cap, stats=stats)
    mc = plan(split, tuned)
    procs = procs or min(os.cpu_count() or 1, 4)
    t0 = time.time()
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
    meta = {"dataset": dataset, "split": split, "n_curves": len(curves), "load_stats": stats,
            "plan": mc, "tuned": tuned, "seconds": round(time.time() - t0, 1)}
    out_dir = out_dir or RESULTS / split / dataset
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / "rows.csv", index=False)
    (out_dir / "meta.json").write_text(json.dumps(meta, indent=2, default=str))
    return df, meta


def tune(df: pd.DataFrame) -> dict:
    """Pick the smoothing level per method on dev (protocol 2.1): lowest mean cluster log-ratio vs c=0."""
    from .stats import cluster_scores, log_ratios
    chosen, table = {}, {}
    for m in C.TUNED_METHODS:
        scores = {}
        for c in C.SMOOTHING_GRID:
            lr = log_ratios(df, m, m, c, 0, C.PRIMARY_METRIC)
            s = cluster_scores(lr)
            scores[c] = float(s.mean()) if len(s) else 0.0
        best = min(C.SMOOTHING_GRID, key=lambda c: (round(scores[c], 9), c))
        chosen[m] = int(best)
        table[m] = scores
    return {"chosen": chosen, "mean_cluster_log_ratio_vs_c0": table}
