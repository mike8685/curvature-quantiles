"""Error metrics against the ORIGINAL polyline (protocol 5.2)."""

from __future__ import annotations

import numpy as np
import shapely

from . import config as C
from .geometry import arclength, eval_samples


def _segments(poly: np.ndarray):
    return shapely.linestrings(np.stack([poly[:-1], poly[1:]], axis=1))


def dist_to_polyline(points: np.ndarray, poly: np.ndarray) -> np.ndarray:
    """Exact Euclidean distance from each point to the nearest segment of ``poly`` (STRtree)."""
    segs = _segments(poly)
    keep = ~shapely.is_empty(segs)
    tree = shapely.STRtree(segs[keep])
    idx, d = tree.query_nearest(shapely.points(points), return_distance=True, all_matches=False)
    out = np.full(len(points), np.nan)
    out[idx[0]] = d
    if np.isnan(out).any():
        raise RuntimeError("nearest-segment query returned no match for some points")
    return out


def _seg_dist(p, a, b):
    ab = b - a
    ab2 = np.maximum((ab ** 2).sum(1), 1e-300)
    t = np.clip(((p - a) * ab).sum(1) / ab2, 0.0, 1.0)
    return np.hypot(*(p - a - t[:, None] * ab).T)


def n_eval(N: int) -> int:
    return int(min(max(C.EVAL_SAMPLES_FACTOR * N, C.EVAL_SAMPLES_MIN), C.EVAL_SAMPLES_MAX))


def errors(ref: np.ndarray, out: np.ndarray, out_s: np.ndarray, vertices_only: bool = False) -> dict:
    """All metrics for one simplification.

    ref    original polyline (N, 2)
    out    simplified polyline (points on ref, in order)
    out_s  arc-length positions of ``out`` on ref (for the shortcut metric)
    """
    if vertices_only:
        E, Es, mask = ref, arclength(ref), np.zeros(len(ref), bool)
    else:
        E, Es, mask = eval_samples(ref, n_eval(len(ref)))
    d = dist_to_polyline(E, out)
    H = float(d.max())
    rmse = float(np.sqrt(np.mean(d[mask] ** 2))) if mask.any() else float(np.sqrt(np.mean(d ** 2)))

    # reverse direction: simplified polyline (vertices + uniform samples) to ref
    so, _, _ = eval_samples(out, n_eval(len(out)))
    H_rev = float(dist_to_polyline(so, ref).max())

    # shortcut: distance to the simplified segment spanning the same arc length
    S = np.maximum.accumulate(np.asarray(out_s, dtype=np.float64))
    k = np.clip(np.searchsorted(S, Es, side="right") - 1, 0, len(out) - 2)
    shortcut = float(_seg_dist(E, out[k], out[k + 1]).max())

    return {"hausdorff": H, "rmse": rmse, "hausdorff_sym": max(H, H_rev), "shortcut": shortcut}
