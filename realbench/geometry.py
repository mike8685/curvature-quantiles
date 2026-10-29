"""Geometric helpers: arc length, smoothing, locating points on a polyline, evaluation samples."""

from __future__ import annotations

import numpy as np


def arclength(pts: np.ndarray) -> np.ndarray:
    return np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(pts, axis=0).T))])


def diameter(pts: np.ndarray) -> float:
    return float(np.ptp(pts, axis=0).max())


def dedupe(pts: np.ndarray, rel_tol: float = 1e-12) -> np.ndarray:
    """Drop consecutive duplicate points (tolerance relative to the curve's extent)."""
    pts = np.asarray(pts, dtype=np.float64)
    if len(pts) < 2:
        return pts
    tol = rel_tol * max(diameter(pts), 1e-300)
    keep = np.r_[True, np.hypot(*np.diff(pts, axis=0).T) > tol]
    return pts[keep]


def is_closed(pts: np.ndarray, rel_tol: float = 1e-9) -> bool:
    return bool(np.hypot(*(pts[0] - pts[-1])) <= rel_tol * diameter(pts))


def gaussian_smooth(pts: np.ndarray, c: float) -> np.ndarray:
    """Gaussian smoothing in arc length, sigma = c * median edge length; endpoints fixed.

    Each vertex becomes the weighted mean of vertices within 4 sigma of arc length,
    weights exp(-ds^2 / 2 sigma^2) times the vertex's dual length (so irregular
    sampling does not bias the average toward dense clusters). Vertex count and
    order are preserved, so vertex i of the output corresponds to vertex i of the input.
    """
    pts = np.asarray(pts, dtype=np.float64)
    if c <= 0:
        return pts.copy()
    s = arclength(pts)
    edges = np.diff(s)
    sigma = c * float(np.median(edges))
    if sigma <= 0:
        return pts.copy()
    dual = np.empty(len(pts))
    dual[0], dual[-1] = edges[0] / 2, edges[-1] / 2
    dual[1:-1] = 0.5 * (edges[:-1] + edges[1:])
    out = np.empty_like(pts)
    lo = np.searchsorted(s, s - 4 * sigma, side="left")
    hi = np.searchsorted(s, s + 4 * sigma, side="right")
    for i in range(len(pts)):
        a, b = lo[i], hi[i]
        w = np.exp(-0.5 * ((s[a:b] - s[i]) / sigma) ** 2) * dual[a:b]
        out[i] = (w[:, None] * pts[a:b]).sum(0) / w.sum()
    out[0], out[-1] = pts[0], pts[-1]
    return out


def locate(poly: np.ndarray, q: np.ndarray, tol: float) -> np.ndarray:
    """Fractional vertex indices of points q that lie on polyline ``poly``, in order.

    Points must be ordered along the curve. For each point the search continues
    forward from the previous point's segment and takes the first segment within
    ``tol``. Raises if a point is not on the curve.
    """
    A, B = poly[:-1], poly[1:]
    AB = B - A
    ab2 = np.maximum((AB ** 2).sum(1), 1e-300)
    out = np.empty(len(q))
    j = 0
    nseg = len(A)
    for k, x in enumerate(q):
        found = False
        # vectorized look-ahead in blocks
        start = j
        while start < nseg:
            stop = min(nseg, start + 256)
            AP = x - A[start:stop]
            t = np.clip((AP * AB[start:stop]).sum(1) / ab2[start:stop], 0.0, 1.0)
            d = np.hypot(*(AP - t[:, None] * AB[start:stop]).T)
            hit = np.nonzero(d <= tol)[0]
            if len(hit):
                jj = start + int(hit[0])
                out[k] = jj + float(t[hit[0]])
                j = jj
                found = True
                break
            start = stop
        if not found:
            raise ValueError(f"point {k} is not on the polyline (tolerance {tol:g})")
    return out


def points_at_index(poly: np.ndarray, f: np.ndarray) -> np.ndarray:
    """Points of ``poly`` at fractional vertex indices f (linear along each edge)."""
    f = np.clip(np.asarray(f, dtype=np.float64), 0, len(poly) - 1)
    i = np.minimum(np.floor(f).astype(int), len(poly) - 2)
    t = (f - i)[:, None]
    return poly[i] * (1 - t) + poly[i + 1] * t


def index_to_s(poly: np.ndarray, f: np.ndarray) -> np.ndarray:
    s = arclength(poly)
    f = np.clip(np.asarray(f, dtype=np.float64), 0, len(poly) - 1)
    i = np.minimum(np.floor(f).astype(int), len(poly) - 2)
    return s[i] + (f - i) * (s[i + 1] - s[i])


def eval_samples(poly: np.ndarray, n_uniform: int):
    """Reference evaluation points: all vertices plus n_uniform points at uniform arc length.

    Returns (points, s_positions, uniform_mask).
    """
    s = arclength(poly)
    su = np.linspace(0.0, s[-1], n_uniform)
    pu = np.column_stack([np.interp(su, s, poly[:, 0]), np.interp(su, s, poly[:, 1])])
    pts = np.vstack([poly, pu])
    spos = np.concatenate([s, su])
    mask = np.r_[np.zeros(len(poly), bool), np.ones(n_uniform, bool)]
    return pts, spos, mask
