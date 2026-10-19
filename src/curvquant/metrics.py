"""Approximation-error metrics between a reference polyline and a simplification."""

from __future__ import annotations

import numpy as np


def point_to_polyline_distances(points: np.ndarray, poly: np.ndarray, chunk: int = 2048) -> np.ndarray:
    """Distance from each point to the nearest segment of ``poly``."""
    A, B = poly[:-1], poly[1:]
    AB = B - A
    ab2 = np.maximum((AB ** 2).sum(1), 1e-30)
    out = np.empty(len(points))
    for s in range(0, len(points), chunk):
        P = points[s:s + chunk]
        AP = P[:, None, :] - A[None, :, :]
        t = np.clip((AP * AB[None]).sum(-1) / ab2[None], 0.0, 1.0)
        d = AP - t[..., None] * AB[None]
        out[s:s + chunk] = np.sqrt((d ** 2).sum(-1)).min(1)
    return out


def approximation_errors(reference: np.ndarray, simplified: np.ndarray) -> dict:
    """One-sided errors from every reference vertex to the simplified polyline.

    L_inf is the directed Hausdorff distance (reference -> simplification),
    the standard objective for polyline simplification.
    """
    d = point_to_polyline_distances(reference, simplified)
    return {"Linf": float(d.max()), "RMSE": float(np.sqrt(np.mean(d ** 2)))}
