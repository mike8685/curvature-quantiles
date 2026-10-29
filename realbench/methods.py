"""Method registry. Wraps the frozen curvquant implementations without changing them.

Each method maps (decision polyline, K, closed) to output points on that polyline.
The runner converts them to fractional vertex indices so they can be placed on the
original polyline when the decision polyline is a smoothed copy (protocol 2.1).
"""

from __future__ import annotations

import numpy as np

from curvquant import optimal_subset_min_eps, quantile_resample, rdp_fixed_k, uniform_arclength, vw_fixed_k

from . import config as C

VERTEX_SUBSET = {"rdp", "vw", "cq_snapped", "oracle", "_paper_snapped_unpinned"}


def _cq(pts, K, closed, **kw):
    return quantile_resample(pts, K, beta=None, p=C.CQ_P, closed=closed,
                             corner_threshold=C.CQ_CORNER_THRESHOLD, **kw)


METHODS = {
    "cq": lambda p, K, closed: _cq(p, K, closed, pin_corners=True),
    "rdp": lambda p, K, closed: rdp_fixed_k(p, K),
    "vw": lambda p, K, closed: vw_fixed_k(p, K),
    "uniform": lambda p, K, closed: uniform_arclength(p, K, closed),
    "cq_unpinned": lambda p, K, closed: _cq(p, K, closed),
    "cq_snapped": lambda p, K, closed: _cq(p, K, closed, pin_corners=True, snap_to_vertices=True),
    "cq_p1": lambda p, K, closed: quantile_resample(p, K, beta=None, p=1.0, closed=closed),
    "oracle": lambda p, K, closed: optimal_subset_min_eps(p, K)[0],
    # the preprint's "snapped" row snapped the UNPINNED variant; kept only so selfcheck can
    # reproduce that published row (not part of any run plan)
    "_paper_snapped_unpinned": lambda p, K, closed: _cq(p, K, closed, snap_to_vertices=True),
}

ABLATIONS = ("uniform", "cq_unpinned", "cq_snapped", "cq_p1")


def vertex_indices(poly: np.ndarray, out: np.ndarray) -> np.ndarray:
    """Indices of output vertices in ``poly``: for each, the first exact match after the previous."""
    idx = np.empty(len(out))
    j = 0
    for k, q in enumerate(out):
        eq = np.nonzero((poly[j:, 0] == q[0]) & (poly[j:, 1] == q[1]))[0]
        if not len(eq):
            raise ValueError("vertex-subset output contains a point that is not an input vertex")
        j = j + int(eq[0])
        idx[k] = j
    return idx
