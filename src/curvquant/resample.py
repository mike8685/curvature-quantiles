"""Polyline simplification: curvature-quantile resampling and classical baselines."""

from __future__ import annotations

import heapq

import numpy as np

from .core import _as_points, profile, turning_angles


def _arclength(pts):
    lengths = np.hypot(*np.diff(pts, axis=0).T)
    return np.concatenate([[0.0], np.cumsum(lengths)])


def _points_at_s(pts, s_vert, s_target):
    return np.column_stack([np.interp(s_target, s_vert, pts[:, 0]), np.interp(s_target, s_vert, pts[:, 1])])


def uniform_arclength(pts, K: int, closed: bool = False) -> np.ndarray:
    pts = _as_points(pts)
    s = _arclength(pts)
    out = _points_at_s(pts, s, np.linspace(0.0, s[-1], K))
    if closed:
        out[-1] = out[0]
    return out


def quantile_resample(
    pts,
    K: int,
    beta: float | None = 1.0,
    p: float = 1.0,
    closed: bool = False,
    pin_corners: bool = False,
    corner_threshold: float = np.pi / 6,
    snap_to_vertices: bool = False,
) -> np.ndarray:
    """Place K points (endpoints included) at equal quantiles of the measure.

    pin_corners: vertices whose turning angle exceeds ``corner_threshold`` are
        kept exactly; the remaining budget is shared among the pieces between
        pinned vertices in proportion to their mass (largest remainder).
    snap_to_vertices: move every output point to the nearest input vertex in
        arc length (turns resampling into decimation; duplicates are dropped,
        so the result can have slightly fewer than K points).
    """
    pts = _as_points(pts)
    if K >= len(pts):
        return pts.copy()
    prof = profile(pts, beta=beta, p=p, closed=closed)
    s = _arclength(pts)
    L = s[-1]
    u_vert = s / L

    q_targets = None
    if pin_corners and K > 4:
        _, w = turning_angles(pts, closed=False)
        corners = np.where(w > corner_threshold)[0] + 1  # vertex indices
        if 0 < len(corners) <= K - 2:
            anchors_u = np.concatenate([[0.0], u_vert[corners], [1.0]])
            anchors_q = prof.cdf(anchors_u)
            dq = np.diff(anchors_q)
            budget = K - len(anchors_u)
            raw = budget * dq / max(dq.sum(), 1e-15)
            alloc = np.floor(raw).astype(int)
            short = budget - int(alloc.sum())
            if short > 0:
                alloc[np.argsort(-(raw - alloc))[:short]] += 1
            u_list = []
            for j in range(len(dq)):
                u_list.append([anchors_u[j]])
                if alloc[j] > 0:
                    qs = np.linspace(anchors_q[j], anchors_q[j + 1], alloc[j] + 2)[1:-1]
                    u_list.append(prof.quantile(qs))
            u_list.append([1.0])
            u_targets = np.concatenate(u_list)
            q_targets = "done"
    if q_targets is None:
        if closed:
            q = np.linspace(0.0, 1.0, K)  # last point closes the loop
        else:
            q = np.linspace(0.0, 1.0, K)
        u_targets = prof.quantile(q)
        u_targets[0], u_targets[-1] = 0.0, 1.0

    s_targets = np.sort(u_targets) * L
    if snap_to_vertices:
        idx = np.clip(np.searchsorted(s, s_targets), 1, len(s) - 1)
        left_closer = (s_targets - s[idx - 1]) < (s[idx] - s_targets)
        idx = np.where(left_closer, idx - 1, idx)
        idx = np.unique(np.concatenate([[0], idx, [len(s) - 1]]))
        out = pts[idx].copy()
    else:
        out = _points_at_s(pts, s, s_targets)
    if closed:
        out[-1] = out[0]
    return out


def _perp_dists(pts, i, j):
    a, b = pts[i], pts[j]
    sub = pts[i + 1:j]
    v = b - a
    n = np.hypot(v[0], v[1])
    if n < 1e-15:
        return np.hypot(*(sub - a).T)
    return np.abs(v[0] * (a[1] - sub[:, 1]) - v[1] * (a[0] - sub[:, 0])) / n


def rdp_fixed_k(pts, K: int) -> np.ndarray:
    """Ramer-Douglas-Peucker driven to exactly K vertices (largest deviation first)."""
    pts = _as_points(pts)
    N = len(pts)
    if K >= N:
        return pts.copy()

    def best_split(i, j):
        if j <= i + 1:
            return -1.0, -1
        d = _perp_dists(pts, i, j)
        k = int(np.argmax(d))
        return float(d[k]), i + 1 + k

    kept = {0, N - 1}
    d0, k0 = best_split(0, N - 1)
    heap = [(-d0, k0, 0, N - 1)]
    while len(kept) < K and heap:
        nd, k, i, j = heapq.heappop(heap)
        if k < 0:
            break
        kept.add(k)
        for a, b in ((i, k), (k, j)):
            d, kk = best_split(a, b)
            if kk >= 0:
                heapq.heappush(heap, (-d, kk, a, b))
    return pts[sorted(kept)]


def vw_fixed_k(pts, K: int) -> np.ndarray:
    """Visvalingam-Whyatt (effective area, min-heap) down to K vertices."""
    pts = _as_points(pts)
    N = len(pts)
    if K >= N:
        return pts.copy()
    prev = np.arange(N) - 1
    nxt = np.arange(N) + 1
    removed = np.zeros(N, dtype=bool)
    area = np.full(N, np.inf)

    def tri(i):
        a, c = prev[i], nxt[i]
        if a < 0 or c >= N:
            return np.inf
        p, q, r = pts[a], pts[i], pts[c]
        return 0.5 * abs(p[0] * (q[1] - r[1]) + q[0] * (r[1] - p[1]) + r[0] * (p[1] - q[1]))

    for i in range(1, N - 1):
        area[i] = tri(i)
    heap = [(area[i], i) for i in range(1, N - 1)]
    heapq.heapify(heap)
    remaining = N
    last = 0.0
    while remaining > K and heap:
        a, i = heapq.heappop(heap)
        if removed[i] or a != area[i]:
            continue
        last = max(last, a)
        removed[i] = True
        remaining -= 1
        l, r = prev[i], nxt[i]
        nxt[l] = r
        prev[r] = l
        for nb in (l, r):
            if 0 < nb < N - 1 and not removed[nb]:
                area[nb] = max(tri(nb), last)
                heapq.heappush(heap, (area[nb], nb))
    return pts[~removed]


def optimal_subset_min_eps(pts, K: int) -> tuple[np.ndarray, float]:
    """Exact min-max vertex-subset approximation with K vertices (endpoints kept).

    Error of a shortcut (i, j) is the largest distance from an intermediate
    vertex to the *segment* P_i P_j.  O(N^2) memory, O(N^3 + K N^2) time:
    intended as a gold-standard reference for N <= ~1000.
    """
    pts = _as_points(pts)
    N = len(pts)
    err = np.zeros((N, N))
    for i in range(N - 2):
        a = pts[i]
        B = pts[i + 2:]                         # candidate endpoints j = i+2..N-1
        M = pts[i + 1:N - 1]                    # intermediate vertices m = i+1..N-2
        AB = B - a
        ab2 = np.maximum((AB ** 2).sum(1), 1e-30)
        AM = M - a
        t = np.clip((AM @ AB.T) / ab2[None, :], 0.0, 1.0)        # (m, j)
        dx = AM[:, 0:1] - t * AB[None, :, 0]
        dy = AM[:, 1:2] - t * AB[None, :, 1]
        D = np.hypot(dx, dy)
        # only m < j counts: row r (m = i+1+r) is valid for column c (j = i+2+c) iff r <= c
        mask = np.triu(np.ones(D.shape, dtype=bool))
        err[i, i + 2:] = np.where(mask, D, 0.0).max(axis=0)
    INF = np.inf
    E = np.full(N, INF)
    E[0] = 0.0
    parents = []
    for _ in range(K - 1):
        cand = np.maximum(E[:, None], err)          # (i, j)
        cand[np.tril_indices(N)] = INF
        par = np.argmin(cand, axis=0)
        E = cand[par, np.arange(N)]
        E[0] = INF
        parents.append(par)
    eps = float(E[N - 1])
    idx = [N - 1]
    for par in reversed(parents):
        idx.append(int(par[idx[-1]]))
    idx = sorted(set(idx))
    return pts[idx], eps
