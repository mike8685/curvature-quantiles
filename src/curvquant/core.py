"""Curvature-mass profiles, quantiles and profile distances for planar polylines.

Discretization
--------------
A polyline P_0, ..., P_n has edges e_i = P_i - P_{i-1} (i = 1..n) with lengths
l_i and total length L.  The unsigned exterior turning angle at interior vertex
P_i (i = 1..n-1) is

    w_i = atan2(|t_i x t_{i+1}|, t_i . t_{i+1}) in [0, pi],

with t_i the unit tangent of edge i.  This is the exact discrete analogue of
|kappa| ds = |d theta|.  The mass w_i is spread uniformly over the dual cell of
P_i, i.e. between the midpoints of edges i and i+1, so the cumulative turning
Phi(u) is piecewise linear in normalized arc length u = s / L with knots at
u = 0, the edge midpoints, and u = 1.

Measures
--------
pure       d nu = |kappa| ds / K,                    K = total absolute curvature
power p    d nu_p  ~ |kappa|^p L^(p-1) ds             (dimensionless, scale invariant)
hybrid     d mu_{p,beta} = ds / L + beta |kappa|^p L^(p-1) ds   (then normalized)

For p = 1 the hybrid weight beta has units of 1/radian: a curve with total
turning K puts the fraction beta K / (1 + beta K) of its mass on turning and the
rest on arc length.  As K -> 0 the hybrid measure degrades continuously to
uniform arc length, which is what keeps it stable on nearly straight curves.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

_EPS = 1e-15


def _as_points(x, y=None) -> np.ndarray:
    pts = np.asarray(x, dtype=np.float64) if y is None else np.column_stack([x, y]).astype(np.float64)
    if pts.ndim != 2 or pts.shape[1] != 2 or len(pts) < 3:
        raise ValueError("need an (n, 2) array of at least 3 points")
    return pts


def turning_angles(pts: np.ndarray, closed: bool = False):
    """Edge lengths and unsigned exterior turning angles of a polyline.

    Returns (lengths, w) where ``lengths`` has one entry per edge and ``w`` one
    entry per interior vertex.  For ``closed=True`` the first and last points are
    assumed equal and the seam angle is appended as the last element of ``w``.
    """
    pts = _as_points(pts)
    edges = np.diff(pts, axis=0)
    lengths = np.hypot(edges[:, 0], edges[:, 1])
    t = edges / np.maximum(lengths[:, None], _EPS)
    t_next = t[1:]
    t_prev = t[:-1]
    if closed:
        t_prev = np.vstack([t_prev, t[-1:]])
        t_next = np.vstack([t_next, t[:1]])
    cross = t_prev[:, 0] * t_next[:, 1] - t_prev[:, 1] * t_next[:, 0]
    dot = np.einsum("ij,ij->i", t_prev, t_next)
    w = np.arctan2(np.abs(cross), dot)
    return lengths, w


def signed_turning_angles(pts: np.ndarray) -> np.ndarray:
    """Signed exterior turning angles at interior vertices (open polyline)."""
    pts = _as_points(pts)
    edges = np.diff(pts, axis=0)
    cross = edges[:-1, 0] * edges[1:, 1] - edges[:-1, 1] * edges[1:, 0]
    dot = np.einsum("ij,ij->i", edges[:-1], edges[1:])
    return np.arctan2(cross, dot)


@dataclass
class Profile:
    """Piecewise-linear cumulative mass on normalized arc length."""

    u: np.ndarray          # knots in [0, 1], nondecreasing
    H: np.ndarray          # normalized cumulative mass at the knots, 0 -> 1
    L: float               # polyline length
    K: float               # total (unsigned) turning in radians, sum of w_i
    total_mass: float      # unnormalized mass of the measure (dimensionless)

    def cdf(self, u) -> np.ndarray:
        return np.interp(u, self.u, self.H)

    def _inverse(self, q: np.ndarray, side: str, tol: float = 1e-12) -> np.ndarray:
        # the tolerance makes plateaus that differ from q only by rounding error
        # count as plateaus at level q
        shifted = q - tol if side == "left" else q + tol
        k = np.clip(np.searchsorted(self.H, shifted, side=side), 1, len(self.H) - 1)
        h0, h1 = self.H[k - 1], self.H[k]
        span = h1 - h0
        t = np.clip(np.where(span > 0, (q - h0) / np.where(span > 0, span, 1.0), 0.0), 0.0, 1.0)
        return self.u[k - 1] + t * (self.u[k] - self.u[k - 1])

    def quantile(self, q, convention: str = "mid") -> np.ndarray:
        """Quantile function P(q).

        ``"left"``  : inf{u : H(u) >= q} (generalized inverse),
        ``"right"`` : sup{u : H(u) <= q},
        ``"mid"``   : their average (default).  The two differ only where H is
        flat at level q, e.g. at a discrete inflection vertex whose turning angle
        is exactly zero; the mid convention keeps symmetric curves symmetric.
        """
        q = np.clip(np.atleast_1d(np.asarray(q, dtype=np.float64)), 0.0, 1.0)
        left = self._inverse(q, "left")
        if convention == "left":
            out = left
        else:
            right = self._inverse(q, "right")
            out = right if convention == "right" else 0.5 * (left + right)
        out = np.where(q <= 0, 0.0, out)
        out = np.where(q >= 1, 1.0, out)
        return out


def profile(
    pts,
    beta: float | None = None,
    p: float = 1.0,
    closed: bool = False,
) -> Profile:
    """Curvature-mass profile of a polyline.

    beta=None (or 0) gives the pure turning measure; beta > 0 the hybrid
    measure d mu = ds/L + beta |kappa|^p L^(p-1) ds.
    """
    pts = _as_points(pts)
    lengths, w = turning_angles(pts, closed=closed)
    L = float(lengths.sum())
    if L <= _EPS:
        raise ValueError("degenerate polyline (zero length)")
    K = float(w.sum())

    seam = 0.0
    if closed:
        w, seam = w[:-1], float(w[-1])

    # dimensionless per-vertex mass for exponent p
    if p == 1.0:
        m = w.copy()
        m_seam = seam
    else:
        dual = 0.5 * (lengths[:-1] + lengths[1:]) / L
        m = np.power(w, p) * np.power(np.maximum(dual, _EPS), 1.0 - p)
        if closed:
            d_seam = 0.5 * (lengths[-1] + lengths[0]) / L
            m_seam = seam ** p * max(d_seam, _EPS) ** (1.0 - p)
        else:
            m_seam = 0.0
    M = float(m.sum() + m_seam)

    s_mid = np.cumsum(lengths) - 0.5 * lengths
    u_mid = s_mid / L
    # cumulative mass at edge midpoints: half the seam mass sits before the first
    # midpoint (and half after the last) on a closed curve
    phi_mid = np.concatenate([[0.0], np.cumsum(m)]) + 0.5 * m_seam
    u = np.concatenate([[0.0], u_mid, [1.0]])
    phi = np.concatenate([[0.0], phi_mid, [M]])

    use_hybrid = beta is not None and beta > 0
    if use_hybrid:
        H = (u + beta * phi) / (1.0 + beta * M)
        total = 1.0 + beta * M
    else:
        if M <= 1e-12:
            raise ValueError(
                "curve has (numerically) zero turning; the pure profile is "
                "undefined. Use beta > 0 (hybrid measure)."
            )
        H = phi / M
        total = M
    H[0], H[-1] = 0.0, 1.0
    H = np.maximum.accumulate(H)
    return Profile(u=u, H=H, L=L, K=K, total_mass=total)


def point_at(pts, u) -> np.ndarray:
    """Point(s) on the polyline at normalized arc length u."""
    pts = _as_points(pts)
    lengths = np.hypot(*np.diff(pts, axis=0).T)
    s = np.concatenate([[0.0], np.cumsum(lengths)])
    target = np.atleast_1d(u) * s[-1]
    return np.column_stack([np.interp(target, s, pts[:, 0]), np.interp(target, s, pts[:, 1])])


def curvature_quantile_points(pts, q, beta=None, p=1.0, closed=False) -> np.ndarray:
    """Points on the curve at curvature-mass quantiles q."""
    prof = profile(pts, beta=beta, p=p, closed=closed)
    return point_at(pts, prof.quantile(q))


def curvature_median(pts, beta=None, p=1.0) -> np.ndarray:
    """Point that splits the (pure or hybrid) curvature mass in half."""
    return curvature_quantile_points(pts, [0.5], beta=beta, p=p)[0]


def quantile_gap(prof: Profile, center: float = 0.5, half_width: float = 0.05) -> float:
    """P(center + h) - P(center - h): large when the median sits in a flat stretch."""
    a, b = prof.quantile([center - half_width, center + half_width])
    return float(b - a)


def profile_distance(prof1: Profile, prof2: Profile, n_grid: int = 4001) -> float:
    """D = integral_0^1 |P1(q) - P2(q)| dq = W1 between the two normalized measures.

    Computed as integral_0^1 |H1(u) - H2(u)| du on the union of both knot sets
    plus a uniform grid (both representations give the same number).
    """
    grid = np.unique(np.concatenate([prof1.u, prof2.u, np.linspace(0.0, 1.0, n_grid)]))
    diff = np.abs(prof1.cdf(grid) - prof2.cdf(grid))
    return float(np.trapezoid(diff, grid))


def profile_distance_quantile_form(prof1: Profile, prof2: Profile, n_grid: int = 20001) -> float:
    """Same distance computed from the quantile functions (used to check W1 identity)."""
    q = np.linspace(0.0, 1.0, n_grid)
    return float(np.trapezoid(np.abs(prof1.quantile(q) - prof2.quantile(q)), q))
