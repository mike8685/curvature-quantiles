"""Deterministic test curves used in the paper."""

from __future__ import annotations

import numpy as np
from scipy.special import fresnel


def normalize(pts: np.ndarray) -> np.ndarray:
    """Center the bounding box at the origin and scale its larger side to 1."""
    lo, hi = pts.min(0), pts.max(0)
    return (pts - 0.5 * (lo + hi)) / max(float((hi - lo).max()), 1e-12)


def graph(f, a: float, b: float, n: int) -> np.ndarray:
    x = np.linspace(a, b, n)
    return np.column_stack([x, f(x)])


def log_spiral(n: int) -> np.ndarray:
    th = np.linspace(0.0, 4 * np.pi, n)
    r = 0.1 * np.exp(0.22 * th)
    return np.column_stack([r * np.cos(th), r * np.sin(th)])


def clothoid(n: int) -> np.ndarray:
    s = np.linspace(-2.2, 2.2, n)
    S, C = fresnel(s)
    return normalize(np.column_stack([C, S]))


def corners_and_arcs(n: int) -> np.ndarray:
    n1 = n // 3
    n2 = n // 3
    n3 = n - n1 - n2
    t1 = np.linspace(np.pi, 0.0, n1, endpoint=False)
    arc = np.column_stack([np.cos(t1), np.sin(t1)])
    t2 = np.linspace(0.0, 1.0, n2, endpoint=False)
    fillet = np.column_stack([1.0 + 1.5 * t2, -0.5 * (1.0 - np.cos(np.pi * t2))])
    ax = np.array([2.5, 3.0, 3.5, 4.0, 4.5, 5.2])
    ay = np.array([-1.0, 0.6, -0.8, 0.7, -0.6, 0.0])
    u = np.linspace(0.0, 1.0, n3)
    ua = np.linspace(0.0, 1.0, len(ax))
    zig = np.column_stack([np.interp(u, ua, ax), np.interp(u, ua, ay)])
    return normalize(np.vstack([arc, fillet, zig]))


def bump_line(n: int) -> np.ndarray:
    x = np.linspace(0.0, 10.0, n)
    y = 0.002 * x ** 2 + 0.35 * np.exp(-0.5 * ((x - 3.0) / 0.12) ** 2) - 0.25 * np.exp(-0.5 * ((x - 7.5) / 0.08) ** 2)
    return normalize(np.column_stack([x, y]))


def triple_sigmoid(n: int) -> np.ndarray:
    x = np.linspace(-12.0, 12.0, n)
    y = 1 / (1 + np.exp(-2.5 * (x + 6))) - 1 / (1 + np.exp(-3.0 * x)) + 1 / (1 + np.exp(-2.5 * (x - 6)))
    return normalize(np.column_stack([x, y]))


def ellipse(n: int) -> np.ndarray:
    t = np.linspace(0.0, 2 * np.pi, n)
    pts = np.column_stack([2.5 * np.cos(t), 0.4 * np.sin(t)])
    pts[-1] = pts[0]
    return normalize(pts)


def fourier_star(n: int) -> np.ndarray:
    t = np.linspace(0.0, 2 * np.pi, n)
    r = 1.0 + 0.35 * np.cos(5 * t) + 0.15 * np.sin(3 * t)
    pts = np.column_stack([r * np.cos(t), r * np.sin(t)])
    pts[-1] = pts[0]
    return normalize(pts)


def lissajous(n: int) -> np.ndarray:
    t = np.linspace(0.0, 2 * np.pi, n)
    pts = np.column_stack([np.sin(2 * t) + 0.5 * np.cos(t), np.cos(2 * t) * np.sin(t)])
    pts[-1] = pts[0]
    return normalize(pts)


def mixed_trajectory(n: int) -> np.ndarray:
    t = np.linspace(0.0, 1.0, n)
    g = np.exp(-(((t - 0.35) / 0.12) ** 2))
    x = 8.0 * t + 0.45 * np.sin(18 * np.pi * t) * g
    y = 0.8 * np.tanh(15 * (t - 0.72)) + 0.45 * np.cos(18 * np.pi * t) * g + 0.15 * np.sin(4 * np.pi * t)
    return normalize(np.column_stack([x, y]))


FAMILIES = {
    "F1 log spiral": (log_spiral, False),
    "F2 clothoid": (clothoid, False),
    "F3 corners+arcs": (corners_and_arcs, False),
    "F4 line+bumps": (bump_line, False),
    "F5 triple sigmoid": (triple_sigmoid, False),
    "F6 ellipse (closed)": (ellipse, True),
    "F7 Fourier star (closed)": (fourier_star, True),
    "F8 Lissajous (closed)": (lissajous, True),
    "F9 mixed trajectory": (mixed_trajectory, False),
}


def similarity(pts: np.ndarray, angle: float, scale: float, shift) -> np.ndarray:
    c, s = np.cos(angle), np.sin(angle)
    R = np.array([[c, -s], [s, c]])
    return scale * pts @ R.T + np.asarray(shift)
