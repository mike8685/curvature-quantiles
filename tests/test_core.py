import numpy as np
import pytest

from curvquant import (
    approximation_errors, curvature_median, optimal_subset_min_eps, profile,
    profile_distance, profile_distance_quantile_form, quantile_resample,
    rdp_fixed_k, uniform_arclength, vw_fixed_k,
)
from curvquant.curves import graph, log_spiral, similarity


def test_total_turning_equals_tangent_angle_variation():
    # convex graph: total |d theta| = arctan f'(b) - arctan f'(a)
    pts = graph(lambda x: x ** 2, 0, 2, 2001)
    e = np.diff(pts, axis=0)
    chord_turn = np.arctan2(e[-1, 1], e[-1, 0]) - np.arctan2(e[0, 1], e[0, 0])
    assert profile(pts).K == pytest.approx(chord_turn, abs=1e-12)        # exact, discrete
    assert profile(pts).K == pytest.approx(np.arctan(4.0), abs=2e-3)     # O(h) to continuum


@pytest.mark.parametrize("f,fp,inv,a,b", [
    (lambda x: x ** 2, lambda x: 2 * x, lambda m: m / 2, 0.0, 2.0),
    (np.sqrt, lambda x: 0.5 / np.sqrt(x), lambda m: (0.5 / m) ** 2, 1.0, 9.0),
    (lambda x: 1 - np.exp(-x), lambda x: np.exp(-x), lambda m: -np.log(m), 0.0, 4.0),
])
def test_theorem1_halfway_angle(f, fp, inv, a, b):
    exact = inv(np.tan(0.5 * (np.arctan(fp(a)) + np.arctan(fp(b)))))
    assert curvature_median(graph(f, a, b, 4001))[0] == pytest.approx(exact, abs=2e-3)


def test_theorem1_closed_form_parabola():
    assert curvature_median(graph(lambda x: x ** 2, 0, 2, 20001))[0] == pytest.approx((np.sqrt(17) - 1) / 8, abs=1e-4)


def test_symmetric_sigmoid_median_is_inflection():
    for n in (20, 21, 200, 201):
        pts = graph(lambda x: 1 / (1 + np.exp(-x)), -6, 6, n)
        assert abs(curvature_median(pts)[0]) < 1e-9


def test_similarity_and_parameterization_invariance():
    x = np.linspace(0, 2, 801)
    a = np.column_stack([x, 1 - np.exp(-4 * x)])
    b = similarity(a, 1.1, 3.7, (2, -5))
    assert profile_distance(profile(a), profile(b)) < 1e-10
    assert profile_distance(profile(a, beta=1.0), profile(b, beta=1.0)) < 1e-10
    xc = np.linspace(0, 2, 211)
    c = np.column_stack([xc, 1 - np.exp(-4 * xc)])
    assert profile_distance(profile(a), profile(c)) < 5e-3


def test_w1_identity_cdf_vs_quantile_form():
    x = np.linspace(0, 2, 601)
    a = profile(np.column_stack([x, 1 - np.exp(-4 * x)]))
    b = profile(np.column_stack([x, np.exp(4 * (x - 2))]))
    assert profile_distance(a, b) == pytest.approx(profile_distance_quantile_form(a, b), abs=1e-5)


def test_convex_arc_reduces_to_normalized_turning_function():
    # Proposition: for a convex arc F(u) = (theta(u) - theta(0)) / K, so D equals the
    # L1 distance of normalized turning functions with fixed start point and rotation
    x = np.linspace(0, 2, 3001)
    a = np.column_stack([x, x ** 2])
    b = np.column_stack([x, np.exp(x)])
    u = np.linspace(0, 1, 3001)

    def turning(pts):
        e = np.diff(pts, axis=0)
        th = np.unwrap(np.arctan2(e[:, 1], e[:, 0]))
        L = np.hypot(*e.T)
        s = (np.cumsum(L) - 0.5 * L) / L.sum()
        T = (np.interp(u, s, th) - th[0]) / (th[-1] - th[0])
        return T

    d_turn = np.trapezoid(np.abs(turning(a) - turning(b)), u)
    assert profile_distance(profile(a), profile(b)) == pytest.approx(d_turn, abs=2e-3)


def test_hybrid_degrades_to_arclength_on_straight_line():
    pts = np.column_stack([np.linspace(0, 1, 101), np.zeros(101)])
    with pytest.raises(ValueError):
        profile(pts)
    pr = profile(pts, beta=1.0)
    u = np.linspace(0, 1, 11)
    assert np.allclose(pr.cdf(u), u, atol=1e-12)


def test_resampling_budget_and_endpoints():
    pts = log_spiral(1500)
    for K in (10, 30, 77):
        for kw in ({}, {"beta": None, "p": 0.5}, {"beta": None, "p": 0.5, "pin_corners": True}):
            out = quantile_resample(pts, K, **kw)
            assert len(out) == K
            assert np.allclose(out[0], pts[0]) and np.allclose(out[-1], pts[-1])
        assert len(rdp_fixed_k(pts, K)) == K
        assert len(vw_fixed_k(pts, K)) == K
        assert len(uniform_arclength(pts, K)) == K


def test_optimal_subset_is_lower_bound_for_decimation():
    pts = log_spiral(300)
    K = 15
    simp, eps = optimal_subset_min_eps(pts, K)
    assert len(simp) == K
    assert approximation_errors(pts, simp)["Linf"] == pytest.approx(eps, rel=1e-9)
    for other in (rdp_fixed_k(pts, K), vw_fixed_k(pts, K)):
        assert approximation_errors(pts, other)["Linf"] >= eps - 1e-12
