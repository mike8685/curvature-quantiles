"""Fast harness tests that need no downloaded data (run in CI)."""

import numpy as np
import pandas as pd
import pytest

from curvquant.curves import FAMILIES, similarity
from realbench import config as C
from realbench.datasets import split_of
from realbench.geometry import gaussian_smooth, index_to_s, locate, points_at_index
from realbench.metrics import dist_to_polyline, errors
from realbench.runner import budgets, run_one
from realbench.stats import compare, decide, holm, log_ratios


def _curve(n=200, seed=0):
    rng = np.random.default_rng(seed)
    t = np.linspace(0, 1, n)
    return np.column_stack([t, 0.2 * np.sin(6 * t) + 0.05 * np.sin(23 * t + 1)]) + rng.normal(0, 1e-3, (n, 2))


def test_budgets_floor_and_cap():
    assert budgets(1000) == [(0.02, 20), (0.05, 50), (0.10, 100), (0.20, 200)]
    assert budgets(60) == [(0.02, 8), (0.20, 12)]
    assert all(K <= 0.5 * 30 for _, K in budgets(30))


def test_locate_roundtrip():
    p = _curve()
    f = np.array([0, 3.25, 10.5, 77.0, 150.9, 199])
    q = points_at_index(p, f)
    assert np.allclose(locate(p, q, 1e-12), f)
    assert np.isclose(index_to_s(p, np.array([0.0]))[0], 0.0)


def test_metrics_zero_for_identity_and_positive_otherwise():
    p = _curve()
    s = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(p, axis=0).T))])
    e = errors(p, p, s)
    assert max(e.values()) < 1e-12
    e2 = errors(p, p[[0, -1]], s[[0, -1]])
    assert e2["hausdorff"] > 0 and e2["shortcut"] >= e2["hausdorff"] - 1e-12


def test_distance_matches_bruteforce():
    p = _curve(80)
    q = np.random.default_rng(2).uniform(-0.2, 1.2, (300, 2))
    A, B = p[:-1], p[1:]
    AB = B - A
    t = np.clip(((q[:, None] - A) * AB).sum(-1) / (AB ** 2).sum(-1), 0, 1)
    brute = np.hypot(*(q[:, None] - A - t[..., None] * AB).transpose(2, 0, 1)).min(1)
    assert np.allclose(dist_to_polyline(q, p), brute, atol=1e-12)


@pytest.mark.parametrize("m", ["cq", "rdp", "vw", "uniform", "cq_snapped", "cq_unpinned", "cq_p1", "oracle"])
def test_methods_on_curve_and_c0_exact(m):
    p = _curve(150)
    r = run_one(p, False, m, 15, 0, p)
    assert not r["failed"]
    assert dist_to_polyline(r["_out"], p).max() < 1e-9
    r4 = run_one(p, False, m, 15, 4, gaussian_smooth(p, 4))
    assert dist_to_polyline(r4["_out"], p).max() < 1e-9


def test_oracle_bound():
    p = _curve(150, seed=3)
    o = run_one(p, False, "oracle", 12, 0, p)["shortcut"]
    for m in ("rdp", "vw", "cq_snapped"):
        assert o <= run_one(p, False, m, 12, 0, p)["shortcut"] * (1 + 1e-9)


def test_reproduces_published_errors_on_two_families():
    pub = pd.read_csv("results/e7_benchmark.csv")
    for fam in ("F1 log spiral", "F6 ellipse (closed)"):
        gen, closed = FAMILIES[fam]
        pts = gen(1000)
        for name, m in [("RDP (fixed $K$)", "rdp"), (r"CQ $p{=}\frac12$, pinned", "cq")]:
            row = pub[(pub.family == fam) & (pub.K == 50) & (pub.method == name)]
            r = run_one(pts, closed, m, 50, 0, pts, vertices_only=True)
            assert r["hausdorff"] == pytest.approx(row.Linf.iloc[0], rel=1e-9)


def test_similarity_invariance():
    p = _curve(200, seed=5)
    T = similarity(p, 1.1, 2.5, [3, 4])
    for m in ("cq", "rdp"):
        a = run_one(p, False, m, 20, 0, p)["hausdorff"]
        b = run_one(T, False, m, 20, 0, T)["hausdorff"]
        assert b == pytest.approx(2.5 * a, rel=1e-6)


def test_stats_identity_symmetry_and_decisions():
    rng = np.random.default_rng(0)
    rows = [{"curve_id": f"c{i}", "cluster": f"g{i % 25}", "K": 10, "method": m, "c": 0, "diam": 1.0,
             "hausdorff": rng.lognormal() * (0.5 if m == "cq" else 1.0)} for i in range(100) for m in ("cq", "rdp")]
    df = pd.DataFrame(rows)
    assert compare(log_ratios(df, "cq", "cq", 0, 0, "hausdorff"))["gmr"] == pytest.approx(1.0)
    ab = compare(log_ratios(df, "cq", "rdp", 0, 0, "hausdorff"))
    ba = compare(log_ratios(df, "rdp", "cq", 0, 0, "hausdorff"))
    assert ab["gmr"] * ba["gmr"] == pytest.approx(1.0)
    assert ab["n_clusters"] == 25
    assert decide({"gmr": 0.5, "ci_bonf_lo": 0.4, "ci_bonf_hi": 0.6, "ci_tost_lo": 0.4, "ci_tost_hi": 0.6, "p": 1e-6}) == "CQ better"
    assert decide({"gmr": 1.0, "ci_bonf_lo": 0.9, "ci_bonf_hi": 1.1, "ci_tost_lo": 0.97, "ci_tost_hi": 1.03, "p": 0.5}) == "practically equivalent"
    assert decide({"gmr": 1.2, "ci_bonf_lo": 0.9, "ci_bonf_hi": 1.5, "ci_tost_lo": 0.95, "ci_tost_hi": 1.4, "p": 0.2}) == "inconclusive"
    h = holm({"a": 0.001, "b": 0.04}, 18)
    assert h["a"] == pytest.approx(0.018) and h["b"] == pytest.approx(0.68)


def test_split_is_deterministic_and_cluster_based():
    assert split_of("gpx", "file1.gpx") == split_of("gpx", "file1.gpx")
    frac = np.mean([split_of("x", f"c{i}") == "dev" for i in range(20000)])
    assert abs(frac - C.DEV_FRACTION) < 0.02


def test_exact_fits_tie_instead_of_rounding_ratios():
    df = pd.DataFrame([
        {"curve_id": "a", "cluster": "a", "K": 8, "method": "cq", "c": 0, "diam": 1.0, "hausdorff": 2e-16},
        {"curve_id": "a", "cluster": "a", "K": 8, "method": "rdp", "c": 0, "diam": 1.0, "hausdorff": 1e-16},
        {"curve_id": "b", "cluster": "b", "K": 8, "method": "cq", "c": 0, "diam": 1.0, "hausdorff": 0.0},
        {"curve_id": "b", "cluster": "b", "K": 8, "method": "rdp", "c": 0, "diam": 1.0, "hausdorff": 3e-16},
    ])
    lr = log_ratios(df, "cq", "rdp", 0, 0, "hausdorff")
    assert (lr.lr == 0).all() and lr.exact_a.all() and lr.exact_b.all()


def test_smoothing_has_no_boundary_curvature():
    from curvquant.core import turning_angles
    t = np.linspace(0, 2 * np.pi, 401)
    circ = np.column_stack([np.cos(t), np.sin(t)])
    _, w = turning_angles(gaussian_smooth(circ, 16, closed=True), closed=True)
    assert w.max() / np.median(w) < 1.1
    a = np.linspace(0, 2, 300)
    _, w = turning_angles(gaussian_smooth(np.column_stack([np.cos(a), np.sin(a)]), 16))
    assert 0.9 < w.min() / np.median(w) and w.max() / np.median(w) < 1.1


def test_uji_ids_unique_and_writer_clusters(tmp_path, monkeypatch):
    (tmp_path / "uji").mkdir()
    lines = []
    for sess in ("trn_UJI_W03-01", "trn_UJI_W04-01"):
        for ch in ("a", "a", "b"):
            lines += [f"WORD {ch} {sess}", "  NUMSTROKES 2",
                      "  POINTS 3 # 0 0 1 1 2 0", "  POINTS 3 # 5 5 6 6 7 5"]
    (tmp_path / "uji" / "ujipenchars2.txt").write_text("\n".join(lines))
    monkeypatch.setenv("REALBENCH_DATA", str(tmp_path))
    from realbench.datasets import DATASETS
    cs = list(DATASETS["uji"].loader())
    ids = [c.curve_id for c in cs]
    assert len(ids) == 12 and len(set(ids)) == 12
    assert {c.cluster for c in cs} == {"W03", "W04"}


def test_omniglot_strokes_split_on_break_and_cluster_by_character(tmp_path, monkeypatch):
    d = tmp_path / "omniglot" / "strokes_background" / "Latin" / "character01"
    d.mkdir(parents=True)
    (d / "0001_01.txt").write_text("START\n0,0,0\n1,0,5\n2,1,9\nBREAK\n5,5,20\n6,6,25\n")
    (d / "0001_02.txt").write_text("START\n0,0,0\n1,1,5\nbad,line\n")
    monkeypatch.setenv("REALBENCH_DATA", str(tmp_path))
    from realbench.datasets import DATASETS
    cs = list(DATASETS["omniglot"].loader())
    assert [len(c.pts) for c in cs] == [3, 2, 2]
    assert {c.cluster for c in cs} == {"Latin/character01"}
    assert len({c.curve_id for c in cs}) == 3


def test_test_split_refuses_without_freeze(monkeypatch):
    from realbench import lock
    monkeypatch.setattr(lock, "REGISTRY", lock.HERE / "no-such-registry.json")
    from realbench.runner import run_test
    with pytest.raises(SystemExit):
        run_test("ne_roads", allow_unpushed=True)
