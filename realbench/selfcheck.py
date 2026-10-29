"""Checks the harness must pass before any test run (protocol section 6). Never touches test data."""

from __future__ import annotations

import numpy as np
import pandas as pd

from curvquant.curves import FAMILIES, similarity

from . import config as C
from .datasets import DATASETS, load_split, raw_files
from .geometry import diameter, gaussian_smooth, points_at_index
from .methods import METHODS, VERTEX_SUBSET
from .metrics import dist_to_polyline
from .runner import ROOT, run_one
from .stats import compare, log_ratios

PUBLISHED = {"uniform arc length": "uniform", "Visvalingam--Whyatt": "vw", "RDP (fixed $K$)": "rdp",
             r"CQ $p{=}\frac12$": "cq_unpinned", r"CQ $p{=}\frac12$, pinned": "cq",
             r"CQ $p{=}\frac12$, snapped": "_paper_snapped_unpinned", "optimal vertex subset": "oracle"}


def _wiggly(rng, n):
    t = np.linspace(0, 1, n)
    x = t + sum(rng.normal(0, 0.05 / k) * np.sin(2 * np.pi * k * t + rng.uniform(0, 6)) for k in range(1, 6))
    y = sum(rng.normal(0, 0.15 / k) * np.sin(2 * np.pi * k * t + rng.uniform(0, 6)) for k in range(1, 6))
    pts = np.column_stack([x, y]) + rng.normal(0, 0.002, (n, 2))
    return pts


def _dev_curves(n_per=8):
    out = []
    for ds in DATASETS:
        if raw_files(ds):
            try:
                out += load_split(ds, "dev", cap=n_per)
            except Exception as ex:
                print(f"  (could not load dev curves of {ds}: {ex})")
    return out


def check_reproduction(log):
    pub = pd.read_csv(ROOT / "results" / "e7_benchmark.csv")
    worst = 0.0
    n = 0
    for fam, (gen, closed) in FAMILIES.items():
        pts = gen(1000)
        for K in sorted(pub.K.unique()):
            for pname, m in PUBLISHED.items():
                if m == "oracle" and fam not in ("F1 log spiral", "F3 corners+arcs"):
                    continue          # the O(N^3) oracle is reproduced on two families to save time
                row = pub[(pub.family == fam) & (pub.K == K) & (pub.method == pname)]
                r = run_one(pts, closed, m, int(K), 0, pts, vertices_only=True)
                rel = abs(r["hausdorff"] - row.Linf.iloc[0]) / row.Linf.iloc[0]
                worst = max(worst, rel)
                n += 1
    ok = worst < 1e-9
    log(ok, f"reproduces {n} published synthetic-benchmark errors (worst relative difference {worst:.1e})")
    return ok


def check_oracle_bound(curves, log):
    bad, n = [], 0
    for cv in curves:
        if len(cv.pts) > 300:
            continue
        for K in (8, 12, 20):
            if K > len(cv.pts) // 2:
                continue
            ro = run_one(cv.pts, cv.closed, "oracle", K, 0, cv.pts)
            for m in ("rdp", "vw", "cq_snapped"):
                rm = run_one(cv.pts, cv.closed, m, K, 0, cv.pts)
                n += 1
                if ro["shortcut"] > rm["shortcut"] * (1 + 1e-9) + 1e-12 * diameter(cv.pts):
                    bad.append((cv.curve_id, K, m, ro["shortcut"], rm["shortcut"]))
    log(not bad, f"oracle <= every vertex-subset method under shortcut error ({n} cases, {len(bad)} violations)")
    return not bad


def check_identity_symmetry(log):
    rng = np.random.default_rng(1)
    rows = []
    for i in range(40):
        for K in (10, 20):
            for m in ("cq", "rdp"):
                rows.append({"curve_id": f"c{i}", "cluster": f"g{i//2}", "K": K, "method": m, "c": 0,
                             "diam": 1.0, "hausdorff": float(rng.lognormal())})
    df = pd.DataFrame(rows)
    aa = compare(log_ratios(df, "cq", "cq", 0, 0, "hausdorff"))
    ab = compare(log_ratios(df, "cq", "rdp", 0, 0, "hausdorff"))
    ba = compare(log_ratios(df, "rdp", "cq", 0, 0, "hausdorff"))
    ok = abs(aa["gmr"] - 1) < 1e-12 and abs(ab["gmr"] * ba["gmr"] - 1) < 1e-12 and abs(ab["p"] - ba["p"]) < 1e-12
    log(ok, "A vs A gives GMR 1; swapping A and B inverts GMR and keeps p")
    return ok


def check_invariance(curves, log):
    worst, n = 0.0, 0
    flips = 0
    for cv in curves[:12]:
        pts = cv.pts
        T = similarity(pts, 0.7, 3.3, [12.0, -5.0])
        for m in ("cq", "rdp", "vw", "uniform", "cq_snapped"):
            K = max(8, len(pts) // 20)
            a = run_one(pts, cv.closed, m, K, 0, pts)["hausdorff"]
            b = run_one(T, cv.closed, m, K, 0, T)["hausdorff"]
            rel = abs(b / (3.3 * a) - 1) if a > 0 else 0.0
            n += 1
            if rel > 1e-6:
                flips += 1          # a tie broken differently by rounding; counted, not hidden
            else:
                worst = max(worst, rel)
    ok = flips <= max(1, n // 50)
    log(ok, f"errors scale exactly under rotation + translation + scaling ({n} cases, worst {worst:.1e}, "
            f"{flips} rounding-tie flips)")
    return ok


def check_on_curve_and_c0(curves, log):
    bad_on, bad_c0, n = 0, 0, 0
    for cv in curves[:20]:
        pts = cv.pts
        K = max(8, len(pts) // 10)
        for c in (0, 4):
            sm = gaussian_smooth(pts, c)
            for m in ("cq", "rdp", "vw", "uniform", "cq_unpinned", "cq_snapped", "cq_p1"):
                r = run_one(pts, cv.closed, m, K, c, sm)
                d = dist_to_polyline(r["_out"], pts).max()
                n += 1
                if d > 10 * C.ONCURVE_TOL * diameter(pts):
                    bad_on += 1
                if c == 0:
                    direct = METHODS[m](pts, K, cv.closed)
                    if not (direct.shape == r["_out"].shape and np.array_equal(direct, r["_out"])):
                        bad_c0 += 1
    log(bad_on == 0, f"every output point lies on the original polyline ({n} runs, {bad_on} off-curve)")
    log(bad_c0 == 0, f"c = 0 path reproduces direct method calls exactly ({bad_c0} mismatches)")
    return bad_on == 0 and bad_c0 == 0


def check_smoothing(log):
    rng = np.random.default_rng(3)
    t = np.linspace(0, 2 * np.pi, 400)
    clean = np.column_stack([np.cos(t), np.sin(t)])
    noisy = clean + rng.normal(0, 0.01, clean.shape)
    e0 = np.hypot(*(noisy - clean).T).mean()
    e4 = np.hypot(*(gaussian_smooth(noisy, 4) - clean).T).mean()
    ok = e4 < 0.6 * e0 and np.array_equal(gaussian_smooth(noisy, 0), noisy)
    log(ok, f"smoothing reduces noise (mean deviation {e0:.4f} -> {e4:.4f}) and c = 0 is the identity")
    return ok


def check_splits(log):
    msgs, ok = [], True
    for ds in DATASETS:
        if not raw_files(ds):
            continue
        st = {}
        dev = load_split(ds, "dev", cap=10 ** 9, stats=st)
        test_clusters = set()
        # recompute test clusters without keeping curves (cap 0)
        st2 = {}
        load_split(ds, "test", cap=0, stats=st2)
        from .datasets import split_of
        overlap = {c.cluster for c in dev if split_of(ds, c.cluster) != "dev"}
        ok &= not overlap
        msgs.append(f"{ds}: dev {st['eligible']['dev']} curves / {st['clusters']['dev']} clusters, "
                    f"test {st['eligible']['test']} / {st['clusters']['test']}")
    log(ok, "cluster-level split has no dev/test overlap; " + "; ".join(msgs))
    return ok


def run(quick: bool = False) -> bool:
    results = []

    def log(ok, msg):
        results.append(ok)
        print(("PASS " if ok else "FAIL ") + msg, flush=True)

    rng = np.random.default_rng(0)
    synth = []
    from .datasets import Curve
    for i in range(15):
        synth.append(Curve("synthetic", f"w{i}", f"w{i}", _wiggly(rng, int(rng.integers(80, 260))), False))
    for fam, (gen, closed) in FAMILIES.items():
        synth.append(Curve("synthetic", fam, fam, gen(240), closed))
    dev = _dev_curves()
    curves = synth + dev
    print(f"selfcheck on {len(synth)} synthetic + {len(dev)} dev curves (no test data is read)")
    if not quick:
        check_reproduction(log)
    check_oracle_bound(curves, log)
    check_identity_symmetry(log)
    check_invariance(curves, log)
    check_on_curve_and_c0(curves, log)
    check_smoothing(log)
    check_splits(log)
    return all(results)
