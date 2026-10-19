"""Experiments E6 (stroke classification) and E7 (fixed-budget polyline simplification)."""

from __future__ import annotations

import time

import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
from scipy.stats import wilcoxon

import matplotlib.pyplot as plt

from curvquant import (
    approximation_errors, optimal_subset_min_eps, profile, profile_distance,
    quantile_resample, rdp_fixed_k, uniform_arclength, vw_fixed_k,
)
from curvquant.core import signed_turning_angles
from curvquant.curves import FAMILIES, log_spiral, similarity

from common import C, RESULTS, macro, save_fig, sci, write_table

# =========================================================================== E6


def _stroke(cls: str, rng, jitter: float):
    j = lambda v: v * (1 + jitter * rng.uniform(-1, 1))
    t = np.linspace(0, 1, 2000)
    if cls == "alpha":                       # nodal cubic: loop with two tails
        a = j(1.6)
        tt = -a + 2 * a * t
        x, y = tt ** 2 - 1, tt * (tt ** 2 - 1)
    elif cls == "S":                         # one period of a sine wave
        amp = j(0.6)
        tt = -np.pi + 2 * np.pi * t
        x, y = tt / np.pi, amp * np.sin(tt)
    elif cls == "C":                         # circular arc
        span = j(1.4 * np.pi)
        tt = span * t
        x, y = np.cos(tt), np.sin(tt)
    elif cls == "spiral":                    # Archimedean spiral (decreasing curvature)
        t0, t1 = j(0.5 * np.pi), j(3.0 * np.pi)
        tt = t0 + (t1 - t0) * t
        x, y = tt * np.cos(tt), tt * np.sin(tt)
    elif cls == "J":                         # straight shaft then a semicircle
        shaft, r = j(1.5), j(0.35)
        n1 = int(2000 * shaft / (shaft + np.pi * r))
        xs, ys = np.zeros(n1), np.linspace(shaft, 0, n1, endpoint=False)
        tt = np.linspace(np.pi, 2 * np.pi, 2000 - n1)
        x = np.concatenate([xs, r + r * np.cos(tt)])
        y = np.concatenate([ys, r * np.sin(tt)])
    else:
        raise ValueError(cls)
    return np.column_stack([x, y])


def _perturb(pts, rng, amp: float, n_out: int):
    """Low-frequency smooth displacement, non-uniform resampling, random similarity."""
    size = np.ptp(pts, axis=0).max()
    s = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(pts, axis=0).T))])
    u = s / s[-1]
    disp = np.zeros_like(pts)
    for k in (1, 2, 3):
        disp += (amp * size / k) * np.outer(np.sin(2 * np.pi * k * u + rng.uniform(0, 2 * np.pi)), rng.normal(size=2))
    pts = pts + disp
    # non-uniform parameterization: warped arc-length positions
    w = np.linspace(0, 1, n_out)
    w = w + 0.08 * np.sin(2 * np.pi * w + rng.uniform(0, 2 * np.pi)) * w * (1 - w) * 4
    w = np.sort(np.clip(w, 0, 1))
    s = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(pts, axis=0).T))])
    tgt = w * s[-1]
    pts = np.column_stack([np.interp(tgt, s, pts[:, 0]), np.interp(tgt, s, pts[:, 1])])
    return similarity(pts, rng.uniform(0, 2 * np.pi), rng.uniform(0.5, 3.0), rng.uniform(-10, 10, 2))


def _resample_uniform(pts, m=64):
    return uniform_arclength(pts, m)


def _procrustes(a, b):
    a = a - a.mean(0); b = b - b.mean(0)
    a /= np.sqrt((a ** 2).sum()); b /= np.sqrt((b ** 2).sum())
    U, _, Vt = np.linalg.svd(a.T @ b)
    d = np.sign(np.linalg.det(U @ Vt))
    R = U @ np.diag([1, d]) @ Vt
    return float(np.sqrt(((a @ R - b) ** 2).sum()))


def _turning_function(pts, m=64):
    """Arkin-style signed turning function sampled at m arc-length positions."""
    r = uniform_arclength(pts, m + 1)
    e = np.diff(r, axis=0)
    th = np.unwrap(np.arctan2(e[:, 1], e[:, 0]))
    return th


def e6_classification(per_class=30, seed=7):
    classes = ["alpha", "S", "C", "spiral", "J"]
    results = []
    for level, amp in [("mild", 0.02), ("strong", 0.06)]:
        rng = np.random.default_rng(seed)
        data, labels = [], []
        for ci, cls in enumerate(classes):
            for _ in range(per_class):
                base = _stroke(cls, rng, jitter=0.15 if level == "mild" else 0.30)
                data.append(_perturb(base, rng, amp, int(rng.integers(80, 400))))
                labels.append(ci)
        labels = np.array(labels)
        n = len(data)
        feats = {
            "coords": [_resample_uniform(p) for p in data],
            "turnfn": [_turning_function(p) for p in data],
            "pure": [profile(p) for p in data],
            "hybrid": [profile(p, beta=1.0) for p in data],
        }

        def dist(kind, i, j):
            if kind == "Euclidean, unaligned":
                a, b = feats["coords"][i], feats["coords"][j]
                return float(np.sqrt(((a - a.mean(0) - b + b.mean(0)) ** 2).sum()))
            if kind == "Procrustes (rotation + scale)":
                return _procrustes(feats["coords"][i], feats["coords"][j])
            if kind == "turning function (Arkin et al.)":
                a, b = feats["turnfn"][i], feats["turnfn"][j]
                dd = a - b
                return float(np.sqrt(np.mean((dd - dd.mean()) ** 2)))
            if kind == r"curvature profile $D$ (pure)":
                return profile_distance(feats["pure"][i], feats["pure"][j], n_grid=513)
            if kind == r"curvature profile $D_\beta$ ($\beta=1$)":
                return profile_distance(feats["hybrid"][i], feats["hybrid"][j], n_grid=513)
            raise ValueError(kind)

        kinds = ["Euclidean, unaligned", "Procrustes (rotation + scale)", "turning function (Arkin et al.)",
                 r"curvature profile $D$ (pure)", r"curvature profile $D_\beta$ ($\beta=1$)"]
        for kind in kinds:
            M = np.full((n, n), np.inf)
            for i in range(n):
                for jj in range(i + 1, n):
                    M[i, jj] = M[jj, i] = dist(kind, i, jj)
            pred = labels[np.argmin(M, axis=1)]
            acc = float((pred == labels).mean())
            results.append({"level": level, "method": kind, "accuracy": acc})
    df = pd.DataFrame(results)
    df.to_csv(RESULTS / "e6_classification.csv", index=False)
    piv = df.pivot(index="method", columns="level", values="accuracy").loc[kinds]
    lines = [r"\begin{tabular}{lcc}", r"\toprule", r"descriptor & mild & strong\\", r"\midrule"]
    for k, r in piv.iterrows():
        lines.append(f"{k} & {100*r['mild']:.1f}\\% & {100*r['strong']:.1f}\\%\\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    write_table("tab_classify", "\n".join(lines) + "\n")
    macro("clsPerClass", per_class)
    macro("clsN", per_class * len(classes))
    for lvl, nm in [("mild", "Mild"), ("strong", "Strong")]:
        sub = df[df.level == lvl].set_index("method")["accuracy"]
        macro(f"clsHybrid{nm}", 100 * sub[kinds[4]], ".1f")
        macro(f"clsPure{nm}", 100 * sub[kinds[3]], ".1f")
        macro(f"clsProcrustes{nm}", 100 * sub[kinds[1]], ".1f")
        macro(f"clsTurning{nm}", 100 * sub[kinds[2]], ".1f")
        macro(f"clsEuclid{nm}", 100 * sub[kinds[0]], ".1f")
    return df


# =========================================================================== E7

KS = [20, 30, 50, 75, 100]
N_BENCH = 1000


def methods(closed):
    return {
        "uniform arc length": lambda p, k: uniform_arclength(p, k, closed),
        "Visvalingam--Whyatt": lambda p, k: vw_fixed_k(p, k),
        "RDP (fixed $K$)": lambda p, k: rdp_fixed_k(p, k),
        r"CQ $p{=}1,\beta{=}1$": lambda p, k: quantile_resample(p, k, beta=1.0, p=1.0, closed=closed),
        r"CQ $p{=}\frac12$": lambda p, k: quantile_resample(p, k, beta=None, p=0.5, closed=closed),
        r"CQ $p{=}\frac12$, pinned": lambda p, k: quantile_resample(p, k, beta=None, p=0.5, closed=closed, pin_corners=True),
        r"CQ $p{=}\frac12$, snapped": lambda p, k: quantile_resample(p, k, beta=None, p=0.5, closed=closed, snap_to_vertices=True),
    }


MAIN = r"CQ $p{=}\frac12$, pinned"
RDP = "RDP (fixed $K$)"
OPT = "optimal vertex subset"


def e7_benchmark():
    rows = []
    for fam, (gen, closed) in FAMILIES.items():
        pts = gen(N_BENCH)
        t0 = time.perf_counter()
        opt = {}
        for K in KS:
            simp, eps = optimal_subset_min_eps(pts, K)
            opt[K] = (simp, eps)
        print(f"  optimal DP for {fam}: {time.perf_counter()-t0:.1f}s", flush=True)
        for K in KS:
            for name, fn in methods(closed).items():
                simp = fn(pts, K)
                e = approximation_errors(pts, simp)
                rows.append({"family": fam, "K": K, "method": name, "n_out": len(simp), **e})
            e = approximation_errors(pts, opt[K][0])
            rows.append({"family": fam, "K": K, "method": OPT, "n_out": len(opt[K][0]), **e})
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "e7_benchmark.csv", index=False)
    return df


def e7_tables(df):
    order = list(methods(False)) + [OPT]
    sub = df[df.K == 50].pivot(index="family", columns="method", values="Linf")[order]
    fams = list(FAMILIES)
    sub = sub.loc[fams]
    short = ["uniform", "VW", "RDP", r"CQ\,$1,\beta$", r"CQ\,$\frac12$", r"CQ\,$\frac12$ pin", r"CQ\,$\frac12$ snap", "optimal"]
    lines = [r"\begin{tabular}{l" + "r" * len(order) + "}", r"\toprule",
             "family & " + " & ".join(short) + r"\\", r"\midrule"]
    for fam, r in sub.iterrows():
        vals = r.values
        best = np.nanmin(vals[:-1])  # best non-oracle
        cells = []
        for v in vals[:-1]:
            s = f"{1e3*v:.2f}"
            cells.append(f"\\textbf{{{s}}}" if np.isclose(v, best) else s)
        cells.append(f"\\textit{{{1e3*vals[-1]:.2f}}}")
        lines.append(f"{fam} & " + " & ".join(cells) + r"\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    write_table("tab_benchmark", "\n".join(lines) + "\n")

    # paired statistics, main CQ variant vs RDP over all (family, K)
    piv = df.pivot_table(index=["family", "K"], columns="method", values="Linf")
    ratio = piv[MAIN] / piv[RDP]
    n_pairs = len(ratio)
    gm = float(np.exp(np.log(ratio).mean()))
    rng = np.random.default_rng(0)
    lr = np.log(ratio.values)
    boots = [np.exp(rng.choice(lr, lr.size).mean()) for _ in range(10000)]
    lo, hi = np.percentile(boots, [2.5, 97.5])
    w = wilcoxon(piv[MAIN].values, piv[RDP].values, method="exact")
    wins = int((ratio < 1).sum()); losses = int((ratio > 1).sum())
    macro("benchN", N_BENCH)
    macro("benchPairs", n_pairs)
    macro("benchFamilies", len(FAMILIES))
    macro("benchGM", gm, ".3f")
    macro("benchGMReduction", 100 * (1 - gm), ".0f")
    macro("benchCILo", lo, ".3f")
    macro("benchCIHi", hi, ".3f")
    macro("benchWilcoxonP", sci(w.pvalue))
    macro("benchWilcoxonW", w.statistic, ".0f")
    macro("benchWins", wins)
    macro("benchLosses", losses)
    # per-family losses
    lost = ratio[ratio > 1]
    lost_fams = sorted({f for f, _ in lost.index})
    macro("benchLossFamilies", ", ".join(lost_fams) if lost_fams else "none")

    # ratio to the optimal subset
    r_opt = piv[MAIN] / piv[OPT]
    macro("optRatioGM", float(np.exp(np.log(r_opt).mean())), ".2f")
    macro("optRatioMax", float(r_opt.max()), ".2f")
    macro("optRatioMin", float(r_opt.min()), ".2f")
    r_rdp_opt = piv[RDP] / piv[OPT]
    macro("rdpOptRatioGM", float(np.exp(np.log(r_rdp_opt).mean())), ".2f")
    beats_opt = int((r_opt < 1).sum())
    macro("cqBeatsOpt", beats_opt)
    snap = piv[r"CQ $p{=}\frac12$, snapped"] / piv[RDP]
    macro("snapGM", float(np.exp(np.log(snap).mean())), ".3f")
    r_p1 = piv[r"CQ $p{=}1,\beta{=}1$"] / piv[RDP]
    macro("cqOneGM", float(np.exp(np.log(r_p1).mean())), ".3f")
    r_unpinned = piv[r"CQ $p{=}\frac12$"] / piv[RDP]
    macro("cqUnpinnedGM", float(np.exp(np.log(r_unpinned).mean())), ".3f")
    uni = piv["uniform arc length"] / piv[RDP]
    macro("uniformGM", float(np.exp(np.log(uni).mean())), ".2f")
    vw = piv["Visvalingam--Whyatt"] / piv[RDP]
    macro("vwGM", float(np.exp(np.log(vw).mean())), ".2f")

    # spiral numbers quoted in the drafts (K = 30 and K = 50)
    sp = df[df.family == "F1 log spiral"].pivot(index="K", columns="method", values="Linf")
    for K, nm in [(30, "Thirty"), (50, "Fifty")]:
        macro(f"spiralUniform{nm}", 1e3 * sp.loc[K, "uniform arc length"], ".2f")
        macro(f"spiralRDP{nm}", 1e3 * sp.loc[K, RDP], ".2f")
        macro(f"spiralCQOne{nm}", 1e3 * sp.loc[K, r"CQ $p{=}1,\beta{=}1$"], ".2f")
        macro(f"spiralCQHalf{nm}", 1e3 * sp.loc[K, MAIN], ".2f")
        macro(f"spiralRedUniformCQOne{nm}", 100 * (1 - sp.loc[K, r"CQ $p{=}1,\beta{=}1$"] / sp.loc[K, "uniform arc length"]), ".0f")
        macro(f"spiralRedRDP{nm}", 100 * (1 - sp.loc[K, MAIN] / sp.loc[K, RDP]), ".0f")
    f3 = df[(df.family == "F3 corners+arcs")].pivot(index="K", columns="method", values="Linf")
    macro("cornerUnpinnedVsRDP", 100 * (f3.loc[20, r"CQ $p{=}\frac12$"] / f3.loc[20, RDP] - 1), "+.0f")
    macro("cornerPinnedVsRDP", 100 * (f3.loc[20, MAIN] / f3.loc[20, RDP] - 1), "+.0f")
    return ratio


def fig_benchmark(df):
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.45))
    pts = log_spiral(N_BENCH)
    K = 30
    for a, (name, simp, col) in zip(ax[:2], [
        ("RDP (fixed $K$)", rdp_fixed_k(pts, K), C["orange"]),
        (r"curvature quantiles, $p=1/2$", quantile_resample(pts, K, beta=None, p=0.5, pin_corners=True), C["blue"]),
    ]):
        e = approximation_errors(pts, simp)["Linf"]
        a.plot(*pts.T, color=C["grey"], lw=0.8)
        a.plot(*simp.T, "-o", color=col, ms=2.2, lw=0.9)
        a.set_aspect("equal"); a.set_xticks([]); a.set_yticks([])
        for sp_ in a.spines.values():
            sp_.set_visible(False)
        a.set_title(f"{name}\n$K={K}$, $L_\\infty={1e3*e:.1f}\\times10^{{-3}}$", fontsize=8)
    a = ax[2]
    fams = list(FAMILIES)
    piv = df.pivot_table(index=["family", "K"], columns="method", values="Linf")
    for name, col, mk in [("uniform arc length", C["grey"], "x"), ("Visvalingam--Whyatt", C["green"], "^"),
                          (RDP, C["orange"], "s"), (MAIN, C["blue"], "o"), (OPT, C["black"], "")]:
        gm = [np.exp(np.log(piv.loc[(slice(None), K), name] / piv.loc[(slice(None), K), RDP]).mean()) for K in KS]
        a.plot(KS, gm, marker=mk, ms=3, color=col, ls="--" if name == OPT else "-",
               label={MAIN: r"CQ $p=1/2$", OPT: "optimal subset"}.get(name, name.replace(" (fixed $K$)", "")))
    a.axhline(1, color=C["grey"], lw=0.5)
    a.set_xscale("log"); a.set_xticks(KS); a.set_xticklabels([str(k) for k in KS]); a.minorticks_off()
    a.set_xlabel("vertex budget $K$"); a.set_ylabel(r"$L_\infty$ / RDP (geom. mean)")
    a.set_title(f"(c) all {len(fams)} families, $N={N_BENCH}$", fontsize=8)
    a.legend(frameon=False, fontsize=6.2)
    fig.tight_layout()
    save_fig(fig, "fig4_resampling")


def e7_noise(trials=20, seed=11):
    rng = np.random.default_rng(seed)
    rows = []
    K = 50
    for fam in ["F1 log spiral", "F4 line+bumps"]:
        gen, closed = FAMILIES[fam]
        clean = gen(N_BENCH)
        scale = np.ptp(clean, axis=0).max()
        for eta in [0.0, 0.002, 0.005, 0.01]:
            acc = {}
            for _ in range(trials if eta > 0 else 1):
                noisy = clean + rng.normal(0, eta * scale, clean.shape)
                win = 31
                sm = np.column_stack([savgol_filter(noisy[:, 0], win, 3), savgol_filter(noisy[:, 1], win, 3)]) if eta > 0 else noisy
                for tag, inp in [("raw", noisy), ("smoothed", sm)]:
                    for name, fn in [("RDP", lambda p: rdp_fixed_k(p, K)),
                                     ("CQ", lambda p: quantile_resample(p, K, beta=None, p=0.5, pin_corners=True))]:
                        e = approximation_errors(clean, fn(inp))["Linf"]
                        acc.setdefault(f"{name} {tag}", []).append(e)
            rows.append({"family": fam, "eta": eta, **{k: float(np.mean(v)) for k, v in acc.items()}})
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "e7_noise.csv", index=False)
    lines = [r"\begin{tabular}{llcccc}", r"\toprule",
             r"family & noise $\eta$ & RDP raw & CQ raw & RDP smoothed & CQ smoothed\\", r"\midrule"]
    for _, r in df.iterrows():
        cells = [r["RDP raw"], r["CQ raw"], r["RDP smoothed"], r["CQ smoothed"]]
        lines.append(f"{r['family']} & {100*r['eta']:.1f}\\% & " + " & ".join(f"{1e3*c:.2f}" for c in cells) + r"\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    write_table("tab_noise_resample", "\n".join(lines) + "\n")
    macro("noiseResampleTrials", trials)
    return df


def e7_timing():
    rows = []
    for N in [1000, 10000, 100000]:
        pts = log_spiral(N)
        K = max(20, N // 50)
        for name, fn in [("uniform", lambda p: uniform_arclength(p, K)),
                         ("CQ", lambda p: quantile_resample(p, K, beta=None, p=0.5)),
                         ("RDP", lambda p: rdp_fixed_k(p, K)),
                         ("VW", lambda p: vw_fixed_k(p, K))]:
            ts = []
            for _ in range(5 if N <= 10000 else 2):
                t0 = time.perf_counter(); fn(pts); ts.append(time.perf_counter() - t0)
            rows.append({"N": N, "K": K, "method": name, "ms": 1e3 * float(np.median(ts))})
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "e7_timing.csv", index=False)
    piv = df.pivot(index="N", columns="method", values="ms")[["uniform", "CQ", "RDP", "VW"]]
    lines = [r"\begin{tabular}{rrrrrr}", r"\toprule", r"$N$ & $K$ & uniform & CQ & RDP & VW\\", r"\midrule"]
    for N, r in piv.iterrows():
        K = max(20, N // 50)
        lines.append(f"{N:,} & {K:,} & " + " & ".join(f"{v:.2f}" for v in r.values) + r"\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    write_table("tab_timing", "\n".join(lines).replace(",", "{,}") + "\n")
    return df
