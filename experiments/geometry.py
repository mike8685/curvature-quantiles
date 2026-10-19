"""Experiments E1-E5: analytic landmarks, noise, truncation, multi-transition, distances."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.signal import savgol_filter

import matplotlib.pyplot as plt

from curvquant import (
    curvature_median, point_at, profile, profile_distance,
    profile_distance_quantile_form, quantile_gap, turning_angles,
)
from curvquant.curves import graph, similarity

from common import C, RESULTS, macro, save_fig, sci, write_table

# --------------------------------------------------------------------------- E1


def exact_curvature_quantile(fp, a, b, crit=None, q=0.5):
    """Exact x_q for a graph whose slope angle is monotone on [a,c] and [c,b]."""
    th = lambda x: float(np.arctan(fp(x)))
    pieces = [(a, b)] if crit is None else [(a, crit), (crit, b)]
    var = [abs(th(r) - th(l)) for l, r in pieces]
    target = q * sum(var)
    acc = 0.0
    for (l, r), v in zip(pieces, var):
        if acc + v >= target - 1e-15:
            sgn = np.sign(th(r) - th(l))
            lvl = th(l) + sgn * (target - acc)
            return brentq(lambda x: th(x) - lvl, l, r, xtol=1e-14)
        acc += v
    return b


BENCH = [
    # label, f, f', a, b, critical point of f' (None if monotone slope)
    (r"$x^2$", lambda x: x ** 2, lambda x: 2 * x, 0.0, 2.0, None),
    (r"$x^3$", lambda x: x ** 3, lambda x: 3 * x ** 2, 0.0, 2.0, None),
    (r"$\sqrt{x}$", np.sqrt, lambda x: 0.5 / np.sqrt(x), 1.0, 9.0, None),
    (r"$1-e^{-x}$", lambda x: 1 - np.exp(-x), lambda x: np.exp(-x), 0.0, 4.0, None),
    (r"logistic $1/(1+e^{-x})$", lambda x: 1 / (1 + np.exp(-x)),
     lambda x: np.exp(-x) / (1 + np.exp(-x)) ** 2, -6.0, 6.0, 0.0),
    (r"Gompertz $e^{-e^{-x}}$", lambda x: np.exp(-np.exp(-x)),
     lambda x: np.exp(-x - np.exp(-x)), -6.0, 8.0, 0.0),
]
RES = [21, 51, 101, 201, 501, 1001, 2001]


def e1_analytic():
    rows, conv = [], {}
    for label, f, fp, a, b, c in BENCH:
        exact = exact_curvature_quantile(fp, a, b, c)
        errs = []
        for n in RES:
            errs.append(abs(curvature_median(graph(f, a, b, n))[0] - exact))
        errs = np.array(errs)
        conv[label] = errs
        mask = np.array(RES) >= 101
        nz = mask & (errs > 1e-12)
        order = -np.polyfit(np.log(np.array(RES)[nz]), np.log(errs[nz]), 1)[0] if nz.sum() >= 3 else np.nan
        rows.append({
            "f": label, "a": a, "b": b,
            "theta_a_deg": np.degrees(np.arctan(fp(a))), "theta_b_deg": np.degrees(np.arctan(fp(b))),
            "exact": exact, **{f"err_n{n}": e for n, e in zip(RES, errs)}, "order": order,
        })
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "e1_analytic.csv", index=False)

    # Theorem 1 closed forms used in the text
    macro("medXsq", (np.sqrt(17) - 1) / 8, ".6f")
    macro("medXcube", np.sqrt(np.tan(0.5 * np.arctan(12.0)) / 3), ".6f")
    gomp = df.iloc[5]
    macro("medGompertz", gomp["exact"], ".4f")
    macro("gompertzSlopeEndDeg", gomp["theta_b_deg"], ".3f")

    # landmarks of x^2 on [0,2]
    pts = graph(lambda x: x ** 2, 0, 2, 20001)
    lens = np.hypot(*np.diff(pts, axis=0).T)
    s = np.concatenate([[0], np.cumsum(lens)])
    macro("xsqArcMid", np.interp(0.5 * s[-1], s, pts[:, 0]), ".3f")

    def fmt_err(e):
        if e < 5e-13:
            return r"$<10^{-12}$"
        m, ex = f"{e:.1e}".split("e")
        return rf"${m}\times10^{{{int(ex)}}}$"

    lines = [r"\begin{tabular}{lcrrcrrrc}", r"\toprule",
             r"$f(x)$ & $[a,b]$ & $\theta(a)$ & $\theta(b)$ & exact $x_{1/2}$ & \multicolumn{3}{c}{$|x^{(n)}_{1/2}-x_{1/2}|$} & order\\",
             r"\cmidrule(lr){6-8}", r" & & (deg) & (deg) & & $n=51$ & $n=201$ & $n=2001$ & \\", r"\midrule"]
    for _, r in df.iterrows():
        order = "--" if not np.isfinite(r["order"]) else f"{r['order']:.2f}"
        lines.append(
            f"{r['f']} & $[{r['a']:g},{r['b']:g}]$ & {r['theta_a_deg']:.2f} & {r['theta_b_deg']:.3f} & "
            f"${r['exact']:.6f}$ & {fmt_err(r['err_n51'])} & {fmt_err(r['err_n201'])} & {fmt_err(r['err_n2001'])} & {order}\\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    write_table("tab_analytic", "\n".join(lines) + "\n")
    orders = df["order"].values[[0, 2, 3]]  # x^2, sqrt x, 1-e^-x (x^3 error changes sign)
    macro("convOrderMin", orders.min(), ".2f")
    macro("convOrderMax", orders.max(), ".2f")
    return df, conv


def fig_concept(conv):
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.3))
    pts = graph(lambda x: x ** 2, 0, 2, 2001)
    prof = profile(pts)
    a = ax[0]
    a.plot(pts[:, 0], pts[:, 1], color=C["black"], lw=1.2)
    qs = [0.25, 0.5, 0.75]
    qp = point_at(pts, prof.quantile(qs))
    lens = np.hypot(*np.diff(pts, axis=0).T)
    s = np.concatenate([[0], np.cumsum(lens)])
    arc_mid = point_at(pts, [0.5])[0]
    a.plot(*qp[[0, 2]].T, "o", ms=3.5, color=C["sky"], label=r"curvature quartiles")
    a.plot(*qp[1], "o", ms=5, color=C["blue"], label=r"curvature median")
    a.plot(*arc_mid, "s", ms=4.5, color=C["orange"], label="arc-length midpoint")
    a.plot(1.0, 1.0, "^", ms=4.5, color=C["green"], label="domain midpoint")
    a.set_xlabel("$x$"); a.set_ylabel("$y=x^2$")
    a.set_title("(a) landmarks on $y=x^2$, $[0,2]$")
    a.legend(frameon=False, loc="upper left", handletextpad=0.3)

    a = ax[1]
    uu = np.linspace(0, 1, 400)
    a.plot(uu, profile(pts).cdf(uu), color=C["blue"], label=r"pure, $\beta=0$")
    a.plot(uu, profile(pts, beta=1.0).cdf(uu), color=C["red"], ls="--", label=r"hybrid, $\beta=1$")
    a.plot(uu, uu, color=C["grey"], lw=0.8, ls=":", label="arc length")
    a.axhline(0.5, color=C["grey"], lw=0.5)
    a.set_xlabel(r"normalized arc length $u=s/L$"); a.set_ylabel(r"cumulative mass $F(u)$")
    a.set_title("(b) curvature-mass profile")
    a.legend(frameon=False, loc="lower right")

    a = ax[2]
    cols = [C["blue"], C["orange"], C["green"], C["red"], C["purple"], C["sky"]]
    for (label, *_), col in zip(BENCH, cols):
        e = np.maximum(conv[label], 1e-16)
        if label.startswith("logistic"):
            continue
        a.loglog(RES, e, "o-", ms=2.5, color=col, label=label)
    a.loglog(RES, 1.0 / np.array(RES), color=C["grey"], ls=":", lw=0.8, label=r"$\propto 1/n$")
    a.set_xlabel("vertices $n$"); a.set_ylabel(r"$|x^{(n)}_{1/2}-x_{1/2}|$")
    a.set_title("(c) discrete convergence")
    a.legend(frameon=False, fontsize=6.0, loc="center left", bbox_to_anchor=(1.0, 0.5))
    fig.tight_layout()
    save_fig(fig, "fig1_concept")


# --------------------------------------------------------------------------- E2


def _kappa_argmax_x(pts):
    lens, w = turning_angles(pts)
    k = w / np.maximum(0.5 * (lens[:-1] + lens[1:]), 1e-15)
    return pts[1 + int(np.argmax(k)), 0]


def e2_noise(trials=200, seed=20261019):
    rng = np.random.default_rng(seed)
    x = np.linspace(0, 4, 201)
    yc = 1 - np.exp(-x)
    clean = np.column_stack([x, yc])
    ref = {
        "median_pure": curvature_median(clean)[0],
        "median_hybrid": curvature_median(clean, beta=1.0)[0],
        "kappa_max": _kappa_argmax_x(clean),
    }
    R = np.ptp(yc)
    rows = []
    for lvl in [0.001, 0.005, 0.01, 0.02]:
        e = {k: [] for k in ["raw_pure", "sg_pure", "sg_hybrid", "sg_kmax"]}
        for _ in range(trials):
            yn = yc + rng.normal(0, lvl * R, x.size)
            raw = np.column_stack([x, yn])
            sm = np.column_stack([x, savgol_filter(yn, 51, 3)])
            e["raw_pure"].append(abs(curvature_median(raw)[0] - ref["median_pure"]))
            e["sg_pure"].append(abs(curvature_median(sm)[0] - ref["median_pure"]))
            e["sg_hybrid"].append(abs(curvature_median(sm, beta=1.0)[0] - ref["median_hybrid"]))
            e["sg_kmax"].append(abs(_kappa_argmax_x(sm) - ref["kappa_max"]))
        rows.append({"noise": lvl, **{k: float(np.mean(v)) for k, v in e.items()},
                     **{k + "_sd": float(np.std(v)) for k, v in e.items()}})
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "e2_noise.csv", index=False)
    lines = [r"\begin{tabular}{rcccc}", r"\toprule",
             r"$\sigma/\operatorname{range}(y)$ & median, raw & median, smoothed & hybrid median, smoothed & $\arg\max\kappa$, smoothed\\",
             r"\midrule"]
    for _, r in df.iterrows():
        lines.append(f"{100*r['noise']:.1f}\\% & {r['raw_pure']:.3f} & {r['sg_pure']:.3f} & {r['sg_hybrid']:.3f} & {r['sg_kmax']:.3f}\\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    write_table("tab_noise", "\n".join(lines) + "\n")
    macro("noiseCleanMedian", ref["median_pure"], ".3f")
    macro("noiseCleanHybrid", ref["median_hybrid"], ".3f")
    macro("noiseCleanKmax", ref["kappa_max"], ".3f")
    macro("noiseTrials", trials)
    ratio = df["sg_kmax"] / df["sg_pure"]
    macro("noiseKmaxRatioMin", ratio.min(), ".1f")
    macro("noiseKmaxRatioMax", ratio.max(), ".1f")
    macro("noiseRawWorst", df["raw_pure"].max(), ".2f")
    return df


# --------------------------------------------------------------------------- E3


def e3_truncation():
    f = lambda x: 1 / (1 + np.exp(-x))
    fp = lambda x: np.exp(-x) / (1 + np.exp(-x)) ** 2
    rows = []
    for a, b in [(-6, 6), (-6, 3), (-3, 6), (-2, 4), (-1, 5)]:
        pts = graph(f, a, b, 4001)
        exact = exact_curvature_quantile(fp, a, b, 0.0)
        rows.append({
            "a": a, "b": b,
            "domain_mid": 0.5 * (a + b),
            "arc_mid": point_at(pts, [0.5])[0, 0],
            "median": curvature_median(pts)[0],
            "median_exact": exact,
            "hybrid": curvature_median(pts, beta=1.0)[0],
            "dtheta_deg": np.degrees(abs(np.arctan(fp(b)) - np.arctan(fp(a)))),
        })
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "e3_truncation.csv", index=False)
    lines = [r"\begin{tabular}{crrrrr}", r"\toprule",
             r"window $[a,b]$ & $|\theta(b)-\theta(a)|$ (deg) & domain midpoint & arc-length midpoint & curvature median & hybrid median ($\beta=1$)\\",
             r"\midrule"]
    for _, r in df.iterrows():
        lines.append(f"$[{r['a']},{r['b']}]$ & {r['dtheta_deg']:.2f} & ${r['domain_mid']:+.3f}$ & ${r['arc_mid']:+.3f}$ & ${r['median']:+.3f}$ & ${r['hybrid']:+.3f}$\\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    write_table("tab_truncation", "\n".join(lines) + "\n")
    r = df.iloc[1]
    macro("truncDomainMid", f"{r['domain_mid']:+.2f}")
    macro("truncArcMid", f"{r['arc_mid']:+.2f}")
    macro("truncMedian", f"{r['median']:+.2f}")
    return df


# --------------------------------------------------------------------------- E4


def e4_multi():
    rows, profs = [], {}
    for d in [1, 2, 3, 4, 6]:
        x = np.linspace(-12, 12, 4001)
        y = 1 / (1 + np.exp(-(x + d))) + 1 / (1 + np.exp(-(x - d)))
        pts = np.column_stack([x, y])
        pr = profile(pts)
        profs[d] = pr
        q = [0.25, 0.45, 0.5, 0.55, 0.75]
        xs = point_at(pts, pr.quantile(q))[:, 0]
        rows.append({"d": d, **{f"x_q{int(100*qq)}": v for qq, v in zip(q, xs)},
                     "gap": quantile_gap(pr), "W25": float(np.diff(pr.quantile([0.25, 0.75]))[0])})
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "e4_multi.csv", index=False)
    lines = [r"\begin{tabular}{rrrrrrcc}", r"\toprule",
             r"$d$ & $x(P_{0.25})$ & $x(P_{0.45})$ & $x(P_{0.50})$ & $x(P_{0.55})$ & $x(P_{0.75})$ & gap $P_{0.55}-P_{0.45}$ & $P_{0.75}-P_{0.25}$\\",
             r"\midrule"]
    for _, r in df.iterrows():
        lines.append(f"{int(r['d'])} & ${r['x_q25']:+.2f}$ & ${r['x_q45']:+.2f}$ & ${r['x_q50']:+.2f}$ & ${r['x_q55']:+.2f}$ & ${r['x_q75']:+.2f}$ & {r['gap']:.3f} & {r['W25']:.3f}\\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    write_table("tab_multi", "\n".join(lines) + "\n")
    macro("gapDOne", df.iloc[0]["gap"], ".3f")
    macro("gapDSix", df.iloc[-1]["gap"], ".3f")
    macro("gapRatio", df.iloc[-1]["gap"] / df.iloc[0]["gap"], ".1f")
    macro("quartileDSix", abs(df.iloc[-1]["x_q25"]), ".2f")
    return df, profs


def fig_multi(profs):
    fig, ax = plt.subplots(1, 2, figsize=(5.6, 2.2))
    cmap = plt.get_cmap("viridis")
    uu = np.linspace(0, 1, 800)
    for i, (d, pr) in enumerate(profs.items()):
        ax[0].plot(uu, pr.cdf(uu), color=cmap(i / (len(profs) - 1) * 0.9), label=f"$d={d}$")
    ax[0].axhspan(0.45, 0.55, color=C["yellow"], alpha=0.35, lw=0)
    ax[0].set_xlabel(r"$u=s/L$"); ax[0].set_ylabel(r"$F(u)$")
    ax[0].set_title("(a) profiles of the double sigmoid")
    ax[0].legend(frameon=False, fontsize=6.5)
    ds = list(profs)
    ax[1].plot(ds, [quantile_gap(profs[d]) for d in ds], "o-", color=C["purple"], ms=3, label=r"$P_{0.55}-P_{0.45}$")
    ax[1].plot(ds, [float(np.diff(profs[d].quantile([0.25, 0.75]))[0]) for d in ds], "s-", color=C["sky"], ms=3, label=r"$P_{0.75}-P_{0.25}$")
    ax[1].set_xlabel("separation $d$"); ax[1].set_ylabel(r"span in $u$")
    ax[1].set_title("(b) quantile gaps")
    ax[1].legend(frameon=False)
    fig.tight_layout()
    save_fig(fig, "fig2_multi")


# --------------------------------------------------------------------------- E5


def archetypes(n=601):
    x = np.linspace(0, 2, n)
    early = np.column_stack([x, 1 - np.exp(-4 * x)])
    out = {
        "early bend": early,
        "early, rotated+scaled": similarity(early, np.radians(37), 2.5, (5.0, -3.0)),
        "early, resampled": np.column_stack([np.linspace(0, 2, 157) ** 1.0, 1 - np.exp(-4 * np.linspace(0, 2, 157))]),
        "late bend": np.column_stack([x, np.exp(4 * (x - 2))]),
        "single sigmoid": np.column_stack([np.linspace(-4, 4, n), 1 / (1 + np.exp(-3 * np.linspace(-4, 4, n)))]),
        "double sigmoid": np.column_stack([np.linspace(-8, 8, n), 0.5 * (1 / (1 + np.exp(-3 * (np.linspace(-8, 8, n) + 3.5))) + 1 / (1 + np.exp(-3 * (np.linspace(-8, 8, n) - 3.5))))]),
        "circular arc": np.column_stack([np.cos(np.linspace(0, np.pi / 2, n)), np.sin(np.linspace(0, np.pi / 2, n))]),
    }
    return out


def curve_from_curvature(kappa, s):
    ds = s[1] - s[0]
    th = np.concatenate([[0], np.cumsum(0.5 * (kappa[1:] + kappa[:-1]) * ds)])
    x = np.concatenate([[0], np.cumsum(0.5 * (np.cos(th[1:]) + np.cos(th[:-1])) * ds)])
    y = np.concatenate([[0], np.cumsum(0.5 * (np.sin(th[1:]) + np.sin(th[:-1])) * ds)])
    return np.column_stack([x, y])


def e5_distances():
    arch = archetypes()
    names = list(arch)
    P = {k: profile(v) for k, v in arch.items()}
    Hb = {k: profile(v, beta=1.0) for k, v in arch.items()}
    Dp = np.array([[profile_distance(P[a], P[b]) for b in names] for a in names])
    Dh = np.array([[profile_distance(Hb[a], Hb[b]) for b in names] for a in names])
    pd.DataFrame(Dp, index=names, columns=names).to_csv(RESULTS / "e5_distance_pure.csv")
    pd.DataFrame(Dh, index=names, columns=names).to_csv(RESULTS / "e5_distance_hybrid.csv")
    # W1 identity: CDF form vs quantile form
    dq = np.array([[profile_distance_quantile_form(P[a], P[b]) for b in names] for a in names])
    macro("wOneIdentityMaxDiff", sci(np.abs(dq - Dp).max()))
    i0, i1, i2 = names.index("early bend"), names.index("early, rotated+scaled"), names.index("early, resampled")
    macro("invarianceSimilarity", sci(Dp[i0, i1]))
    macro("invarianceResample", Dp[i0, i2], ".4f")
    il = names.index("late bend")
    macro("distEarlyLate", Dp[i0, il], ".3f")
    macro("distEarlyLateHybrid", Dh[i0, il], ".3f")

    # U vs S: identical unsigned curvature, different shapes
    s = np.linspace(0, 2 * np.pi, 2001)
    U = curve_from_curvature(0.45 * np.abs(np.sin(s)), s)
    S = curve_from_curvature(0.45 * np.sin(s), s)
    dUS = profile_distance(profile(U), profile(S))
    macro("distUS", sci(dUS))

    # nearly straight line with a 0.15% dimple
    x = np.linspace(0, 1, 401)
    base = np.column_stack([x, 0.002 * x ** 2])
    bump = np.column_stack([x, 0.002 * x ** 2 + 0.0015 * np.exp(-((x - 0.15) / 0.03) ** 2)])
    dpure = profile_distance(profile(base), profile(bump))
    dhyb = profile_distance(profile(base, beta=1.0), profile(bump, beta=1.0))
    macro("dimpleTurningBase", profile(base).K, ".4f")
    macro("dimplePure", dpure, ".3f")
    macro("dimpleHybrid", dhyb, ".4f")
    macro("dimpleReduction", 100 * (1 - dhyb / dpure), ".0f")

    # table of the distance matrix (pure / hybrid), main archetypes only
    show = ["early bend", "early, rotated+scaled", "early, resampled", "late bend", "single sigmoid", "double sigmoid", "circular arc"]
    idx = [names.index(k) for k in show]
    head = " & ".join([f"({i+1})" for i in range(len(show))])
    lines = [r"\begin{tabular}{l" + "c" * len(show) + "}", r"\toprule", f" & {head}\\\\", r"\midrule"]
    for r, ir in enumerate(idx):
        cells = []
        for c, ic in enumerate(idx):
            if c < r:
                cells.append(f"{Dh[ir, ic]:.3f}")
            elif c == r:
                cells.append("--")
            else:
                cells.append(f"\\textbf{{{Dp[ir, ic]:.3f}}}")
        lines.append(f"({r+1}) {show[r]} & " + " & ".join(cells) + r"\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    write_table("tab_distance", "\n".join(lines) + "\n")
    return {"U": U, "S": S, "base": base, "bump": bump, "P": P}


def fig_distance(ctx):
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.3))
    q = np.linspace(0, 1, 500)
    P = ctx["P"]
    style = {"early bend": (C["blue"], "-"), "late bend": (C["red"], "-"), "single sigmoid": (C["green"], "--"),
             "double sigmoid": (C["purple"], "-."), "circular arc": (C["black"], ":")}
    for k, (col, ls) in style.items():
        ax[0].plot(q, P[k].quantile(q), color=col, ls=ls, label=k)
    ax[0].fill_between(q, P["early bend"].quantile(q), P["late bend"].quantile(q), color=C["sky"], alpha=0.25, lw=0)
    ax[0].set_xlabel("quantile $q$"); ax[0].set_ylabel(r"$P_\gamma(q)$")
    ax[0].set_title(r"(a) profiles; shaded area $=D$")
    ax[0].legend(frameon=False, fontsize=6.2)
    ax[1].plot(*ctx["U"].T, color=C["blue"], label=r"$\kappa=0.45|\sin s|$")
    ax[1].plot(*ctx["S"].T, color=C["red"], ls="--", label=r"$\kappa=0.45\sin s$")
    ax[1].set_aspect("equal"); ax[1].set_title("(b) same unsigned profile, $D=0$")
    ax[1].legend(frameon=False, fontsize=6.5, loc="lower center")
    for kind, beta, col in [("pure", None, C["red"]), ("hybrid", 1.0, C["green"])]:
        pb, pk = profile(ctx["base"], beta=beta), profile(ctx["bump"], beta=beta)
        ax[2].plot(q, pb.quantile(q), color=col, lw=1.0, label=f"{kind}: base")
        ax[2].plot(q, pk.quantile(q), color=col, ls="--", lw=1.0, label=f"{kind}: dimpled")
    ax[2].set_xlabel("quantile $q$"); ax[2].set_ylabel(r"$P(q)$")
    ax[2].set_title("(c) near-straight line, 0.15% dimple")
    ax[2].legend(frameon=False, fontsize=6.2)
    fig.tight_layout()
    save_fig(fig, "fig3_distance")
