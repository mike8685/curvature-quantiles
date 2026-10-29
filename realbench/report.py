"""Automated report: every pre-registered result, wins and losses alike (protocol 5)."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from . import config as C
from .lock import load_registry
from .runner import RESULTS
from .stats import compare, decide, holm, log_ratios


def _load(split: str, dataset: str):
    d = RESULTS / split / dataset
    if not (d / "rows.csv").exists():
        return None, None
    return pd.read_csv(d / "rows.csv"), json.loads((d / "meta.json").read_text())


def _fmt(r: dict) -> str:
    if not np.isfinite(r.get("gmr", np.nan)):
        return "n/a | n/a | n/a | n/a | 0 / 0"
    return (f"{r['gmr']:.3f} | [{r['ci_lo']:.3f}, {r['ci_hi']:.3f}] | {r['p']:.2g} | "
            f"{100*r['win_rate']:.0f}% | {r['n_clusters']} / {r['n_curves']}")


def analyse(df: pd.DataFrame, tuned: dict) -> dict:
    out = {"confirmatory": {}, "secondary": {}, "ablations": {}, "per_K": {}, "oracle": {}, "failures": {},
           "runtime_ms": {}, "shortfall": {}, "correlates": {}}
    for hid, a, b, mode in C.CONFIRMATORY_COMPARISONS:
        c = 0 if mode == "c0" else "tuned"
        out["confirmatory"][hid] = compare(log_ratios(df, a, b, c, c, C.PRIMARY_METRIC, tuned)) | {"a": a, "b": b, "mode": mode}
        for met in C.SECONDARY_METRICS:
            out["secondary"][f"{hid}/{met}"] = compare(log_ratios(df, a, b, c, c, met, tuned))
        feats = df[["curve_id", "K", "N", "total_turning", "straight_frac"]].drop_duplicates(["curve_id", "K"])
        lr = log_ratios(df, a, b, c, c, C.PRIMARY_METRIC, tuned).merge(feats, on=["curve_id", "K"])
        # actual budget fraction K/N (the K_MIN floor makes the nominal r misleading on short curves)
        lr["kn"] = pd.cut(lr.K / lr.N, [0, 0.03, 0.07, 0.15, 1.0], labels=["<=3%", "3-7%", "7-15%", ">15%"])
        out["per_K"][hid] = {str(k): (float(np.exp(g.groupby("cluster").lr.mean().mean())), int(g.cluster.nunique()))
                             for k, g in lr.groupby("kn", observed=True)}
        if hid == "H1":
            cur = lr.groupby("curve_id").agg(lr=("lr", "mean"), N=("N", "first"), turning=("total_turning", "first"),
                                             straight=("straight_frac", "first"), kn=("K", "mean"))
            cur["kn"] = cur.kn / cur.N
            out["correlates"] = {f: float(spearmanr(cur.lr, cur[f]).statistic) for f in ("N", "turning", "straight", "kn")
                                 if cur[f].nunique() > 1}
    from .methods import ABLATIONS
    for m in ("vw",) + tuple(ABLATIONS):
        out["ablations"][m] = compare(log_ratios(df, m, "rdp", 0, 0, C.PRIMARY_METRIC))
    for m in ("cq", "rdp", "vw", "cq_snapped"):
        lr = log_ratios(df, m, "oracle", 0, 0, "shortcut")
        if len(lr):
            out["oracle"][m] = compare(lr)
    f = df.groupby(["method", "c"]).failed.mean()
    out["failures"] = {f"{m}@c{c}": float(v) for (m, c), v in f.items() if v > 0}
    sf = (df.K - df.n_out).clip(lower=0).groupby([df.method, df.c]).mean()
    out["shortfall"] = {f"{m}@c{c}": float(v) for (m, c), v in sf.items() if v > 0}
    out["runtime_ms"] = {f"{m}@c{c}": float(v) for (m, c), v in df.groupby(["method", "c"]).ms.median().items()}
    return out


def build_report(split: str = "test") -> str:
    reg = load_registry()
    lines = [f"# Real-data benchmark report ({split} split)", "",
             f"Protocol {C.PROTOCOL_VERSION}. Primary metric: directed Hausdorff error. "
             f"GMR < 1 means CQ has lower error. Per-test alpha (Bonferroni over {C.N_PLANNED_TESTS} planned tests) = "
             f"{C.ALPHA_PER_TEST:.5f}.", ""]
    if split != "test":
        lines += ["**DEV SPLIT. Used only for tuning and debugging. Not evidence.**", ""]
    results, pvals = {}, {}
    for ds in sorted({p.name for p in (RESULTS / split).glob("*") if p.is_dir()}):
        df, meta = _load(split, ds)
        if df is None or not len(df):
            continue
        tuned = (meta.get("tuned") or {}) if split == "test" else {}
        res = analyse(df, tuned)
        conf = reg.get(ds, {}).get("confirmatory", False) and split == "test"
        res["confirmatory_dataset"] = conf
        res["power_ok"] = all(r["n_clusters"] >= C.MIN_TEST_CLUSTERS for r in res["confirmatory"].values())
        results[ds] = (res, meta)
        if conf:
            for hid, r in res["confirmatory"].items():
                pvals[f"{ds}/{hid}"] = r["p"]
    adj = holm(pvals, C.N_PLANNED_TESTS)
    for ds, (res, meta) in results.items():
        status = ("confirmatory" if res["confirmatory_dataset"] and res["power_ok"]
                  else "exploratory (power gate: < 30 clusters)" if res["confirmatory_dataset"] else "exploratory")
        lines += [f"## {ds} ({status})", "",
                  f"{meta['n_curves']} curves; eligible dev/test = {meta['load_stats'].get('eligible')}; "
                  f"clusters dev/test = {meta['load_stats'].get('clusters')}; tuned c = {meta.get('tuned')}", "",
                  "| test | comparison | GMR | 95% CI | p | win rate | clusters / curves | Holm p | decision |",
                  "| --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
        for hid, r in res["confirmatory"].items():
            dec = decide(r) if status == "confirmatory" else f"({decide(r)}; exploratory)"
            hp = adj.get(f"{ds}/{hid}")
            lines.append(f"| {hid} | {r['a']} vs {r['b']} ({r['mode']}) | {_fmt(r)} | "
                         f"{'' if hp is None else f'{hp:.2g}'} | {dec} |")
        lines += ["", "Secondary metrics (same comparisons):", "",
                  "| comparison / metric | GMR | 95% CI | p | win rate | clusters / curves |", "| --- | --- | --- | --- | --- | --- |"]
        lines += [f"| {k} | {_fmt(r)} |" for k, r in res["secondary"].items()]
        lines += ["", "Ablations and baselines vs RDP (c = 0):", "",
                  "| method / rdp | GMR | 95% CI | p | win rate | clusters / curves |", "| --- | --- | --- | --- | --- | --- |"]
        lines += [f"| {k} | {_fmt(r)} |" for k, r in res["ablations"].items()]
        if res["oracle"]:
            lines += ["", "Ratio to optimal vertex subset (shortcut error, N <= 800):", "",
                      "| method / oracle | GMR | 95% CI | p | win rate | clusters / curves |", "| --- | --- | --- | --- | --- | --- |"]
            lines += [f"| {k} | {_fmt(r)} |" for k, r in res["oracle"].items()]
        lines += ["", "GMR by actual budget fraction K/N (clusters): " + "; ".join(
            f"{h}: " + ", ".join(f"{k} {v:.3f} ({n})" for k, (v, n) in d.items()) for h, d in res["per_K"].items()), "",
                  "Spearman correlation of H1 per-curve log ratio with curve properties (positive = CQ relatively "
                  "worse as the property grows): " + ", ".join(f"{k} {v:+.2f}" for k, v in res.get("correlates", {}).items()), "",
                  f"Failure rates (fallback to uniform): {res['failures'] or 'none'}", "",
                  f"Mean point shortfall below K: {res['shortfall'] or 'none'}", "",
                  "Median runtime ms (Python, indicative): " + ", ".join(
                      f"{k} {v:.2f}" for k, v in res["runtime_ms"].items()), ""]
    missing = [d for d in C.FIXED_CONFIRMATORY if d not in results]
    if missing and split == "test":
        lines += [f"Planned confirmatory datasets not yet run: {', '.join(missing)} (their tests count as p = 1 in Holm).", ""]
    return "\n".join(lines)
