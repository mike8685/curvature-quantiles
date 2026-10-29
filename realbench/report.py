"""Automated report: every pre-registered result, wins and losses alike (protocol 5)."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from . import config as C
from .runner import RESULTS
from .stats import compare, decide, holm, log_ratios


def _dev_run(dataset: str):
    d = RESULTS / "dev" / dataset
    if not (d / "rows.csv").exists():
        return None, None, []
    meta = json.loads((d / "meta.json").read_text())
    tf = d / "tuned.json"
    meta["tuned"] = json.loads(tf.read_text())["chosen"] if tf.exists() else {}
    return pd.read_csv(d / "rows.csv"), meta, []


def _test_run(dataset: str):
    """The reported test run is the FIRST logged, completed run whose rows match their logged hash.

    Rows that were never logged are refused; every other logged run is listed (protocol 4).
    """
    from .lock import runlog
    from .runner import sha_file
    log = [r for r in runlog() if r.get("dataset") == dataset and r.get("split") == "test"]
    starts = {r["run_id"]: r for r in log if r.get("event") == "start"}
    dones = {r["run_id"]: r for r in log if r.get("event") == "done"}
    notes, chosen = [], None
    for run_id, st in starts.items():
        d = RESULTS / "test" / dataset / run_id
        dn = dones.get(run_id)
        if dn is None or not (d / "rows.csv").exists():
            notes.append(f"{run_id}: started {st['at']}, no completed output")
            continue
        if sha_file(d / "rows.csv") != dn.get("rows_sha256"):
            notes.append(f"{run_id}: rows.csv does not match its logged hash; REFUSED")
            continue
        flag = (" (allow-unpushed; public-check issues: " + "; ".join(st.get("public_issues") or []) + ")"
                if st.get("allow_unpushed") else "")
        if chosen is None:
            chosen = (run_id, st)
            notes.append(f"{run_id}: REPORTED run, started {st['at']}{flag}")
        else:
            notes.append(f"{run_id}: additional run started {st['at']}{flag}; not the reported run")
    on_disk = {p.name for p in (RESULTS / "test" / dataset).glob("run-*")} if (RESULTS / "test" / dataset).exists() else set()
    for extra in sorted(on_disk - set(starts)):
        notes.append(f"{extra}: output on disk with no run-log record; REFUSED")
    if chosen is None:
        return None, None, notes
    run_id, st = chosen
    d = RESULTS / "test" / dataset / run_id
    meta = json.loads((d / "meta.json").read_text())
    meta["confirmatory"] = bool(st.get("confirmatory"))          # as recorded at test start
    return pd.read_csv(d / "rows.csv"), meta, notes


def _fmt(r: dict) -> str:
    if not np.isfinite(r.get("gmr", np.nan)):
        return "n/a | n/a | n/a | n/a | n/a | 0 / 0 | n/a"
    return (f"{r['gmr']:.3f} | [{r['ci_lo']:.3f}, {r['ci_hi']:.3f}] | [{r['ci_bonf_lo']:.3f}, {r['ci_bonf_hi']:.3f}] | "
            f"{r['p']:.2g} | {100*r['win_rate']:.0f}% / {100*r['loss_rate']:.0f}% | "
            f"{r['n_clusters']} / {r['n_curves']} ({r['n_zero_clusters']} tied) | "
            f"{r['exact_both']} / {r['exact_a_only']} / {r['exact_b_only']}")


HEAD = ("GMR | 95% CI | Bonferroni CI | p | win / loss | clusters / curves | exact fits both / A only / B only")


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
    lines = [f"# Real-data benchmark report ({split} split)", "",
             f"Protocol {C.PROTOCOL_VERSION}. Primary metric: directed Hausdorff error. GMR < 1 means the first "
             f"method has lower error. Per-test alpha (Bonferroni over {C.N_PLANNED_TESTS} planned tests) = "
             f"{C.ALPHA_PER_TEST:.5f}; better/worse decisions use the Bonferroni-level CI, equivalence uses the "
             f"1 - 2 alpha CI within [1/{C.EQUIVALENCE_MARGIN}, {C.EQUIVALENCE_MARGIN}]. Errors below "
             f"{C.ERROR_FLOOR:g} x diameter count as exact fits and tie.", ""]
    if split != "test":
        lines += ["**DEV SPLIT. Used only for tuning and debugging. Not evidence.**", ""]
    base = RESULTS / split
    names = sorted({p.name for p in base.glob("*") if p.is_dir()}) if base.exists() else []
    results, pvals = {}, {}
    for ds in names:
        df, meta, notes = _dev_run(ds) if split == "dev" else _test_run(ds)
        if df is None or not len(df):
            if notes:
                results[ds] = (None, None, notes)
            continue
        res = analyse(df, meta.get("tuned") or {})
        conf = bool(meta.get("confirmatory")) and split == "test"
        res["confirmatory_dataset"] = conf
        res["power_ok"] = all(r["n_clusters"] >= C.MIN_TEST_CLUSTERS for r in res["confirmatory"].values())
        res["pin_inactive"] = float(1 - df[(df.method == "cq") & (df.c == 0)].pin_active.mean()) if "pin_active" in df else float("nan")
        results[ds] = (res, meta, notes)
        if conf:
            for hid, r in res["confirmatory"].items():
                pvals[f"{ds}/{hid}"] = r["p"]
    adj = holm(pvals, C.N_PLANNED_TESTS)
    for ds, (res, meta, notes) in results.items():
        if res is None:
            lines += [f"## {ds}", "", "Runs: " + "; ".join(notes), ""]
            continue
        status = ("confirmatory" if res["confirmatory_dataset"] and res["power_ok"]
                  else "exploratory (power gate: < 30 clusters)" if res["confirmatory_dataset"] else "exploratory")
        ls = meta.get("load_stats", {})
        lines += [f"## {ds} ({status})", ""]
        if notes:
            lines += ["Runs: " + "; ".join(notes), ""]
        lines += [f"{meta['n_curves']} curves; eligible dev/test = {ls.get('eligible')}; clusters dev/test = "
                  f"{ls.get('clusters')}; too short {ls.get('too_short')}, duplicates {ls.get('duplicates')}; "
                  f"loader skips {ls.get('skips') or 'none'}; tuned c = {meta.get('tuned')}", "",
                  f"| test | comparison | {HEAD} | Holm p | decision |",
                  "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
        for hid, r in res["confirmatory"].items():
            dec = decide(r) if status == "confirmatory" else f"({decide(r)}; exploratory)"
            hp = adj.get(f"{ds}/{hid}")
            lines.append(f"| {hid} | {r['a']} vs {r['b']} ({r['mode']}) | {_fmt(r)} | "
                         f"{'' if hp is None else f'{hp:.2g}'} | {dec} |")
        for title, key in [("Secondary metrics (same comparisons)", "secondary"),
                           ("Ablations and baselines vs RDP (c = 0)", "ablations"),
                           ("Ratio to optimal vertex subset (shortcut error, N <= 800)", "oracle")]:
            if res[key]:
                lines += ["", title + ":", "", f"| comparison | {HEAD} |", "| --- " * 8 + "|"]
                lines += [f"| {k} | {_fmt(r)} |" for k, r in res[key].items()]
        lines += ["", "GMR by actual budget fraction K/N (clusters): " + "; ".join(
            f"{h}: " + ", ".join(f"{k} {v:.3f} ({n})" for k, (v, n) in d.items()) for h, d in res["per_K"].items()), "",
                  "Spearman correlation of H1 per-curve log ratio with curve properties (positive = CQ relatively "
                  "worse as the property grows): " + ", ".join(f"{k} {v:+.2f}" for k, v in res.get("correlates", {}).items()), "",
                  f"Failure rates (fallback to uniform): {res['failures'] or 'none'}", "",
                  f"Mean point shortfall below K: {res['shortfall'] or 'none'}", "",
                  f"Share of cq runs (c = 0) where corner pinning switched off (no corner, or more than K-2): "
                  f"{res['pin_inactive']:.1%}", "",
                  "Median runtime ms (Python, indicative): " + ", ".join(
                      f"{k} {v:.2f}" for k, v in res["runtime_ms"].items()), ""]
    missing = [d for d in C.CONFIRMATORY_DATASETS if d not in results or results[d][0] is None]
    if missing and split == "test":
        lines += [f"Planned confirmatory datasets without a reported run: {', '.join(missing)} "
                  "(their tests count as p = 1 in Holm).", ""]
    return "\n".join(lines)
