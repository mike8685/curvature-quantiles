"""Paired, cluster-level statistics (protocol 5.3 and 5.4)."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from . import config as C


def log_ratios(df: pd.DataFrame, a: str, b: str, ca: str | int, cb: str | int, metric: str,
               tuned: dict | None = None) -> pd.DataFrame:
    """Per (curve, K) log error ratio A/B. ca/cb: smoothing level or 'tuned'."""
    def pick(m, c):
        if c == "tuned":
            c = int((tuned or {}).get(m, 0))
        return df[(df.method == m) & (df.c == int(c))]
    A = pick(a, ca).set_index(["curve_id", "K"])
    B = pick(b, cb).set_index(["curve_id", "K"])
    j = A[[metric, "diam", "cluster"]].join(B[[metric]], rsuffix="_b", how="inner")
    floor = C.ZERO_FLOOR * j["diam"]
    ea = j[metric].where(j[metric] > 0, floor)
    eb = j[metric + "_b"].where(j[metric + "_b"] > 0, floor)
    both0 = (j[metric] <= 0) & (j[metric + "_b"] <= 0)
    lr = np.log(ea / eb).where(~both0, 0.0)
    return pd.DataFrame({"cluster": j["cluster"], "lr": lr}).reset_index()


def cluster_scores(lr: pd.DataFrame) -> pd.Series:
    curve = lr.groupby(["cluster", "curve_id"])["lr"].mean()
    return curve.groupby("cluster").mean()


def compare(lr: pd.DataFrame, seed: int = C.SEED) -> dict:
    s = cluster_scores(lr).values
    n = len(s)
    out = {"n_clusters": n, "n_curves": int(lr.curve_id.nunique()), "n_pairs": len(lr)}
    if n == 0:
        return out | {"gmr": np.nan, "ci_lo": np.nan, "ci_hi": np.nan, "p": np.nan, "win_rate": np.nan}
    gmr = float(np.exp(s.mean()))
    rng = np.random.default_rng(seed)
    boots = np.exp(s[rng.integers(0, n, size=(C.BOOTSTRAP_RESAMPLES, n))].mean(1))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    nz = s[s != 0]
    p = 1.0 if len(nz) == 0 else float(wilcoxon(nz, zero_method="wilcox", alternative="two-sided").pvalue)
    return out | {"gmr": gmr, "ci_lo": float(lo), "ci_hi": float(hi), "p": p,
                  "win_rate": float((s < 0).mean()), "loss_rate": float((s > 0).mean())}


def decide(r: dict, alpha: float = C.ALPHA_PER_TEST) -> str:
    if not np.isfinite(r.get("gmr", np.nan)):
        return "no data"
    m = C.EQUIVALENCE_MARGIN
    if r["p"] < alpha and r["ci_hi"] < 1:
        return "CQ better"
    if r["p"] < alpha and r["ci_lo"] > 1:
        return "CQ worse"
    if r["ci_lo"] >= 1 / m and r["ci_hi"] <= m:
        return "practically equivalent"
    return "inconclusive"


def holm(pvals: dict, m_total: int) -> dict:
    """Holm step-down over m_total planned tests; tests not run count as p = 1."""
    items = sorted(pvals.items(), key=lambda kv: kv[1])
    adj, running = {}, 0.0
    for i, (k, p) in enumerate(items):
        running = max(running, min(1.0, (m_total - i) * p))
        adj[k] = running
    return adj
