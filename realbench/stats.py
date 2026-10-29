"""Paired, cluster-level statistics (protocol 5.3 and 5.4)."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from . import config as C


def log_ratios(df: pd.DataFrame, a: str, b: str, ca: str | int, cb: str | int, metric: str,
               tuned: dict | None = None) -> pd.DataFrame:
    """Per (curve, K) log error ratio A/B. ca/cb: smoothing level or 'tuned'.

    Both errors are first raised to the floor tau = ERROR_FLOOR x diameter, so two exact
    fits (errors at rounding level) tie at lr = 0 instead of producing a huge ratio of
    rounding noise. ``exact_a`` / ``exact_b`` flag errors at or below the floor.
    """
    def pick(m, c):
        if c == "tuned":
            c = int((tuned or {}).get(m, 0))
        return df[(df.method == m) & (df.c == int(c))]
    A = pick(a, ca).set_index(["curve_id", "K"])
    B = pick(b, cb).set_index(["curve_id", "K"])
    j = A[[metric, "diam", "cluster"]].join(B[[metric]], rsuffix="_b", how="inner")
    tau = C.ERROR_FLOOR * j["diam"]
    ea = np.maximum(j[metric], tau)
    eb = np.maximum(j[metric + "_b"], tau)
    return pd.DataFrame({"cluster": j["cluster"], "lr": np.log(ea / eb),
                         "exact_a": j[metric] <= tau, "exact_b": j[metric + "_b"] <= tau}).reset_index()


def cluster_scores(lr: pd.DataFrame) -> pd.Series:
    curve = lr.groupby(["cluster", "curve_id"])["lr"].mean()
    return curve.groupby("cluster").mean()


def _ci(boots: np.ndarray, level: float) -> tuple[float, float]:
    lo, hi = np.percentile(boots, [100 * (1 - level) / 2, 100 * (1 + level) / 2])
    return float(lo), float(hi)


def compare(lr: pd.DataFrame, seed: int = C.SEED) -> dict:
    """GMR with three bootstrap intervals over clusters:
    ci_*      95% (for display),
    ci_bonf_* 1 - alpha_per_test (used for the better/worse decisions),
    ci_tost_* 1 - 2 alpha_per_test (two one-sided tests, used for equivalence).
    """
    s = cluster_scores(lr).values
    n = len(s)
    out = {"n_clusters": n, "n_curves": int(lr.curve_id.nunique()) if len(lr) else 0, "n_pairs": len(lr),
           "n_zero_clusters": int((s == 0).sum()),
           "exact_both": int((lr.exact_a & lr.exact_b).sum()) if len(lr) else 0,
           "exact_a_only": int((lr.exact_a & ~lr.exact_b).sum()) if len(lr) else 0,
           "exact_b_only": int((~lr.exact_a & lr.exact_b).sum()) if len(lr) else 0}
    nan = float("nan")
    if n == 0:
        return out | {k: nan for k in ("gmr", "ci_lo", "ci_hi", "ci_bonf_lo", "ci_bonf_hi", "ci_tost_lo",
                                       "ci_tost_hi", "p", "win_rate", "loss_rate")}
    rng = np.random.default_rng(seed)
    boots = np.exp(s[rng.integers(0, n, size=(C.BOOTSTRAP_RESAMPLES, n))].mean(1))
    a = C.ALPHA_PER_TEST
    lo, hi = _ci(boots, 0.95)
    blo, bhi = _ci(boots, 1 - a)
    tlo, thi = _ci(boots, 1 - 2 * a)
    nz = s[s != 0]
    p = 1.0 if len(nz) == 0 else float(wilcoxon(nz, zero_method="wilcox", alternative="two-sided").pvalue)
    return out | {"gmr": float(np.exp(s.mean())), "ci_lo": lo, "ci_hi": hi, "ci_bonf_lo": blo, "ci_bonf_hi": bhi,
                  "ci_tost_lo": tlo, "ci_tost_hi": thi, "p": p,
                  "win_rate": float((s < 0).mean()), "loss_rate": float((s > 0).mean())}


def decide(r: dict, alpha: float = C.ALPHA_PER_TEST) -> str:
    if not np.isfinite(r.get("gmr", np.nan)):
        return "no data"
    m = C.EQUIVALENCE_MARGIN
    if r["p"] < alpha and r["ci_bonf_hi"] < 1:
        return "CQ better"
    if r["p"] < alpha and r["ci_bonf_lo"] > 1:
        return "CQ worse"
    if r["ci_tost_lo"] >= 1 / m and r["ci_tost_hi"] <= m:
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
