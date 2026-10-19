"""Regenerate every table, figure and number in the paper.

    python experiments/run_all.py            # full run (a few minutes)

Outputs: results/*.csv, results/numbers.json, paper/figures/*, paper/generated/*.tex
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import geometry as G  # noqa: E402
import resampling as R  # noqa: E402
from common import write_macros  # noqa: E402


def main() -> None:
    t0 = time.perf_counter()
    print("E1 analytic landmarks"); _, conv = G.e1_analytic(); G.fig_concept(conv)
    print("E2 noise"); G.e2_noise()
    print("E3 truncation"); G.e3_truncation()
    print("E4 multi-transition"); _, profs = G.e4_multi(); G.fig_multi(profs)
    print("E5 profile distances"); ctx = G.e5_distances(); G.fig_distance(ctx)
    print("E6 stroke classification"); R.e6_classification()
    print("E7 simplification benchmark"); df = R.e7_benchmark(); R.e7_tables(df); R.fig_benchmark(df)
    print("E7 noise"); R.e7_noise()
    print("E7 timing"); R.e7_timing()
    write_macros()
    print(f"done in {time.perf_counter() - t0:.0f}s")


if __name__ == "__main__":
    main()
