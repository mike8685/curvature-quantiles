"""Command line: python -m realbench {datasets,selfcheck,dev,freeze,test,report}.

Order of operations (PROTOCOL.md section 4):
    selfcheck  ->  dev DATASET  ->  freeze DATASET  ->  commit + push  ->  test DATASET  ->  report
"""

from __future__ import annotations

import argparse
import json
import sys

from . import config as C


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m realbench")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("datasets", help="list datasets, raw files found, eligible curves and clusters per split")
    s = sub.add_parser("selfcheck", help="protocol section 6 checks (synthetic + dev data only)")
    s.add_argument("--quick", action="store_true", help="skip the slow reproduction of the published benchmark")
    d = sub.add_parser("dev", help="run all methods on the dev split and tune smoothing")
    d.add_argument("dataset")
    d.add_argument("--procs", type=int)
    f = sub.add_parser("freeze", help="pre-register a dataset: hash protocol, code, data, tuned settings")
    f.add_argument("dataset")
    f.add_argument("--amend", metavar="REASON")
    t = sub.add_parser("test", help="run the held-out test split (refuses unless frozen, unchanged and pushed)")
    t.add_argument("dataset")
    t.add_argument("--procs", type=int)
    t.add_argument("--allow-unpushed", action="store_true", help="skip the public-commit check (recorded in the run log)")
    r = sub.add_parser("report", help="write the report for a split")
    r.add_argument("--split", default="test", choices=["dev", "test"])
    a = ap.parse_args(argv)

    if a.cmd == "datasets":
        from .datasets import DATASETS, load_split, raw_files
        for ds, spec in DATASETS.items():
            files = raw_files(ds)
            line = f"{ds:14s} {spec.category:12s} {len(files):5d} raw files"
            if files:
                st = {}
                load_split(ds, "test", cap=0, stats=st)
                line += (f" | parsed {st['raw_curves']}, too short {st['too_short']}, duplicates {st['duplicates']}"
                         f" | eligible dev/test {st['eligible']['dev']}/{st['eligible']['test']}"
                         f" | clusters dev/test {st['clusters']['dev']}/{st['clusters']['test']}")
                if st.get("skips"):
                    line += f" | skips {st['skips']}"
                if ds in C.CONFIRMATORY_DATASETS:
                    line += "  [confirmatory]"
                if st["clusters"]["test"] < C.MIN_TEST_CLUSTERS:
                    line += "  [below power gate]"
            print(line)

    elif a.cmd == "selfcheck":
        from .selfcheck import run
        sys.exit(0 if run(quick=a.quick) else 1)

    elif a.cmd == "dev":
        from .runner import run_dev
        df, meta = run_dev(a.dataset, procs=a.procs)
        print(f"dev: {meta['n_curves']} curves, {len(df)} rows, {meta['seconds']} s; chosen smoothing c = {meta['tuned']}")

    elif a.cmd == "freeze":
        from .lock import freeze
        e = freeze(a.dataset, a.amend)
        print(f"frozen {a.dataset}: confirmatory={e['confirmatory']} code={e['code_hash'][:12]} "
              f"data={e['data_hash'][:12]} tuned={e['tuned']}")
        print("Next: commit realbench/registry.json and push to the public repo BEFORE running the test split.")

    elif a.cmd == "test":
        from .runner import run_test
        df, meta = run_test(a.dataset, procs=a.procs, allow_unpushed=a.allow_unpushed)
        print(f"test: {meta['n_curves']} curves, {len(df)} rows ({meta['run_id']}). Run `python -m realbench report`.")

    elif a.cmd == "report":
        from .report import build_report
        from .runner import RESULTS
        txt = build_report(a.split)
        out = RESULTS / f"report_{a.split}.md"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(txt)
        print(txt)
        print(f"\nwritten to {out}")


if __name__ == "__main__":
    main()
