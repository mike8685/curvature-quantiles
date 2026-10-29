# Real-data benchmark: pre-registered protocol

Version 1.0 · drafted 2026-10-28 · status: **DRAFT (not frozen)**

This document fixes, before any held-out result exists, how curvature-quantile
resampling (CQ) is tested against established polyline simplifiers on real
curves. Every number the harness produces follows from this file and
`realbench/config.py`. Anything not specified here is **exploratory** and is
labelled that way wherever it is reported.

## 1. Question

At an equal vertex budget K, does CQ approximate real-world curves (maps, GPS
tracks, pen strokes, CAD outlines) more accurately than Ramer-Douglas-Peucker
(RDP) and Visvalingam-Whyatt (VW)?

The synthetic benchmark in the preprint (9 clean analytic families) found a
geometric-mean error ratio CQ/RDP of 0.72. That is the claim under test. It
may not transfer: real data is noisy, irregularly sampled, often already
simplified, and full of straight runs and corners.

## 2. Methods compared

All methods receive the same polyline and the same budget K, counted with
endpoints included. Output points must lie on the **original** input polyline.
Errors are always measured against the original input.

| id | method | role |
| --- | --- | --- |
| `cq` | CQ exactly as published: density proportional to sqrt(curvature) (p = 1/2, pure measure), corners above pi/6 pinned | method under test |
| `rdp` | RDP driven to exactly K vertices (largest deviation first) | primary baseline |
| `vw` | Visvalingam-Whyatt (effective area) down to K vertices | primary baseline |
| `uniform` | K points at equal arc length | floor |
| `cq_unpinned` | CQ without corner pinning | ablation |
| `cq_snapped` | `cq` (pinned) with each point moved to the nearest input vertex (vertex-subset version, directly comparable to RDP/VW). The preprint's "snapped" row snapped the unpinned variant | ablation |
| `cq_p1` | points at equal quantiles of plain turning (p = 1), the "naive" version | ablation |
| `oracle` | exact min-max optimal vertex subset (only when N <= 800) | descriptive bound |

Implementations are the frozen ones in `src/curvquant/resample.py` (tag
v0.1.0). The harness wraps them and does not modify them.

Comparison class. RDP, VW, `cq_snapped` and `oracle` keep a subset of input
vertices. `cq`, `cq_unpinned`, `cq_p1` and `uniform` may place points anywhere
on the input curve. That freedom is part of the method. `cq_snapped` is reported
to show how much of any advantage comes from it.

### 2.1 Preprocessing (applied equally)

Real data is noisy, and curvature estimates amplify noise, so every
shape-aware method (`cq`, `rdp`, `vw`) may **decide** point positions on a
smoothed copy of the curve:

- Gaussian smoothing in arc length, sigma = c x (median edge length), with
  c in {0, 1, 2, 4, 8, 16}. Endpoints are held fixed. Vertex i of the smoothed
  copy corresponds to vertex i of the original.
- The method runs on the smoothed copy and returns fractional vertex indices.
  The output points are then placed at those indices **on the original
  polyline**.
- c is chosen per method and per dataset on the development split only
  (Section 4), and every method gets the same grid. c = 0 (no smoothing) is
  always in the grid, so tuning can never make a method worse on dev.

No other parameter is tuned. p = 1/2, the pinning threshold pi/6 and every
baseline setting stay fixed at their published values.

## 3. Data

### 3.1 Categories and planned datasets

| dataset id | category | source | status at drafting |
| --- | --- | --- | --- |
| `ne_coastline` | maps | Natural Earth 1:10m coastline | in hand |
| `ne_rivers` | maps | Natural Earth 1:10m river centerlines | in hand |
| `ne_roads` | maps | Natural Earth 1:10m roads | in hand |
| `gpx` | GPS | GPX tracks (gpxpy test files; more to be added) | 8 tracks, too few |
| `quickdraw_raw` | pen strokes | Google Quick Draw, **raw** ndjson only | not yet obtained |
| `inkml` | pen strokes | CROHME / MathWriting InkML | 1 sample file |
| `uji` | pen strokes | UJI Pen Characters v2 | not yet obtained |
| `fusion360` | CAD | Fusion 360 Gallery reconstruction JSON (profile loops) | parser checked on line-only samples |
| `dxf` | CAD | any DXF corpus (fallback) | not yet obtained |

Six datasets are confirmatory: `ne_coastline`, `ne_rivers`, `ne_roads`,
`gpx`, one pen-stroke dataset and one CAD dataset. Which pen and CAD dataset
counts is fixed at freeze time, before their test splits are run. Others are
exploratory.

### 3.2 Unit of analysis and clusters

A **curve** is one polyline: a coastline linestring, one GPX track segment,
one pen stroke, one closed CAD profile loop. Curves from the same source are
not independent (strokes from one drawing, segments from one GPX file, parts of
one river), so each curve carries a **cluster** id: the GeoJSON feature, the
GPX file, the drawing / InkML file / writer, the CAD design. Statistics treat
clusters as the independent units (Section 5).

### 3.3 Inclusion (fixed, blind to method performance)

- After consecutive duplicate points are removed, N >= 40 vertices.
- Finite coordinates, positive length.
- Geographic data is projected to a local equirectangular plane about each
  curve's centroid, in km for maps and m for GPS.
- A curve is closed when its first and last points coincide. Closed status is
  passed to the methods that support it (`uniform`, CQ variants). RDP and VW
  treat the ring as a path from its stored start vertex, as in the preprint.
- CAD loops are sampled at uniform arc-length spacing (loop length / 1000), with
  every primitive endpoint inserted. Curvature-dependent sampling would hand
  vertex-subset methods free corner vertices and CQ a pre-made density, so it is
  not used.
- **Quick Draw `simplified` files are excluded**: Google produced them with RDP,
  which would bias the test toward RDP. Raw files only.

There is no filter on curvature, shape or "difficulty". Curves on which a
method fails stay in (Section 5.4).

### 3.4 Splits

Each **cluster** is assigned to dev or test by
`sha256("realbench-v1:" + dataset + ":" + cluster) / 2^256 < 0.30`, which gives
roughly 30% dev and 70% test. Assignment depends only on the ids, never on
results. Adding data later never moves an existing cluster.

Exact duplicate curves (identical after normalizing position and scale) are
removed across the whole dataset before splitting, keeping the first in file
order, so the same curve can never sit in both dev and test.

If a split holds more than 400 curves (test) or 150 curves (dev), a simple
random sample is drawn **before** any method runs: the curves with the smallest
`sha256(seed 20261028, dataset, curve id)`. This is a uniform random sample that
can be drawn while streaming large files.

## 4. Procedure and order of operations

1. Commit this protocol and the harness code.
2. **Dev phase (per dataset):** run every method across the smoothing grid on
   dev and pick c per method (Section 2.1). Dev results may be inspected freely
   and used to find bugs. They are never reported as evidence.
3. **Freeze:** `python -m realbench freeze` records SHA-256 hashes of this file,
   `config.py`, all harness and `curvquant` source, the tuned c values and every
   raw data file used. Then commit and **push the freeze commit to the public
   repo before step 4.** The public commit timestamp is the pre-registration.
4. **Test phase:** `python -m realbench run DATASET --split test`. The harness
   refuses to run if any frozen hash has changed. Every test run is appended to
   `results/realbench/runlog.jsonl`, including failed or repeated runs.
5. **Report:** `python -m realbench report` writes every confirmatory and
   secondary result, wins and losses alike.

A dataset's test split is analysed once. If a bug in the harness (not an
unwelcome result) forces a rerun, both runs are reported, with the bug
described. If more data arrives for a dataset that has already been tested, it
becomes a new dataset version, reported separately as a replication, never
merged into the original.

## 5. Analysis

### 5.1 Budgets

For each curve with N vertices, K = max(8, round(r x N)) for
r in {0.02, 0.05, 0.10, 0.20}. Any K > N/2 is dropped, as are duplicate K
values on that curve. A curve with no remaining K is excluded and counted.

### 5.2 Error metrics

The reference is the original polyline, evaluated at all of its vertices plus
max(4N, 2000) points at uniform arc length (capped at 200,000), so metrics are
not dominated by dense vertex clusters such as a GPS unit idling at a stop.

- **Primary: directed Hausdorff error** `H`, the largest distance from the
  reference to the simplified polyline (nearest segment). This is the
  preprint's metric, extended from vertices to vertices plus arc-length
  samples.
- Secondary: arc-length RMS distance to the simplified polyline; symmetric
  Hausdorff; "shortcut" error (distance to the simplified segment spanning the
  same stretch of arc length, the objective the oracle minimizes); runtime.

### 5.3 Confirmatory comparisons (per confirmatory dataset)

| id | comparison | smoothing |
| --- | --- | --- |
| H1 | `cq` vs `rdp` | c = 0 for both (the published methods as-is) |
| H2 | `cq` vs `vw` | c = 0 for both |
| H3 | `cq` vs `rdp` | each at its own dev-tuned c |

For each curve and K: `lr = log(H_A / H_B)`. When both errors are 0, lr = 0.
When only one is, the zero is floored at 1e-12 x curve diameter. Curve score =
mean lr over that curve's K values. Cluster score = mean over the cluster's
curves.

- **Test:** two-sided Wilcoxon signed-rank on cluster scores (zero
  differences dropped, as in Wilcoxon's original procedure).
- **Effect size:** geometric-mean ratio `GMR = exp(mean cluster score)`, with a
  95% percentile bootstrap CI over clusters (10,000 resamples, seed 20261028).
- **Win rate:** share of clusters where A beats B.
- **Multiplicity:** 6 datasets x 3 comparisons = 18 planned tests. Each test
  passes at Bonferroni alpha = 0.05 / 18 = 0.00278, applied as each dataset
  completes. The final report also gives Holm-adjusted p-values. Tests that are
  never run stay in the denominator.

**Power gate:** a dataset with fewer than 30 test clusters is reported as
exploratory, whatever its p-value.

### 5.4 Decision rules (fixed now)

For each confirmatory comparison:

- **CQ better:** p < alpha and the CI upper bound < 1.
- **CQ worse:** p < alpha and the CI lower bound > 1.
- **Practically equivalent:** the CI lies entirely inside [1/1.05, 1.05].
- **Inconclusive:** anything else.

Failures: if a method raises an error on a curve, its output is replaced by
`uniform` for that curve and K (intention-to-treat), and failures per method
are reported. A method returning fewer than K points (`cq_snapped`, rarely
RDP/VW) is scored as returned, and the shortfall is reported.

### 5.5 Secondary analyses (reported in full, no confirmatory claims)

- H1 to H3 under each secondary metric.
- Breakdown by actual budget fraction K/N (bins at 3%, 7%, 15%; the K >= 8
  floor makes the nominal r misleading on short curves), and Spearman
  correlations of per-curve H1 log ratio with N, total turning, straight-run
  fraction and K/N.
- Ablations against `rdp`: `cq_unpinned`, `cq_snapped`, `cq_p1`, `uniform`.
- Ratio to `oracle` (shortcut error, N <= 800).
- The chosen c per method, and the dev-to-test change.
- Runtime (Python; indicative only).

## 6. Checks the harness must pass before any test run

`python -m realbench selfcheck` runs:

1. **Reproduction:** on the preprint's synthetic families, the harness, with
   vertex-only evaluation, reproduces the published per-curve errors exactly.
2. **Oracle bound:** under the shortcut metric, `oracle` <= every vertex-subset
   method on every curve and K. A violation means a bug.
3. **Identity and symmetry:** A vs A gives GMR = 1; swapping A and B inverts
   GMR.
4. **Invariance:** rotating, translating and uniformly scaling a curve leaves
   every error ratio unchanged.
5. **On-curve:** every output point lies on the original polyline (tolerance
   1e-9 x diameter).
6. **Smoothing mapping:** with c = 0, the smoothed path reproduces the
   unsmoothed outputs exactly.

## 7. Known threats to validity (stated in advance)

- Natural Earth geometry was already generalized by its makers. Pre-simplified
  vertices may favour vertex-subset methods. Raw survey data would be a
  stronger test.
- The GPS set is tiny until more tracks are added, so it stays exploratory under
  the power gate.
- The primary metric is the objective RDP greedily targets. RMS is reported
  because an equal-error density such as CQ's may do better there. The primary
  metric was chosen to test the preprint's claim, not to favour either side.
- The method's author (M. J. Wish) runs the study. The safeguards are a public
  freeze before test, an append-only run log, and an automated report that
  prints every pre-registered result.
- Timings come from Python implementations of unequal optimization and are
  indicative only.

## 8. Changes and disclosures

Any change after freeze is recorded here with its date and reason. A change
that touches a frozen hash re-locks the protocol, and any dataset already tested
under the old lock is reported under both versions.

- 2026-10-28 (before any freeze), disclosure: while the harness was being
  built, the `ne_roads` **dev** split was run once as a smoke test. On dev, CQ
  had higher error than RDP (GMR about 1.7). After that, only reporting was
  changed: the per-K table was re-binned by actual K/N, and the
  curve-property correlations already listed in 5.5 were implemented. No
  hypothesis, metric, budget, method, inclusion rule or decision rule was
  changed after any result had been seen.
