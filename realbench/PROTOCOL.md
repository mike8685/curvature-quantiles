# Real-data benchmark: pre-registered protocol

Version 1.1 · 2026-10-29 · status: **DRAFT (not frozen)**

This document fixes, before any held-out result exists, how curvature-quantile
resampling (CQ) is tested against established polyline simplifiers on real
curves. Every number the harness produces follows from this file and
`realbench/config.py`. Anything not specified here is **exploratory** and is
labelled that way wherever it is reported. Version 1.0 was revised after an
independent adversarial review; Section 8 lists every change and every result
that was seen before this version.

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
endpoints included. On a closed curve every method returns K points whose
first and last points coincide. Output points must lie on the **original** input
polyline, and errors are always measured against the original input.

| id | method | role |
| --- | --- | --- |
| `cq` | CQ exactly as published: density proportional to sqrt(curvature) (p = 1/2, pure measure), corners above pi/6 pinned | method under test |
| `rdp` | RDP driven to exactly K vertices (largest deviation first) | primary baseline |
| `vw` | Visvalingam-Whyatt (effective area) down to K vertices | primary baseline |
| `uniform` | K points at equal arc length | floor |
| `cq_unpinned` | CQ without corner pinning | ablation |
| `cq_snapped` | `cq` (pinned) with each point moved to the nearest input vertex. This is the vertex-subset version, directly comparable to RDP/VW. The preprint's "snapped" row snapped the unpinned variant | ablation |
| `cq_p1` | the preprint's naive row: equal quantiles of turning (p = 1) blended with arc length (beta = 1) | ablation |
| `oracle` | exact min-max optimal vertex subset (only when N <= 800) | descriptive bound |

Implementations are the frozen ones in `src/curvquant/resample.py` (tag
v0.1.0). The harness wraps them and does not modify them.

Comparison class: RDP, VW, `cq_snapped` and `oracle` keep a subset of input
vertices. `cq`, `cq_unpinned`, `cq_p1` and `uniform` may place points anywhere
on the input curve; that freedom is part of the method. `cq_snapped` is reported
to show how much of any advantage comes from it.

### 2.1 Preprocessing (applied equally)

Real data is noisy and curvature estimates amplify noise, so every shape-aware
method (`cq`, `rdp`, `vw`) may **decide** point positions on a smoothed copy of
the curve:

- **Smoothing:** Gaussian smoothing in arc length, with sigma = c x (median edge
  length) and c in {0, 1, 2, 4, 8, 16}.
- **Boundaries do not add curvature.**
  - Closed curves are smoothed periodically.
  - Open curves are padded by mirroring across the line through each endpoint
    perpendicular to its tangent (a second-order one-sided estimate), so
    curvature continues evenly past the ends.
  - Vertex i of the smoothed copy corresponds to vertex i of the original.
- **Placement on the original:** the method runs on the smoothed copy and
  returns fractional vertex indices. The output points are then placed at
  those indices **on the original polyline**.
- **Tuning:** c is chosen per method and per dataset on the development split
  only (Section 4). Every method gets the same grid. c = 0 (no smoothing) is
  always in the grid.

No other parameter is tuned. p = 1/2, the pinning threshold pi/6 and every
baseline setting stay fixed at their published values.

## 3. Data

### 3.1 Datasets (all six confirmatory datasets are named now)

| dataset id | category | source | role |
| --- | --- | --- | --- |
| `ne_coastline` | maps | Natural Earth 1:10m coastline | confirmatory |
| `ne_rivers` | maps | Natural Earth 1:10m river centerlines | confirmatory |
| `ne_roads` | maps | Natural Earth 1:10m roads | confirmatory |
| `geolife` | GPS | Microsoft GeoLife Trajectories 1.3, every trajectory of every user | confirmatory |
| `quickdraw_raw` | pen strokes | Google Quick Draw **raw** ndjson, 10 categories fixed by hash (below) | confirmatory |
| `fusion360` | CAD | Fusion 360 Gallery reconstruction dataset, all designs, profile loops | confirmatory |
| `gpx` | GPS | gpxpy's bundled test tracks (parser fixtures, not a sample) | exploratory |
| `inkml` | pen strokes | CROHME / MathWriting InkML | exploratory |
| `uji` | pen strokes | UJI Pen Characters v2 | exploratory |
| `dxf` | CAD | any DXF corpus | exploratory |

**Quick Draw categories.** These are the 10 of the 345 official categories with
the smallest `sha256("realbench-v1:quickdraw:" + name)`: moustache, sun,
binoculars, stop sign, paint can, strawberry, crocodile, bandage, ambulance,
nose. They were chosen by hash, not by eye, and before any pen data was seen.
The loader reads no other category.

**Whole corpora only.** GeoLife, Quick Draw and Fusion 360 are used in full, as
published, with no hand selection. A confirmatory dataset whose data never
arrives stays in the multiplicity denominator and is reported as not run.

### 3.2 Unit of analysis and clusters

A **curve** is one polyline: a coastline linestring, one GeoLife trajectory, one
pen stroke, one closed CAD profile loop. Curves that share a source are not
independent, so each curve carries a **cluster**, and statistics treat clusters
as the independent units (Section 5):

| data | cluster |
| --- | --- |
| maps | Features with the same name (rivers: `name`; roads: country + `name`) form one group. Each group goes to the 10 x 10 degree cell of its centroid, and the cell is the cluster. Nearby and same-named features always share a cluster. |
| GeoLife | the user |
| Quick Draw | the drawing |
| Fusion 360 | the project |
| InkML | the writer annotation, else the file |
| UJI | the writer |
| GPX / DXF | the file |

### 3.3 Inclusion (fixed, blind to method performance)

- After consecutive duplicate points are removed, N >= 40 vertices.
- Finite coordinates and positive length.
- Geographic data is projected to a local equirectangular plane about each
  curve's centroid, in km for maps and m for GPS.
- A curve is closed when its ends are within 1e-6 x diameter, and it is then
  snapped exactly shut. Closed status is passed to the methods that support it
  (`uniform`, CQ variants) and to smoothing. RDP and VW treat a ring as a path
  from its stored start vertex, as in the preprint.
- Exact duplicate curves are removed across the whole dataset before splitting,
  whatever their position, scale, direction or start point. The first in file
  order is kept.
- GPX: only track segments are used; routes and waypoints are plans, not
  recordings.
- CAD loops are sampled at uniform arc-length spacing (loop length / 1000) with
  every primitive endpoint inserted. Curvature-dependent sampling is not used,
  because it would hand vertex-subset methods free corner vertices and CQ a
  ready-made density.
- Fusion 360 loops containing curve types the loader cannot sample exactly
  (splines, ellipses), loops that do not close, and circles inside multi-curve
  loops are skipped. Skips are counted by reason and reported. This makes the
  CAD set line-and-arc loops only (Section 7).
- Quick Draw `simplified` files are refused: Google made them with RDP.

There is no filter on curvature, shape or difficulty. Curves on which a method
fails stay in (Section 5.4). Everything a loader drops is counted and reported.

### 3.4 Splits

Each **cluster** is assigned to dev or test by
`sha256("realbench-v1:" + dataset + ":" + cluster) / 2^256 < 0.30`, which puts
about 30% in dev and 70% in test. Assignment depends only on the ids. Adding data
later never moves an existing cluster.

If a split holds more than 400 curves (test) or 150 curves (dev), a simple
random sample is drawn **before** any method runs: the curves with the smallest
`sha256(seed 20261028, dataset, curve id)`.

## 4. Procedure and order of operations

1. Commit this protocol and the harness code.
2. `python -m realbench selfcheck` must pass (Section 6).
3. **Dev phase (per dataset):** `python -m realbench dev DATASET` runs every
   method across the smoothing grid on dev and picks c per method. Dev results
   may be inspected freely and used to find bugs. They are never evidence. Any
   dev result seen before the freeze is disclosed in Section 8.
4. **Freeze:** `python -m realbench freeze DATASET` records SHA-256 hashes of
   this file, `config.py`, all harness and `curvquant` source, the package
   versions, the dev rows and tuned c values, and every raw data file. It
   refuses if the dev run used different code. Then commit and **push the
   freeze commit to the public GitHub repo before step 5.** The public commit
   timestamp is the pre-registration.
5. **Test phase:** `python -m realbench test DATASET`. The harness refuses if any
   frozen hash or package version has changed, or if HEAD is not the tip of a
   branch on the real remote (checked with `git ls-remote`). It appends a start
   record to `results/realbench/runlog.jsonl`, then commits and pushes that
   record **before** computing anything. Output goes to a fresh `run-<UTC
   time>` directory, and the hash of the rows is logged when the run finishes.
   `--allow-unpushed` skips the public steps; that is logged and printed in the
   report.
6. **Report:** `python -m realbench report` uses, for each dataset, the
   **first** logged and completed test run whose rows match their logged hash.
   It lists every other run, and refuses outputs that are absent from the log or
   altered. Confirmatory status is taken from the start record.

A dataset's test split is analysed once. If a harness bug (not an unwelcome
result) forces a rerun, both runs are reported, with the bug described. More data
for an already-tested dataset becomes a new dataset version, reported separately
as a replication.

## 5. Analysis

### 5.1 Budgets

For each curve with N vertices, K = max(8, round(r x N)) for
r in {0.02, 0.05, 0.10, 0.20}. Any K > N/2 is dropped, as are duplicate K values
on that curve.

### 5.2 Error metrics

The reference is the original polyline, evaluated at all of its vertices plus
max(4N, 2000) points at uniform arc length (capped at 200,000).

- **Primary: directed Hausdorff error** `H`, the largest distance from the
  reference to the simplified polyline. This is the preprint's metric, extended
  from vertices to vertices plus arc-length samples.
- **Secondary:**
  - arc-length RMS distance;
  - "shortcut" error: the distance to the simplified segment spanning the same
    stretch of arc length, which is the objective the oracle minimizes.
- Symmetric Hausdorff equals H whenever outputs lie on the input, so it is not
  reported.
- Runtime is recorded.

**Exact fits.** Any error below tau = 1e-9 x curve diameter is an exact fit.
Before ratios are taken, errors are raised to tau, so two exact fits tie
(lr = 0) instead of yielding a ratio of rounding noise. Exact-fit counts are
reported for every comparison.

### 5.3 Confirmatory comparisons (per confirmatory dataset)

| id | comparison | smoothing |
| --- | --- | --- |
| H1 | `cq` vs `rdp` | c = 0 for both (the published methods as-is) |
| H2 | `cq` vs `vw` | c = 0 for both |
| H3 | `cq` vs `rdp` | each at its own dev-tuned c |

**Scores.** For each curve and K, `lr = log(max(H_A, tau) / max(H_B, tau))`.
The curve score is the mean lr over that curve's K values, and the cluster score
is the mean over the cluster's curves.

- **Test:** two-sided Wilcoxon signed-rank on cluster scores. Zero scores (ties)
  are dropped, as in Wilcoxon's procedure, and the number of tied clusters is
  reported.
- **Effect size:** geometric-mean ratio `GMR = exp(mean cluster score)`, with
  percentile bootstrap intervals over clusters (20,000 resamples, seed
  20261028) at 95% for display, 1 - alpha for decisions, and 1 - 2 alpha for
  equivalence.
- **Win and loss rates:** the share of clusters where A beats B, and where it
  loses.
- **Multiplicity:**
  - 6 datasets x 3 comparisons = 18 planned tests.
  - Bonferroni alpha = 0.05 / 18 = 0.00278 per test, applied as each dataset
    completes.
  - The final report also gives Holm-adjusted p-values, with tests not run
    counted as p = 1.

**Power gate:** a confirmatory dataset with fewer than 30 test clusters is
reported as exploratory, whatever its p-value.

### 5.4 Decision rules

For each confirmatory comparison:

| decision | condition |
| --- | --- |
| CQ better | p < alpha and the upper bound of the 1 - alpha CI < 1 |
| CQ worse | p < alpha and the lower bound of the 1 - alpha CI > 1 |
| practically equivalent | the 1 - 2 alpha CI lies inside [1/1.05, 1.05] (two one-sided tests at alpha) |
| inconclusive | anything else |

**Failures:** any exception, in a method or in recovering its output positions,
replaces that output with `uniform` for that curve and K (intention-to-treat).
Failure rates and reasons are reported. A method that returns fewer than K
points is scored as returned, and the shortfall is reported. How often corner
pinning switches off (no corner, or more than K-2 corners) is reported.

### 5.5 Secondary analyses (reported in full, no confirmatory claims)

- H1 to H3 under each secondary metric.
- Breakdown by actual budget fraction K/N, in bins at 3%, 7% and 15%.
- Spearman correlations of per-curve H1 log ratio with N, total turning,
  straight-run fraction and K/N.
- Ablations against `rdp`: `vw`, `uniform`, `cq_unpinned`, `cq_snapped`,
  `cq_p1`.
- Ratio to `oracle` (shortcut error, N <= 800).
- The chosen c per method.
- Runtime (Python; indicative only).
- All exploratory datasets.

## 6. Checks the harness must pass before any test run

`python -m realbench selfcheck` runs these checks. No method is ever run on test
curves. Test geometry is parsed only to count clusters.

1. **Reproduction:** with vertex-only evaluation, the harness reproduces the
   preprint's published synthetic-benchmark errors to 1e-9.
2. **Oracle bound:** under the shortcut metric, `oracle` <= every vertex-subset
   method on every curve and K.
3. **Identity and symmetry:** A vs A gives GMR = 1; swapping A and B inverts the
   GMR and keeps p.
4. **Invariance:** rotating, translating and scaling a curve scales every error
   by the scale factor. Rare tie-breaks flipped by floating-point rounding are
   counted and capped at 2%.
5. **On-curve:** every output point lies on the original polyline (tolerance
   1e-8 x diameter), at c = 0 and c = 4.
6. **c = 0 path:** with c = 0, outputs equal the direct method calls exactly.
7. **Smoothing:** it reduces noise, adds no turning at a closed seam or at open
   ends (max/median turning within 10%), and is the identity at c = 0.
8. **Exact-fit ties:** polygons that RDP and VW both fit exactly give GMR = 1.
9. **Split integrity:** no cluster sits in both dev and test.

## 7. Known threats to validity (stated in advance)

- **Pre-generalized maps.** Natural Earth geometry was already generalized by
  its makers. Coarse polylines with sparse corners favour vertex-subset methods,
  and the dev results (Section 8) are consistent with that.
- **Line-and-arc CAD only.** The Fusion 360 set excludes loops with splines and
  ellipses. Such loops are fitted exactly by vertex subsets at small K, which
  favours RDP/VW. Many exact-fit ties are expected, and they are reported.
- **Metric choice.** The primary metric is the objective RDP greedily targets.
  RMS is reported because an equal-error density such as CQ's may do better
  there. The primary metric was chosen to test the preprint's claim.
- **Projection.** Local equirectangular projection distorts the few map curves
  that span tens of degrees. The distortion is identical for all methods.
- **Repeated coordinates.** On tracks that retrace themselves exactly, recovering
  a vertex index can pick an earlier identical vertex. This affects only the
  shortcut metric (secondary).
- **Author conflict.** The method's author (M. J. Wish) runs the study. The
  safeguards are named confirmatory datasets, a public freeze, a pushed
  start record, an append-only run log, and an automated report that prints
  every pre-registered result.
- **Timings** come from Python implementations of unequal optimization and are
  indicative only.

## 8. Changes and disclosures

Any change after freeze is recorded here with its date and reason. A change that
touches a frozen hash re-locks the protocol, and any dataset already tested
under the old lock is reported under both versions.

**Results seen before this version (all dev split; no test split has ever been
run):**

- 2026-10-28 22:04 UTC-4, `ne_roads` dev, protocol v1.0 draft (written but not
  yet committed). H1 GMR 1.73 [1.62, 1.90] (CQ worse than RDP). H2 GMR 1.45
  (CQ worse than VW). Tuned H3 GMR 1.46. Chosen c: cq 2, rdp 0, vw 0.
- 2026-10-28, `ne_rivers` dev, protocol v1.0. H1 GMR 1.88 [1.81, 1.95] (CQ
  worse). H2 GMR 1.46. Tuned H3 GMR 1.51. Chosen c: cq 4, rdp 0, vw 1.

Those runs used v1.0's per-feature clusters, edge-pinned smoothing and
zero-only error floor. They are superseded by the v1.1 dev runs below, which are
committed alongside this version.

- 2026-10-29, v1.1 dev runs, 150 curves each, GMR CQ/RDP. In all three, CQ was
  worse than RDP on every comparison.

  | dataset | H1 | H2 (vs VW) | tuned H3 | chosen c (cq / rdp / vw) | corner pinning switched off |
  | --- | --- | --- | --- | --- | --- |
  | `ne_coastline` | 1.80 | 1.12 | 1.66 | 4 / 0 / 0 | 84% of cq runs |
  | `ne_rivers` | 1.94 | 1.45 | 1.55 | 4 / 0 / 1 | 70% of cq runs |
  | `ne_roads` | 1.74 | 1.52 | 1.39 | 2 / 0 / 0 | 60% of cq runs |

  Pinning switches off mostly because these coarse polylines have more than K-2
  corners above pi/6. The method is **not** being changed in response. Any
  modified CQ variant would be a new method, needing new data or a new
  pre-registration.

**What changed after those results were seen (v1.0 to v1.1).** Every change
below came from an independent adversarial review of the harness. None is a
change to the hypotheses, the primary metric, the budgets or the methods.

1. **Error floor applied to all near-zero errors**, so exact fits tie. Under v1.0,
   rounding noise produced huge spurious ratios, mostly against CQ.
2. **All six confirmatory datasets named in advance:** the GeoLife, Quick Draw
   (with hash-chosen categories) and Fusion 360 corpora. v1.0 allowed the
   pen-stroke and CAD datasets to be picked at freeze time, after dev results.
   The gpxpy fixtures were demoted to exploratory: they were parser test files,
   not a sample.
3. **Map clusters** are now name groups within 10-degree cells. Single features
   were not independent.
4. **Smoothing boundaries** are periodic or mirrored. v1.0 created artificial
   curvature at ends and seams, which penalized CQ.
5. **Test-run protection:**
   - the test split runs only through the verified path, with a pushed start
     record;
   - each run gets its own directory;
   - outputs are hash-logged;
   - the report reads only logged runs;
   - confirmatory status is taken from the start record;
   - package versions are frozen.
6. **Decision intervals** are now at the Bonferroni level (TOST for
   equivalence), not at 95%.
7. **Data handling:**
   - `cq_p1` now matches the preprint's naive row;
   - symmetric Hausdorff is dropped (identical to H);
   - dedupe is invariant to direction and start point;
   - failures in recovering positions are caught;
   - UJI ids are unique per stroke, InkML is clustered by writer, GPX routes are
     excluded;
   - loader skips and pinning switch-off are counted.
