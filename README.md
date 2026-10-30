# curvature-quantiles

Code, data and paper for **"Curvature Quantiles: Turning as a Distribution Along Planar Curves, with an Application to Fixed-Budget Polyline Simplification"** (M. J. Wish, preprint, October 2026). DOI: [10.5281/zenodo.23483923](https://doi.org/10.5281/zenodo.23483923).

The idea: for a planar curve, `|κ| ds = |dθ|`. Normalize it and you have a probability distribution of *turning* along the curve, with a median, quantiles, quantile gaps and a Wasserstein-1 distance. Placing `K` points at equal quantiles of `|κ|^(1/2) ds` is also a fast, single-pass polyline simplifier.

## Install

```bash
pip install -e .                 # library only (numpy, scipy)
pip install -e ".[experiments,test]"
```

## Use

```python
import numpy as np
from curvquant import profile, curvature_median, profile_distance, quantile_resample

x = np.linspace(0, 2, 2001)
pts = np.column_stack([x, x**2])

curvature_median(pts)                    # ~ (0.3904, 0.1524); exact x = (sqrt(17)-1)/8
prof = profile(pts)                      # pure turning distribution
prof.quantile([0.25, 0.5, 0.75])         # quantiles in normalized arc length
hyb = profile(pts, beta=1.0)             # hybrid: ds/L + beta*|kappa| ds

profile_distance(prof, profile(pts[::-1]))        # W1 between turning distributions

simplified = quantile_resample(pts, K=30, beta=None, p=0.5, pin_corners=True)
```

Baselines in `curvquant.resample`: `uniform_arclength`, `rdp_fixed_k` (Ramer–Douglas–Peucker to exactly K vertices), `vw_fixed_k` (Visvalingam–Whyatt), and `optimal_subset_min_eps` (exact dynamic-programming optimum over vertex subsets, O(N³), for N up to ~1000).

## Reproduce the paper

```bash
python experiments/run_all.py     # ~5 min; writes results/, paper/figures/, paper/generated/
python -m pytest -q
cd paper && latexmk -pdf main.tex
```

Every number in the paper is read from `paper/generated/numbers.tex`, which `run_all.py` writes. Nothing is typed in by hand.

## Real-data benchmark (pre-registered, in progress)

`realbench/` tests the resampling claim on real curves: Natural Earth maps, GPX tracks, pen strokes and CAD outlines. The plan is fixed in [`realbench/PROTOCOL.md`](realbench/PROTOCOL.md) before any held-out result exists. Raw data layout is described in [`data/README.md`](data/README.md).

```bash
pip install -e ".[realdata]"
python -m realbench datasets          # what data is present; eligible curves and clusters per split
python -m realbench selfcheck         # protocol section 6 checks (synthetic + dev data only)
python -m realbench dev ne_rivers     # dev split: all methods, tune smoothing
python -m realbench freeze ne_rivers  # hash protocol, code, data, tuned settings -> realbench/registry.json
git commit -am "freeze ne_rivers" && git push   # public timestamp = the pre-registration
python -m realbench test ne_rivers    # refuses if anything changed since the freeze or it was not pushed
python -m realbench report            # every pre-registered result, wins and losses
```

## What holds up and what doesn't

Verified by this code:
- Closed-form curvature medians on convex graphs; the median at an inflection is square-root sensitive to boundary-slope mismatch.
- Profile distance = Wasserstein-1, exactly similarity- and reparameterization-invariant; on convex arcs it equals the normalized turning-function distance.
- On nine **synthetic**, noise-free curve families of 1000 points each, quantile resampling with p = 1/2 beats fixed-K RDP (geometric-mean error ratio ≈ 0.72 over 45 family/budget pairs) and lands within ≈1.3× of the optimal vertex subset. This is a statement about those synthetic curves only.

Not supported:
- Under measurement noise, smoothed RDP beats quantile resampling in every case tested.
- Any claim about real data. Preliminary development-split runs of the real-data benchmark (`realbench/`, branch `realbench`) show quantile resampling losing to RDP on maps, pen strokes and one GPS run. The pre-registered test is not yet run; see `realbench/PROTOCOL.md` section 8.
- As a stroke classifier the unsigned profile distance is weaker than the signed turning function under strong perturbation.
- The profile distance ignores the sign of curvature (an S and a U can be at distance 0) and depends on the start point for closed curves.

## Prior work

The cumulative distribution here coincides with the curvature-based canonical parameterizations of Tumpach (GSI 2023, doi:10.1007/978-3-031-38271-0_4). The |κ|^(1/2) density is classical (Fejes Tóth 1948; McClure & Vitale 1975). See the paper for the full list.

## Cite

See `CITATION.cff`. Archived on Zenodo: [10.5281/zenodo.23483923](https://doi.org/10.5281/zenodo.23483923) (concept DOI, always the latest version; v0.1.0 is 10.5281/zenodo.23483924).

## License

MIT (code). The paper text is © the author, shared under CC BY 4.0.
