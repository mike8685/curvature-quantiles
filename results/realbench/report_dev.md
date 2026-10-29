# Real-data benchmark report (dev split)

Protocol 1.1. Primary metric: directed Hausdorff error. GMR < 1 means the first method has lower error. Per-test alpha (Bonferroni over 18 planned tests) = 0.00278; better/worse decisions use the Bonferroni-level CI, equivalence uses the 1 - 2 alpha CI within [1/1.05, 1.05]. Errors below 1e-09 x diameter count as exact fits and tie.

**DEV SPLIT. Used only for tuning and debugging. Not evidence.**

## ne_coastline (exploratory)

150 curves; eligible dev/test = {'dev': 308, 'test': 661}; clusters dev/test = {'dev': 61, 'test': 154}; too short 3164, duplicates 0; loader skips none; tuned c = {'cq': 4, 'rdp': 0, 'vw': 0}

| test | comparison | GMR | 95% CI | Bonferroni CI | p | win / loss | clusters / curves | exact fits both / A only / B only | Holm p | decision |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H1 | cq vs rdp (c0) | 1.804 | [1.672, 1.958] | [1.612, 2.051] | 2.4e-10 | 0% / 100% | 53 / 150 (0 tied) | 0 / 0 / 0 |  | (CQ worse; exploratory) |
| H2 | cq vs vw (c0) | 1.124 | [1.020, 1.242] | [0.961, 1.310] | 0.0065 | 32% / 68% | 53 / 150 (0 tied) | 0 / 0 / 0 |  | (inconclusive; exploratory) |
| H3 | cq vs rdp (tuned) | 1.663 | [1.556, 1.783] | [1.508, 1.857] | 2.8e-10 | 2% / 98% | 53 / 150 (0 tied) | 0 / 0 / 0 |  | (CQ worse; exploratory) |

Secondary metrics (same comparisons):

| comparison | GMR | 95% CI | Bonferroni CI | p | win / loss | clusters / curves | exact fits both / A only / B only |
| --- | --- | --- | --- | --- | --- | --- | --- |
| H1/rmse | 1.721 | [1.634, 1.830] | [1.593, 1.895] | 2.4e-10 | 0% / 100% | 53 / 150 (0 tied) | 0 / 0 / 0 |
| H1/shortcut | 1.815 | [1.684, 1.969] | [1.625, 2.059] | 2.4e-10 | 0% / 100% | 53 / 150 (0 tied) | 0 / 0 / 0 |
| H2/rmse | 1.405 | [1.303, 1.521] | [1.254, 1.588] | 9.9e-09 | 6% / 94% | 53 / 150 (0 tied) | 0 / 0 / 0 |
| H2/shortcut | 1.103 | [0.995, 1.225] | [0.939, 1.295] | 0.027 | 36% / 64% | 53 / 150 (0 tied) | 0 / 0 / 0 |
| H3/rmse | 1.518 | [1.440, 1.602] | [1.405, 1.651] | 2.5e-10 | 2% / 98% | 53 / 150 (0 tied) | 0 / 0 / 0 |
| H3/shortcut | 1.688 | [1.585, 1.803] | [1.537, 1.875] | 2.5e-10 | 2% / 98% | 53 / 150 (0 tied) | 0 / 0 / 0 |

Ablations and baselines vs RDP (c = 0):

| comparison | GMR | 95% CI | Bonferroni CI | p | win / loss | clusters / curves | exact fits both / A only / B only |
| --- | --- | --- | --- | --- | --- | --- | --- |
| vw | 1.604 | [1.455, 1.798] | [1.391, 1.917] | 3e-10 | 2% / 98% | 53 / 150 (0 tied) | 0 / 0 / 0 |
| uniform | 1.793 | [1.685, 1.910] | [1.637, 1.974] | 2.4e-10 | 0% / 100% | 53 / 150 (0 tied) | 0 / 0 / 0 |
| cq_unpinned | 1.747 | [1.639, 1.866] | [1.588, 1.934] | 2.4e-10 | 0% / 100% | 53 / 150 (0 tied) | 0 / 0 / 0 |
| cq_snapped | 1.836 | [1.706, 1.989] | [1.651, 2.079] | 2.4e-10 | 0% / 100% | 53 / 150 (0 tied) | 0 / 0 / 0 |
| cq_p1 | 2.173 | [1.986, 2.402] | [1.909, 2.532] | 2.4e-10 | 0% / 100% | 53 / 150 (0 tied) | 0 / 0 / 0 |

Ratio to optimal vertex subset (shortcut error, N <= 800):

| comparison | GMR | 95% CI | Bonferroni CI | p | win / loss | clusters / curves | exact fits both / A only / B only |
| --- | --- | --- | --- | --- | --- | --- | --- |
| cq | 2.438 | [2.262, 2.639] | [2.186, 2.768] | 3.6e-15 | 0% / 100% | 49 / 141 (0 tied) | 0 / 0 / 0 |
| rdp | 1.401 | [1.362, 1.443] | [1.343, 1.466] | 3.6e-15 | 0% / 100% | 49 / 141 (0 tied) | 0 / 0 / 0 |
| vw | 2.031 | [1.902, 2.172] | [1.838, 2.253] | 3.6e-15 | 0% / 100% | 49 / 141 (0 tied) | 0 / 0 / 0 |
| cq_snapped | 2.485 | [2.323, 2.675] | [2.253, 2.795] | 3.6e-15 | 0% / 100% | 49 / 141 (0 tied) | 0 / 0 / 0 |

GMR by actual budget fraction K/N (clusters): H1: <=3% 1.925 (16), 3-7% 2.222 (26), 7-15% 1.872 (46), >15% 1.877 (53); H2: <=3% 0.923 (16), 3-7% 1.013 (26), 7-15% 1.147 (46), >15% 1.077 (53); H3: <=3% 2.247 (16), 3-7% 2.004 (26), 7-15% 1.538 (46), >15% 1.774 (53)

Spearman correlation of H1 per-curve log ratio with curve properties (positive = CQ relatively worse as the property grows): N +0.41, turning +0.28, straight +0.26, kn -0.42

Failure rates (fallback to uniform): none

Mean point shortfall below K: none

Share of cq runs (c = 0) where corner pinning switched off (no corner, or more than K-2): 97.6%

Median runtime ms (Python, indicative): cq@c0 0.36, cq@c1 0.34, cq@c2 0.44, cq@c4 0.44, cq@c8 0.36, cq@c16 0.35, cq_p1@c0 0.27, cq_snapped@c0 0.40, cq_unpinned@c0 0.34, oracle@c0 9.01, rdp@c0 0.28, rdp@c1 0.29, rdp@c2 0.27, rdp@c4 0.26, rdp@c8 0.28, rdp@c16 0.27, uniform@c0 0.10, vw@c0 0.82, vw@c1 0.80, vw@c2 0.79, vw@c4 0.79, vw@c8 0.75, vw@c16 0.79

## ne_rivers (exploratory)

150 curves; eligible dev/test = {'dev': 441, 'test': 967}; clusters dev/test = {'dev': 49, 'test': 132}; too short 1034, duplicates 0; loader skips none; tuned c = {'cq': 4, 'rdp': 0, 'vw': 1}

| test | comparison | GMR | 95% CI | Bonferroni CI | p | win / loss | clusters / curves | exact fits both / A only / B only | Holm p | decision |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H1 | cq vs rdp (c0) | 1.936 | [1.846, 2.034] | [1.809, 2.089] | 9.1e-13 | 0% / 100% | 41 / 150 (0 tied) | 0 / 0 / 0 |  | (CQ worse; exploratory) |
| H2 | cq vs vw (c0) | 1.447 | [1.377, 1.520] | [1.343, 1.560] | 2.7e-12 | 2% / 98% | 41 / 150 (0 tied) | 0 / 0 / 0 |  | (CQ worse; exploratory) |
| H3 | cq vs rdp (tuned) | 1.548 | [1.493, 1.606] | [1.464, 1.640] | 9.1e-13 | 0% / 100% | 41 / 150 (0 tied) | 0 / 0 / 0 |  | (CQ worse; exploratory) |

Secondary metrics (same comparisons):

| comparison | GMR | 95% CI | Bonferroni CI | p | win / loss | clusters / curves | exact fits both / A only / B only |
| --- | --- | --- | --- | --- | --- | --- | --- |
| H1/rmse | 1.932 | [1.843, 2.028] | [1.806, 2.077] | 9.1e-13 | 0% / 100% | 41 / 150 (0 tied) | 0 / 0 / 0 |
| H1/shortcut | 1.938 | [1.849, 2.035] | [1.811, 2.089] | 9.1e-13 | 0% / 100% | 41 / 150 (0 tied) | 0 / 0 / 0 |
| H2/rmse | 1.747 | [1.670, 1.829] | [1.638, 1.869] | 9.1e-13 | 0% / 100% | 41 / 150 (0 tied) | 0 / 0 / 0 |
| H2/shortcut | 1.442 | [1.370, 1.515] | [1.336, 1.556] | 2.7e-12 | 2% / 98% | 41 / 150 (0 tied) | 0 / 0 / 0 |
| H3/rmse | 1.487 | [1.441, 1.533] | [1.415, 1.558] | 9.1e-13 | 0% / 100% | 41 / 150 (0 tied) | 0 / 0 / 0 |
| H3/shortcut | 1.551 | [1.497, 1.608] | [1.469, 1.641] | 9.1e-13 | 0% / 100% | 41 / 150 (0 tied) | 0 / 0 / 0 |

Ablations and baselines vs RDP (c = 0):

| comparison | GMR | 95% CI | Bonferroni CI | p | win / loss | clusters / curves | exact fits both / A only / B only |
| --- | --- | --- | --- | --- | --- | --- | --- |
| vw | 1.338 | [1.292, 1.387] | [1.268, 1.411] | 1.8e-12 | 2% / 98% | 41 / 150 (0 tied) | 0 / 0 / 0 |
| uniform | 1.738 | [1.665, 1.819] | [1.631, 1.865] | 9.1e-13 | 0% / 100% | 41 / 150 (0 tied) | 0 / 0 / 0 |
| cq_unpinned | 1.735 | [1.678, 1.797] | [1.651, 1.834] | 9.1e-13 | 0% / 100% | 41 / 150 (0 tied) | 0 / 0 / 0 |
| cq_snapped | 1.996 | [1.915, 2.085] | [1.878, 2.135] | 9.1e-13 | 0% / 100% | 41 / 150 (0 tied) | 0 / 0 / 0 |
| cq_p1 | 2.160 | [2.039, 2.295] | [1.981, 2.373] | 9.1e-13 | 0% / 100% | 41 / 150 (0 tied) | 0 / 0 / 0 |

Ratio to optimal vertex subset (shortcut error, N <= 800):

| comparison | GMR | 95% CI | Bonferroni CI | p | win / loss | clusters / curves | exact fits both / A only / B only |
| --- | --- | --- | --- | --- | --- | --- | --- |
| cq | 2.754 | [2.624, 2.895] | [2.565, 2.970] | 9.1e-13 | 0% / 100% | 41 / 149 (0 tied) | 0 / 0 / 0 |
| rdp | 1.416 | [1.386, 1.448] | [1.368, 1.465] | 9.1e-13 | 0% / 100% | 41 / 149 (0 tied) | 0 / 0 / 0 |
| vw | 1.902 | [1.841, 1.969] | [1.810, 2.005] | 9.1e-13 | 0% / 100% | 41 / 149 (0 tied) | 0 / 0 / 0 |
| cq_snapped | 2.838 | [2.717, 2.967] | [2.662, 3.045] | 9.1e-13 | 0% / 100% | 41 / 149 (0 tied) | 0 / 0 / 0 |

GMR by actual budget fraction K/N (clusters): H1: <=3% 1.611 (19), 3-7% 1.751 (37), 7-15% 2.052 (40), >15% 2.161 (41); H2: <=3% 1.337 (19), 3-7% 1.488 (37), 7-15% 1.536 (40), >15% 1.437 (41); H3: <=3% 1.628 (19), 3-7% 1.429 (37), 7-15% 1.442 (40), >15% 1.791 (41)

Spearman correlation of H1 per-curve log ratio with curve properties (positive = CQ relatively worse as the property grows): N +0.21, turning +0.10, straight +0.04, kn -0.21

Failure rates (fallback to uniform): none

Mean point shortfall below K: none

Share of cq runs (c = 0) where corner pinning switched off (no corner, or more than K-2): 83.0%

Median runtime ms (Python, indicative): cq@c0 0.37, cq@c1 0.39, cq@c2 0.49, cq@c4 0.34, cq@c8 0.31, cq@c16 0.32, cq_p1@c0 0.26, cq_snapped@c0 0.40, cq_unpinned@c0 0.32, oracle@c0 48.73, rdp@c0 0.32, rdp@c1 0.33, rdp@c2 0.34, rdp@c4 0.33, rdp@c8 0.32, rdp@c16 0.32, uniform@c0 0.10, vw@c0 1.51, vw@c1 1.53, vw@c2 1.53, vw@c4 1.55, vw@c8 1.54, vw@c16 1.55

## ne_roads (exploratory)

150 curves; eligible dev/test = {'dev': 750, 'test': 1675}; clusters dev/test = {'dev': 59, 'test': 130}; too short 54220, duplicates 0; loader skips none; tuned c = {'cq': 2, 'rdp': 0, 'vw': 0}

| test | comparison | GMR | 95% CI | Bonferroni CI | p | win / loss | clusters / curves | exact fits both / A only / B only | Holm p | decision |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H1 | cq vs rdp (c0) | 1.741 | [1.635, 1.849] | [1.584, 1.907] | 1.4e-14 | 0% / 100% | 47 / 150 (0 tied) | 0 / 0 / 0 |  | (CQ worse; exploratory) |
| H2 | cq vs vw (c0) | 1.523 | [1.426, 1.629] | [1.380, 1.689] | 2e-13 | 6% / 94% | 47 / 150 (0 tied) | 0 / 0 / 0 |  | (CQ worse; exploratory) |
| H3 | cq vs rdp (tuned) | 1.390 | [1.333, 1.447] | [1.302, 1.478] | 1.4e-14 | 0% / 100% | 47 / 150 (0 tied) | 0 / 0 / 0 |  | (CQ worse; exploratory) |

Secondary metrics (same comparisons):

| comparison | GMR | 95% CI | Bonferroni CI | p | win / loss | clusters / curves | exact fits both / A only / B only |
| --- | --- | --- | --- | --- | --- | --- | --- |
| H1/rmse | 1.827 | [1.732, 1.928] | [1.685, 1.986] | 1.4e-14 | 0% / 100% | 47 / 150 (0 tied) | 0 / 0 / 0 |
| H1/shortcut | 1.741 | [1.636, 1.850] | [1.584, 1.909] | 1.4e-14 | 0% / 100% | 47 / 150 (0 tied) | 0 / 0 / 0 |
| H2/rmse | 1.760 | [1.662, 1.866] | [1.614, 1.929] | 1.4e-14 | 0% / 100% | 47 / 150 (0 tied) | 0 / 0 / 0 |
| H2/shortcut | 1.523 | [1.426, 1.630] | [1.380, 1.690] | 2e-13 | 6% / 94% | 47 / 150 (0 tied) | 0 / 0 / 0 |
| H3/rmse | 1.395 | [1.343, 1.449] | [1.315, 1.478] | 1.4e-14 | 0% / 100% | 47 / 150 (0 tied) | 0 / 0 / 0 |
| H3/shortcut | 1.390 | [1.333, 1.447] | [1.302, 1.478] | 1.4e-14 | 0% / 100% | 47 / 150 (0 tied) | 0 / 0 / 0 |

Ablations and baselines vs RDP (c = 0):

| comparison | GMR | 95% CI | Bonferroni CI | p | win / loss | clusters / curves | exact fits both / A only / B only |
| --- | --- | --- | --- | --- | --- | --- | --- |
| vw | 1.143 | [1.100, 1.189] | [1.078, 1.218] | 3e-08 | 17% / 81% | 47 / 150 (1 tied) | 0 / 0 / 0 |
| uniform | 1.856 | [1.738, 1.980] | [1.682, 2.047] | 2.8e-14 | 2% / 98% | 47 / 150 (0 tied) | 0 / 0 / 0 |
| cq_unpinned | 1.620 | [1.534, 1.710] | [1.490, 1.758] | 1.4e-14 | 0% / 100% | 47 / 150 (0 tied) | 0 / 0 / 0 |
| cq_snapped | 1.776 | [1.672, 1.882] | [1.619, 1.937] | 1.4e-14 | 0% / 100% | 47 / 150 (0 tied) | 0 / 0 / 0 |
| cq_p1 | 1.834 | [1.729, 1.948] | [1.677, 2.006] | 1.4e-14 | 0% / 100% | 47 / 150 (0 tied) | 0 / 0 / 0 |

Ratio to optimal vertex subset (shortcut error, N <= 800):

| comparison | GMR | 95% CI | Bonferroni CI | p | win / loss | clusters / curves | exact fits both / A only / B only |
| --- | --- | --- | --- | --- | --- | --- | --- |
| cq | 2.338 | [2.203, 2.477] | [2.130, 2.545] | 1.4e-14 | 0% / 100% | 47 / 150 (0 tied) | 0 / 0 / 0 |
| rdp | 1.343 | [1.309, 1.379] | [1.293, 1.397] | 1.4e-14 | 0% / 100% | 47 / 150 (0 tied) | 0 / 0 / 0 |
| vw | 1.535 | [1.465, 1.608] | [1.430, 1.648] | 1.4e-14 | 0% / 100% | 47 / 150 (0 tied) | 0 / 0 / 0 |
| cq_snapped | 2.385 | [2.261, 2.514] | [2.198, 2.582] | 1.4e-14 | 0% / 100% | 47 / 150 (0 tied) | 0 / 0 / 0 |

GMR by actual budget fraction K/N (clusters): H1: 3-7% 1.330 (6), 7-15% 1.701 (31), >15% 1.763 (47); H2: 3-7% 1.160 (6), 7-15% 1.520 (31), >15% 1.519 (47); H3: 3-7% 1.382 (6), 7-15% 1.374 (31), >15% 1.434 (47)

Spearman correlation of H1 per-curve log ratio with curve properties (positive = CQ relatively worse as the property grows): N +0.09, turning +0.15, straight +0.27, kn -0.09

Failure rates (fallback to uniform): none

Mean point shortfall below K: {'cq_snapped@c0': 0.003367003367003367}

Share of cq runs (c = 0) where corner pinning switched off (no corner, or more than K-2): 71.7%

Median runtime ms (Python, indicative): cq@c0 0.33, cq@c1 0.34, cq@c2 0.29, cq@c4 0.28, cq@c8 0.28, cq@c16 0.28, cq_p1@c0 0.24, cq_snapped@c0 0.37, cq_unpinned@c0 0.29, oracle@c0 4.77, rdp@c0 0.22, rdp@c1 0.22, rdp@c2 0.22, rdp@c4 0.23, rdp@c8 0.23, rdp@c16 0.22, uniform@c0 0.09, vw@c0 0.49, vw@c1 0.48, vw@c2 0.48, vw@c4 0.47, vw@c8 0.47, vw@c16 0.46

## omniglot (exploratory)

150 curves; eligible dev/test = {'dev': 14912, 'test': 33547}; clusters dev/test = {'dev': 494, 'test': 1126}; too short 30983, duplicates 21; loader skips none; tuned c = {'cq': 0, 'rdp': 0, 'vw': 0}

| test | comparison | GMR | 95% CI | Bonferroni CI | p | win / loss | clusters / curves | exact fits both / A only / B only | Holm p | decision |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H1 | cq vs rdp (c0) | 2.916 | [2.041, 4.453] | [1.712, 5.730] | 9.3e-19 | 4% / 89% | 128 / 150 (9 tied) | 23 / 2 / 11 |  | (CQ worse; exploratory) |
| H2 | cq vs vw (c0) | 1.720 | [1.074, 2.803] | [0.832, 3.697] | 2.1e-09 | 20% / 73% | 128 / 150 (9 tied) | 20 / 5 / 10 |  | (inconclusive; exploratory) |
| H3 | cq vs rdp (tuned) | 2.916 | [2.041, 4.453] | [1.712, 5.730] | 9.3e-19 | 4% / 89% | 128 / 150 (9 tied) | 23 / 2 / 11 |  | (CQ worse; exploratory) |

Secondary metrics (same comparisons):

| comparison | GMR | 95% CI | Bonferroni CI | p | win / loss | clusters / curves | exact fits both / A only / B only |
| --- | --- | --- | --- | --- | --- | --- | --- |
| H1/rmse | 3.221 | [2.298, 4.791] | [1.975, 6.044] | 7e-19 | 2% / 91% | 128 / 150 (9 tied) | 23 / 2 / 11 |
| H1/shortcut | 3.959 | [2.617, 6.444] | [2.255, 8.761] | 1.6e-20 | 3% / 91% | 128 / 150 (8 tied) | 19 / 0 / 15 |
| H2/rmse | 2.308 | [1.488, 3.628] | [1.175, 4.685] | 4.8e-14 | 10% / 83% | 128 / 150 (9 tied) | 20 / 5 / 10 |
| H2/shortcut | 2.340 | [1.515, 3.863] | [1.258, 5.313] | 1.5e-10 | 20% / 74% | 128 / 150 (8 tied) | 18 / 1 / 12 |
| H3/rmse | 3.221 | [2.298, 4.791] | [1.975, 6.044] | 7e-19 | 2% / 91% | 128 / 150 (9 tied) | 23 / 2 / 11 |
| H3/shortcut | 3.959 | [2.617, 6.444] | [2.255, 8.761] | 1.6e-20 | 3% / 91% | 128 / 150 (8 tied) | 19 / 0 / 15 |

Ablations and baselines vs RDP (c = 0):

| comparison | GMR | 95% CI | Bonferroni CI | p | win / loss | clusters / curves | exact fits both / A only / B only |
| --- | --- | --- | --- | --- | --- | --- | --- |
| vw | 1.695 | [1.337, 2.305] | [1.255, 2.833] | 3.9e-18 | 7% / 82% | 128 / 150 (14 tied) | 30 / 0 / 4 |
| uniform | 9.683 | [4.730, 21.239] | [3.568, 34.953] | 2.1e-22 | 1% / 98% | 128 / 150 (2 tied) | 6 / 0 / 28 |
| cq_unpinned | 7.785 | [4.252, 15.128] | [3.326, 22.318] | 4.5e-21 | 2% / 96% | 128 / 150 (2 tied) | 7 / 0 / 27 |
| cq_snapped | 2.193 | [1.755, 2.769] | [1.557, 3.147] | 2.7e-18 | 2% / 88% | 128 / 150 (12 tied) | 30 / 2 / 4 |
| cq_p1 | 8.830 | [4.629, 18.075] | [3.578, 27.846] | 6.3e-21 | 2% / 97% | 128 / 150 (2 tied) | 6 / 0 / 28 |

Ratio to optimal vertex subset (shortcut error, N <= 800):

| comparison | GMR | 95% CI | Bonferroni CI | p | win / loss | clusters / curves | exact fits both / A only / B only |
| --- | --- | --- | --- | --- | --- | --- | --- |
| cq | 5.873 | [3.902, 9.525] | [3.277, 12.760] | 2e-21 | 0% / 94% | 128 / 150 (8 tied) | 19 / 0 / 17 |
| rdp | 1.483 | [1.353, 1.688] | [1.324, 1.832] | 6.1e-21 | 0% / 91% | 128 / 150 (11 tied) | 34 / 0 / 2 |
| vw | 2.510 | [1.930, 3.471] | [1.755, 4.210] | 6.1e-21 | 0% / 91% | 128 / 150 (11 tied) | 30 / 0 / 6 |
| cq_snapped | 3.267 | [2.740, 4.008] | [2.551, 4.557] | 6.1e-21 | 0% / 91% | 128 / 150 (11 tied) | 32 / 0 / 4 |

GMR by actual budget fraction K/N (clusters): H1: <=3% 1.121 (2), 3-7% 1.788 (38), 7-15% 2.544 (97), >15% 3.850 (128); H2: <=3% 0.811 (2), 3-7% 1.190 (38), 7-15% 1.503 (97), >15% 1.722 (128); H3: <=3% 1.121 (2), 3-7% 1.788 (38), 7-15% 2.544 (97), >15% 3.850 (128)

Spearman correlation of H1 per-curve log ratio with curve properties (positive = CQ relatively worse as the property grows): N +0.28, turning +0.22, straight -0.16, kn -0.28

Failure rates (fallback to uniform): {'cq@c0': 0.015665796344647518, 'cq_snapped@c0': 0.015665796344647518, 'cq_unpinned@c0': 0.015665796344647518}

Mean point shortfall below K: {'cq_snapped@c0': 0.3289817232375979}

Share of cq runs (c = 0) where corner pinning switched off (no corner, or more than K-2): 91.4%

Median runtime ms (Python, indicative): cq@c0 0.41, cq@c1 0.54, cq@c2 0.43, cq@c4 0.40, cq@c8 0.39, cq@c16 0.37, cq_p1@c0 0.31, cq_snapped@c0 0.45, cq_unpinned@c0 0.35, oracle@c0 15.99, rdp@c0 0.30, rdp@c1 0.32, rdp@c2 0.34, rdp@c4 0.32, rdp@c8 0.34, rdp@c16 0.34, uniform@c0 0.12, vw@c0 1.08, vw@c1 1.06, vw@c2 1.08, vw@c4 1.04, vw@c8 1.10, vw@c16 1.08
