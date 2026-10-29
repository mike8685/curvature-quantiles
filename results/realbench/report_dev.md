# Real-data benchmark report (dev split)

Protocol 1.0. Primary metric: directed Hausdorff error. GMR < 1 means CQ has lower error. Per-test alpha (Bonferroni over 18 planned tests) = 0.00278.

**DEV SPLIT. Used only for tuning and debugging. Not evidence.**

## ne_roads (exploratory)

150 curves; eligible dev/test = {'dev': 766, 'test': 1659}; clusters dev/test = {'dev': 766, 'test': 1659}; tuned c = None

| test | comparison | GMR | 95% CI | p | win rate | clusters / curves | Holm p | decision |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H1 | cq vs rdp (c0) | 1.734 | [1.615, 1.899] | 5.7e-26 | 4% | 150 / 150 |  | (CQ worse; exploratory) |
| H2 | cq vs vw (c0) | 1.454 | [1.354, 1.589] | 5.5e-22 | 13% | 150 / 150 |  | (CQ worse; exploratory) |
| H3 | cq vs rdp (tuned) | 1.734 | [1.615, 1.899] | 5.7e-26 | 4% | 150 / 150 |  | (CQ worse; exploratory) |

Secondary metrics (same comparisons):

| comparison / metric | GMR | 95% CI | p | win rate | clusters / curves |
| --- | --- | --- | --- | --- | --- |
| H1/rmse | 1.845 | [1.724, 2.016] | 2.6e-26 | 1% | 150 / 150 |
| H1/hausdorff_sym | 1.734 | [1.615, 1.899] | 5.7e-26 | 4% | 150 / 150 |
| H1/shortcut | 1.734 | [1.615, 1.899] | 5.7e-26 | 4% | 150 / 150 |
| H2/rmse | 1.739 | [1.626, 1.903] | 2.3e-26 | 1% | 150 / 150 |
| H2/hausdorff_sym | 1.454 | [1.354, 1.589] | 5.5e-22 | 13% | 150 / 150 |
| H2/shortcut | 1.454 | [1.354, 1.589] | 5.6e-22 | 13% | 150 / 150 |
| H3/rmse | 1.845 | [1.724, 2.016] | 2.6e-26 | 1% | 150 / 150 |
| H3/hausdorff_sym | 1.734 | [1.615, 1.899] | 5.7e-26 | 4% | 150 / 150 |
| H3/shortcut | 1.734 | [1.615, 1.899] | 5.7e-26 | 4% | 150 / 150 |

Ablations and baselines vs RDP (c = 0):

| method / rdp | GMR | 95% CI | p | win rate | clusters / curves |
| --- | --- | --- | --- | --- | --- |
| vw | 1.193 | [1.155, 1.233] | 7.8e-17 | 19% | 150 / 150 |
| uniform | 1.881 | [1.750, 2.081] | 2.3e-26 | 1% | 150 / 150 |
| cq_unpinned | 1.602 | [1.500, 1.754] | 1e-25 | 3% | 150 / 150 |
| cq_snapped | 1.759 | [1.663, 1.870] | 3.2e-26 | 2% | 150 / 150 |
| cq_p1 | 1.908 | [1.779, 2.092] | 2.8e-26 | 1% | 150 / 150 |

Ratio to optimal vertex subset (shortcut error, N <= 800):

| method / oracle | GMR | 95% CI | p | win rate | clusters / curves |
| --- | --- | --- | --- | --- | --- |
| cq | 2.346 | [2.186, 2.569] | 2.3e-26 | 1% | 150 / 150 |
| rdp | 1.353 | [1.327, 1.378] | 4.9e-26 | 0% | 150 / 150 |
| vw | 1.614 | [1.563, 1.666] | 2.3e-26 | 0% | 150 / 150 |
| cq_snapped | 2.380 | [2.252, 2.528] | 2.3e-26 | 0% | 150 / 150 |

GMR by budget ratio r: H1: r=0.02 1.730, r=0.1 2.038, r=0.2 1.756; H2: r=0.02 1.480, r=0.1 1.496, r=0.2 1.465; H3: r=0.02 1.730, r=0.1 2.038, r=0.2 1.756

Failure rates (fallback to uniform): none

Mean point shortfall below K: {'cq_snapped@c0': 0.020689655172413793}

Median runtime ms (Python, indicative): cq@c0 0.36, cq@c1 0.42, cq@c2 0.33, cq@c4 0.33, cq@c8 0.33, cq@c16 0.32, cq_p1@c0 0.29, cq_snapped@c0 0.40, cq_unpinned@c0 0.31, oracle@c0 4.63, rdp@c0 0.24, rdp@c1 0.23, rdp@c2 0.23, rdp@c4 0.22, rdp@c8 0.21, rdp@c16 0.21, uniform@c0 0.10, vw@c0 0.50, vw@c1 0.48, vw@c2 0.48, vw@c4 0.49, vw@c8 0.48, vw@c16 0.47
