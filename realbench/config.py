"""Frozen parameters of the real-data benchmark. See PROTOCOL.md; every value here is quoted there."""

import math

PROTOCOL_VERSION = "1.0"
SEED = 20261028
SPLIT_SALT = "realbench-v1"
DEV_FRACTION = 0.30
CAP = {"dev": 150, "test": 400}           # max curves per split (seeded simple random sample)

# inclusion
MIN_VERTICES = 40
CAD_SAMPLES_PER_LOOP = 1000              # uniform arc-length spacing = loop length / this

# budgets
BUDGET_RATIOS = (0.02, 0.05, 0.10, 0.20)
K_MIN = 8
K_MAX_FRACTION = 0.5

# preprocessing (the only tuned parameter, same grid for every shape-aware method)
SMOOTHING_GRID = (0, 1, 2, 4, 8, 16)     # sigma in median edge lengths
TUNED_METHODS = ("cq", "rdp", "vw")

# methods (published settings)
CQ_P = 0.5
CQ_CORNER_THRESHOLD = math.pi / 6
ORACLE_MAX_N = 800

# metrics
EVAL_SAMPLES_FACTOR = 4
EVAL_SAMPLES_MIN = 2000
EVAL_SAMPLES_MAX = 200_000
PRIMARY_METRIC = "hausdorff"
SECONDARY_METRICS = ("rmse", "hausdorff_sym", "shortcut")
ZERO_FLOOR = 1e-12                       # x curve diameter
ONCURVE_TOL = 1e-9                       # x curve diameter

# statistics
CONFIRMATORY_COMPARISONS = (
    ("H1", "cq", "rdp", "c0"),
    ("H2", "cq", "vw", "c0"),
    ("H3", "cq", "rdp", "tuned"),
)
PLANNED_CONFIRMATORY_DATASETS = 6
ALPHA_FAMILY = 0.05
N_PLANNED_TESTS = PLANNED_CONFIRMATORY_DATASETS * len(CONFIRMATORY_COMPARISONS)
ALPHA_PER_TEST = ALPHA_FAMILY / N_PLANNED_TESTS
BOOTSTRAP_RESAMPLES = 10_000
EQUIVALENCE_MARGIN = 1.05
MIN_TEST_CLUSTERS = 30

# confirmatory datasets: the maps and gpx are fixed now; exactly one pen-stroke and one
# CAD dataset are named with `freeze DATASET --confirmatory`, before their test split runs
FIXED_CONFIRMATORY = ("ne_coastline", "ne_rivers", "ne_roads", "gpx")
LATE_CONFIRMATORY_CATEGORIES = ("pen strokes", "CAD")
