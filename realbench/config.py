"""Frozen parameters of the real-data benchmark. See PROTOCOL.md; every value here is quoted there."""

import math

PROTOCOL_VERSION = "1.1"
SEED = 20261028
SPLIT_SALT = "realbench-v1"
DEV_FRACTION = 0.30
CAP = {"dev": 150, "test": 400}           # max curves per split (hash-ranked simple random sample)

# inclusion
MIN_VERTICES = 40
CAD_SAMPLES_PER_LOOP = 1000              # uniform arc-length spacing = loop length / this
CLOSED_TOL = 1e-6                        # x diameter: a loop whose ends are this close is closed (and snapped shut)
MAP_CELL_DEG = 10.0                      # map clusters: 10 x 10 degree cell of the (name-group) centroid

# Quick Draw categories, fixed before any pen data was seen: the 10 of the 345 official
# categories with the smallest sha256("realbench-v1:quickdraw:" + name)
QUICKDRAW_CATEGORIES = ("moustache", "sun", "binoculars", "stop sign", "paint can",
                        "strawberry", "crocodile", "bandage", "ambulance", "nose")

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
SECONDARY_METRICS = ("rmse", "shortcut")
ERROR_FLOOR = 1e-9                       # x curve diameter: errors below this count as exact fits (ties)
ONCURVE_TOL = 1e-9                       # x curve diameter

# statistics
CONFIRMATORY_COMPARISONS = (
    ("H1", "cq", "rdp", "c0"),
    ("H2", "cq", "vw", "c0"),
    ("H3", "cq", "rdp", "tuned"),
)
# all six confirmatory datasets are named now, before any of their data is in hand
CONFIRMATORY_DATASETS = ("ne_coastline", "ne_rivers", "ne_roads", "geolife", "quickdraw_raw", "fusion360")
ALPHA_FAMILY = 0.05
N_PLANNED_TESTS = len(CONFIRMATORY_DATASETS) * len(CONFIRMATORY_COMPARISONS)
ALPHA_PER_TEST = ALPHA_FAMILY / N_PLANNED_TESTS
BOOTSTRAP_RESAMPLES = 20_000
EQUIVALENCE_MARGIN = 1.05
MIN_TEST_CLUSTERS = 30

# software versions recorded at freeze and required at test time
TRACKED_PACKAGES = ("numpy", "scipy", "pandas", "shapely", "gpxpy", "ezdxf")
