"""curvquant: curvature quantiles, curvature-mass profiles and quantile resampling."""

from .core import (
    Profile,
    curvature_median,
    curvature_quantile_points,
    point_at,
    profile,
    profile_distance,
    profile_distance_quantile_form,
    quantile_gap,
    signed_turning_angles,
    turning_angles,
)
from .metrics import approximation_errors
from .resample import (
    optimal_subset_min_eps,
    quantile_resample,
    rdp_fixed_k,
    uniform_arclength,
    vw_fixed_k,
)

__version__ = "0.1.0"
