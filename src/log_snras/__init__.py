"""Log-SNRAS package initialization."""
from .core import (
    calculate_log_snras,
    compute_shape_corrected_sigma_in,
    get_calibrated_tier,
    compute_baseline_metrics
)

__version__ = "2.0.0"
__all__ = [
    "calculate_log_snras",
    "compute_shape_corrected_sigma_in",
    "get_calibrated_tier",
    "compute_baseline_metrics"
]
