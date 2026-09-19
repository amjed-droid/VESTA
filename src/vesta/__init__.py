"""
VESTA: Variance-Evaluated Stability for Transit Analysis
========================================================
A Depth-Invariant, Variance-Stabilized Framework for
Population-Scale Exoplanet Transit Vetting.

Author: Ahmed Sattar Jabbar (ahmed.state.me@gmail.com)
Journal: Astronomy and Computing (ASCOM)
"""

from log_snras.core import (
    calculate_log_snras as calculate_vesta,
    calculate_log_snras,
    compute_shape_corrected_sigma_in,
    get_calibrated_tier,
    compute_baseline_metrics
)
from log_snras.masking import create_ephemeris_mask
from log_snras.stats import compute_auc, stratified_bootstrap_ci, delong_test

__version__ = "3.0.0"
__all__ = [
    "calculate_vesta",
    "calculate_log_snras",
    "compute_shape_corrected_sigma_in",
    "get_calibrated_tier",
    "compute_baseline_metrics",
    "create_ephemeris_mask",
    "compute_auc",
    "stratified_bootstrap_ci",
    "delong_test"
]
