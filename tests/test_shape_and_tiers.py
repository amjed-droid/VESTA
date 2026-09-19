"""
Unit Tests and Synthetic Benchmark Validation for Log-SNRAS
"""

import sys
import os
sys.path.insert(0, os.path.abspath('src'))

import numpy as np
from log_snras import (
    calculate_log_snras,
    compute_shape_corrected_sigma_in,
    get_calibrated_tier,
    compute_baseline_metrics
)


def test_homoscedastic_limit():
    """Homoscedastic Gaussian noise: Psi should approach 0, penalty -> 0, Tier 1."""
    rng = np.random.default_rng(42)
    n_points = 5000
    flux = 1.0 + rng.normal(0, 0.001, size=n_points)
    mask = np.zeros(n_points, dtype=bool)
    mask[2450:2550] = True  # N_in = 100

    res = calculate_log_snras(flux, mask, apply_shape_correction=False)
    tier = get_calibrated_tier(res['penalty'], res['n_in'], res['n_out'])
    
    print(f"Homoscedastic test: Psi = {res['psi']:.4f}, Penalty = {res['penalty']:.4f}, Tier = {tier}")
    assert res['psi'] < 0.20, f"Expected small psi, got {res['psi']}"
    assert 'Tier 1' in tier, f"Expected Tier 1, got {tier}"


def test_deep_transit_shape_correction():
    """Deep transit (e.g. TOI-201, 5000 ppm) with and without shape correction."""
    rng = np.random.default_rng(42)
    n_points = 5000
    noise_sigma = 0.0005
    flux = 1.0 + rng.normal(0, noise_sigma, size=n_points)
    
    # Add deep U-shaped / box transit with curved ingress/egress
    transit_slice = slice(2400, 2550) # N_in = 150
    depth = 0.0053 # 5300 ppm
    x = np.linspace(-1, 1, 150)
    transit_shape = -depth * (1.0 - 0.3 * x**2)
    flux[transit_slice] += transit_shape

    mask = np.zeros(n_points, dtype=bool)
    mask[transit_slice] = True

    # Without shape correction: geometric gradient inflates sigma_in
    res_uncorrected = calculate_log_snras(flux, mask, apply_shape_correction=False)
    tier_uncorrected = get_calibrated_tier(res_uncorrected['penalty'], res_uncorrected['n_in'])

    # With shape correction: geometric gradient is removed
    res_corrected = calculate_log_snras(flux, mask, apply_shape_correction=True)
    tier_corrected = get_calibrated_tier(res_corrected['penalty'], res_corrected['n_in'])

    print(f"Deep transit uncorrected: sigma_in = {res_uncorrected['sigma_in']:.6f}, Penalty = {res_uncorrected['penalty']:.4f}, Tier = {tier_uncorrected}")
    print(f"Deep transit corrected  : sigma_in = {res_corrected['sigma_in']:.6f}, Penalty = {res_corrected['penalty']:.4f}, Tier = {tier_corrected}")

    assert res_uncorrected['penalty'] > res_corrected['penalty'], "Shape correction should reduce penalty on deep transits"
    assert 'Tier 1' in tier_corrected, f"Shape corrected deep transit should be Tier 1, got {tier_corrected}"


def test_heteroscedastic_noise_burst():
    """Noise burst (artifact): in-transit variance is 5x higher than baseline -> Tier 3 Veto."""
    rng = np.random.default_rng(42)
    n_points = 5000
    flux = 1.0 + rng.normal(0, 0.001, size=n_points)
    mask = np.zeros(n_points, dtype=bool)
    mask[2450:2550] = True
    flux[mask] += rng.normal(0, 0.005, size=100) # 5x variance burst

    res = calculate_log_snras(flux, mask, apply_shape_correction=True)
    tier = get_calibrated_tier(res['penalty'], res['n_in'], res['n_out'])
    print(f"Heteroscedastic burst: Psi = {res['psi']:.4f}, Penalty = {res['penalty']:.4f}, Tier = {tier}")
    assert 'Tier 3' in tier, f"Expected Tier 3 (Veto), got {tier}"


if __name__ == '__main__':
    test_homoscedastic_limit()
    test_deep_transit_shape_correction()
    test_heteroscedastic_noise_burst()
    print("ALL TESTS PASSED SUCCESSFULLY!")
