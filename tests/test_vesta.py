"""
Unit Tests for VESTA (Variance-Evaluated Stability for Transit Analysis)
=======================================================================
Validates core statistical invariance, non-parametric shape detrending,
and calibrated F-distribution tiers across synthetic and archival regimes.
"""

import os
import sys
import unittest
import numpy as np

sys.path.insert(0, os.path.abspath('src'))
from vesta import (
    calculate_vesta,
    compute_shape_corrected_sigma_in,
    get_calibrated_tier,
    compute_baseline_metrics
)


class TestVESTAFramework(unittest.TestCase):

    def test_homoscedastic_gaussian_transit(self):
        """Under pure Gaussian homoscedastic noise with a real transit, VESTA preserves matched-filter SNR."""
        rng = np.random.default_rng(42)
        n_points = 5000
        sigma_0 = 0.0002
        flux = 1.0 + rng.normal(0, sigma_0, size=n_points)
        mask = np.zeros(n_points, dtype=bool)
        mask[2450:2550] = True  # N_in = 100
        flux[2450:2550] -= 0.002  # 2,000 ppm genuine transit

        res = calculate_vesta(flux, mask, apply_shape_correction=True)
        tier = get_calibrated_tier(res['penalty'], res['n_in'], res['n_out'])

        self.assertLess(res['psi'], 0.20, f"Expected small contrast, got {res['psi']}")
        self.assertLess(res['penalty'], 0.15, f"Expected low penalty, got {res['penalty']}")
        self.assertIn('Tier 1', tier, f"Expected Tier 1 Clean, got {tier}")
        # Verify score preservation
        score_ratio = res['vesta_score'] / (res['snr_trad'] + 1e-12)
        self.assertGreater(score_ratio, 0.85, f"Expected high score preservation, got {score_ratio}")

    def test_deep_transit_shape_correction(self):
        """Deep transits (e.g. 30,000 ppm) must not trigger false penalties when corrected."""
        rng = np.random.default_rng(42)
        n_points = 5000
        noise_sigma = 0.001
        flux = 1.0 + rng.normal(0, noise_sigma, size=n_points)

        # Ingress/egress transit dip
        transit_slice = slice(2400, 2600)  # N_in = 200
        depth = 0.030  # 30,000 ppm
        x = np.linspace(-1, 1, 200)
        transit_profile = -depth * (1.0 - 0.2 * x**2)
        flux[transit_slice] += transit_profile

        mask = np.zeros(n_points, dtype=bool)
        mask[transit_slice] = True

        res_raw = calculate_vesta(flux, mask, apply_shape_correction=False)
        res_corr = calculate_vesta(flux, mask, apply_shape_correction=True)

        tier_raw = get_calibrated_tier(res_raw['penalty'], res_raw['n_in'], res_raw['n_out'])
        tier_corr = get_calibrated_tier(res_corr['penalty'], res_corr['n_in'], res_corr['n_out'])

        self.assertGreater(res_raw['penalty'], 0.30, "Raw penalty should be elevated due to geometry")
        self.assertLess(res_corr['penalty'], 0.15, f"Corrected penalty should be near zero, got {res_corr['penalty']}")
        self.assertIn('Tier 1', tier_corr, f"Corrected transit should be Tier 1 Clean, got {tier_corr}")

    def test_heteroscedastic_noise_burst_veto(self):
        """Local non-stationary noise bursts must trigger automated Tier 3 Veto."""
        rng = np.random.default_rng(42)
        n_points = 5000
        flux = 1.0 + rng.normal(0, 0.0005, size=n_points)

        # Coincident severe noise burst inside candidate window
        transit_slice = slice(2450, 2550)  # N_in = 100
        flux[transit_slice] += rng.normal(-0.002, 0.0035, size=100)  # 7x variance burst

        mask = np.zeros(n_points, dtype=bool)
        mask[transit_slice] = True

        res = calculate_vesta(flux, mask, apply_shape_correction=True)
        tier = get_calibrated_tier(res['penalty'], res['n_in'], res['n_out'])

        self.assertGreater(res['penalty'], 0.60, f"Expected heavy penalty, got {res['penalty']}")
        self.assertIn('Tier 3', tier, f"Expected Tier 3 Veto, got {tier}")
        # Score must be heavily attenuated
        self.assertLess(res['vesta_score'], res['snr_trad'] / 1.5)

    def test_depth_invariance_across_scales(self):
        """Verify consistent near-zero penalties across 3 decades of synthetic depths."""
        rng = np.random.default_rng(100)
        n_points = 5000
        noise_sigma = 0.0003

        depths = [0.030, 0.003, 0.0003]  # 30,000 ppm, 3,000 ppm, 300 ppm
        penalties = []

        for d in depths:
            flux = 1.0 + rng.normal(0, noise_sigma, size=n_points)
            transit_slice = slice(2400, 2550)  # N_in = 150
            flux[transit_slice] -= d
            mask = np.zeros(n_points, dtype=bool)
            mask[transit_slice] = True

            res = calculate_vesta(flux, mask, apply_shape_correction=True)
            penalties.append(res['penalty'])

        for p, d in zip(penalties, depths):
            self.assertLess(p, 0.15, f"Depth {d*1e6:.0f} ppm produced excessive penalty {p:.4f}")


if __name__ == '__main__':
    unittest.main(verbosity=2)
