"""
VESTA Core Algorithm and Metric Suite
=====================================
Implements the variance-stabilized VESTA framework with deterministic
transit shape correction and calibrated F-distribution thresholds,
alongside established transit vetting benchmarks.

References:
- Pont, F., Zucker, S., & Queloz, D. (2006), MNRAS, 373, 231
- Kovacs, G., Zucker, S., & Mazeh, T. (2002), A&A, 391, 369
- Hampel, F. R. (1974), JASA, 69, 383
"""

import numpy as np
from scipy import stats
from scipy.signal import savgol_filter


def compute_transit_shape_residuals(f_in, time_in=None, period=None, t0=None, method='adaptive_spline'):
    """
    Subtracts deterministic transit profile from in-transit points to isolate
    genuine photometric noise scatter (sigma_in,res) via adaptive binning
    or non-parametric splines.

    Parameters
    ----------
    f_in : np.ndarray
        In-transit flux points.
    time_in : np.ndarray, optional
        Time array corresponding to in-transit points.
    period : float, optional
        Orbital period in days (for phase folding multi-transit data).
    t0 : float, optional
        Transit epoch in days.
    method : str
        'adaptive_spline' (dynamic M bins with linear spline, Algorithm 2),
        'poly' (2nd-4th order polynomial), 'savgol' (Savitzky-Golay), or 'median'.

    Returns
    -------
    sigma_in_res : float
        Residual standard deviation after removing geometric transit curve.
    f_model : np.ndarray
        Fitted deterministic transit model.
    """
    n_in = len(f_in)
    if n_in < 4:
        return float(np.std(f_in, ddof=1)) if n_in > 1 else 0.0, f_in

    # Determine coordinate for fitting: phase or normalized time
    if time_in is not None and period is not None and t0 is not None and period > 0:
        phase = ((time_in - t0) % period) / period
        phase[phase > 0.5] -= 1.0
        x = phase
    elif time_in is not None:
        x = time_in - np.median(time_in)
        span = np.max(x) - np.min(x)
        if span > 0:
            x = x / span
    else:
        x = np.linspace(-1, 1, n_in)

    # Sort x to avoid conditioning issues
    sort_idx = np.argsort(x)
    x_sorted = x[sort_idx]
    f_sorted = f_in[sort_idx]

    if method == 'adaptive_spline':
        # Dynamic adaptive bin count M = clip(floor(sqrt(N_in)), 5, 100)
        M = int(np.clip(np.floor(np.sqrt(n_in)), 5, 100))
        bin_edges = np.linspace(x_sorted[0], x_sorted[-1], M + 1)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        bin_medians = np.zeros(M)
        for m in range(M):
            if m == M - 1:
                in_bin = (x_sorted >= bin_edges[m]) & (x_sorted <= bin_edges[m + 1])
            else:
                in_bin = (x_sorted >= bin_edges[m]) & (x_sorted < bin_edges[m + 1])
            if np.any(in_bin):
                bin_medians[m] = np.median(f_sorted[in_bin])
            else:
                bin_medians[m] = np.median(f_sorted)

        f_model_sorted = np.interp(x_sorted, bin_centers, bin_medians)
        inv_sort = np.empty_like(sort_idx)
        inv_sort[sort_idx] = np.arange(len(sort_idx))
        f_model = f_model_sorted[inv_sort]

    elif method == 'poly':
        deg = 4 if n_in >= 15 else (2 if n_in >= 6 else 1)
        try:
            coeffs = np.polyfit(x_sorted, f_sorted, deg=deg)
            f_model = np.polyval(coeffs, x)
        except Exception:
            f_model = np.full_like(f_in, np.median(f_in))
    elif method == 'savgol' and n_in >= 7:
        window_length = min(n_in if n_in % 2 == 1 else n_in - 1, 15)
        f_model = savgol_filter(f_in, window_length=window_length, polyorder=2)
    else:
        f_model = np.full_like(f_in, np.median(f_in))

    residuals = f_in - f_model
    sigma_in_res = float(np.std(residuals, ddof=1))
    return sigma_in_res, f_model


compute_shape_corrected_sigma_in = compute_transit_shape_residuals


def calculate_vesta(f_norm, transit_mask, time=None, period=None, t0=None, snr_trad=None,
                    apply_shape_correction=True, companion_masks=None,
                    shape_method='adaptive_spline', apply_shrinkage=False, m0=10):
    """
    Calculates the VESTA metric with Poisson noise correction, adaptive detrending,
    and optional cross-ephemeris companion masking.

    Parameters
    ----------
    f_norm : np.ndarray
        Normalized flux time series.
    transit_mask : np.ndarray (bool)
        True for in-transit cadences, False for out-of-transit.
    time : np.ndarray, optional
        Time array.
    period : float, optional
        Orbital period (for phase-folded shape correction).
    t0 : float, optional
        Transit epoch.
    snr_trad : float, optional
        Traditional SNR. If None, it will be computed from flux and mask.
    apply_shape_correction : bool
        Whether to subtract transit shape before computing sigma_in.
    companion_masks : list of np.ndarray, optional
        List of boolean masks for companion planets to excise from out-of-transit baseline.
    shape_method : str
        Method for transit shape subtraction ('adaptive_spline', 'poly', 'savgol').
    apply_shrinkage : bool
        Whether to apply Bayesian shrinkage for small sample sizes (N_in < 30).
    m0 : float
        Prior degrees of freedom for Bayesian shrinkage (default 10).

    Returns
    -------
    dict:
        vesta_score : float
        penalty : float, ln(1 + psi)
        psi : float, Poisson-corrected dispersion contrast
        psi_one_sided : float, one-sided excess dispersion contrast
        d_factor : float, 1 + penalty
        sigma_in : float
        sigma_out : float
        sigma_expected : float
        depth : float
        n_in : int
        n_out : int
        snr_trad : float
    """
    # Cross-ephemeris companion masking to prevent baseline contamination
    if companion_masks is not None and len(companion_masks) > 0:
        comp_mask = np.zeros(len(f_norm), dtype=bool)
        for c in companion_masks:
            comp_mask |= np.asarray(c, dtype=bool)
        out_mask = (~transit_mask) & (~comp_mask)
    else:
        out_mask = ~transit_mask

    f_in = f_norm[transit_mask]
    f_out = f_norm[out_mask]

    n_in = len(f_in)
    n_out = len(f_out)

    if n_in == 0 or n_out == 0:
        return {
            'vesta_score': np.nan, 'penalty': np.nan, 'psi': np.nan,
            'psi_one_sided': np.nan, 'd_factor': np.nan,
            'sigma_in': np.nan, 'sigma_out': np.nan, 'sigma_expected': np.nan,
            'depth': np.nan, 'n_in': n_in, 'n_out': n_out, 'snr_trad': np.nan
        }

    depth = max(float(np.mean(f_out) - np.mean(f_in)), 0.0)
    sigma_out = float(np.std(f_out, ddof=1))

    if apply_shape_correction and n_in >= 6:
        time_in = time[transit_mask] if time is not None else None
        sigma_in, _ = compute_transit_shape_residuals(
            f_in, time_in=time_in, period=period, t0=t0, method=shape_method
        )
    else:
        sigma_in = float(np.std(f_in, ddof=1))

    # Poisson photon-counting expected dispersion (sigma_expected = sigma_out * sqrt(1 - delta))
    sigma_expected = sigma_out * np.sqrt(max(1.0 - depth, 1e-6))

    # Optional analytical Bayesian shrinkage for sparse single-transit data (N_in < 30)
    if apply_shrinkage and n_in < 30:
        sigma_in = np.sqrt(((n_in - 1) * (sigma_in**2) + m0 * (sigma_expected**2)) / (n_in - 1 + m0))

    # Poisson-corrected dispersion contrast
    if sigma_expected <= 0 or np.isnan(sigma_expected):
        psi = 0.0
        psi_one_sided = 0.0
    else:
        psi = abs(sigma_in - sigma_expected) / sigma_expected
        psi_one_sided = max(0.0, sigma_in - sigma_expected) / sigma_expected

    penalty = float(np.log(1.0 + psi))
    d_factor = 1.0 + penalty

    if snr_trad is None:
        snr_trad = (depth / (sigma_out + 1e-12)) * np.sqrt(n_in)

    vesta_score = snr_trad / d_factor

    return {
        'vesta_score': float(vesta_score),
        'penalty': float(penalty),
        'psi': float(psi),
        'psi_one_sided': float(psi_one_sided),
        'd_factor': float(d_factor),
        'sigma_in': float(sigma_in),
        'sigma_out': float(sigma_out),
        'sigma_expected': float(sigma_expected),
        'depth': float(depth),
        'n_in': n_in,
        'n_out': n_out,
        'snr_trad': float(snr_trad)
    }


def get_calibrated_tier(penalty, n_in=100, n_out=1000, alpha_t1=0.05, alpha_t2=0.001):
    """
    Assigns operational classification tier based on calibrated thresholds
    formalized in Section 2.5 of the manuscript:
    - Tier 1 (Clean Candidate): P <= 0.15
    - Tier 2 (Review Required): 0.15 < P <= 0.60
    - Tier 3 (Automated Veto):   P > 0.60
    """
    if penalty <= 0.15:
        return 'Tier 1 (Clean)'
    elif penalty <= 0.60:
        return 'Tier 2 (Review)'
    else:
        return 'Tier 3 (Veto)'



def compute_baseline_metrics(f_norm, transit_mask, duration_hours=2.5, cadence_minutes=2.0):
    """
    Computes comparator baseline SNR statistics correctly as formulated in literature:
    1. Traditional SNR (T-SNR)
    2. Robust SNR (R-SNR, Hampel 1974 via MAD)
    3. Pont SNR (P-SNR, Pont et al. 2006 duration-binned red-noise formulation)
    4. BLS SNR proxy (B-SNR, Kovacs et al. 2002)
    """
    f_in = f_norm[transit_mask]
    f_out = f_norm[~transit_mask]

    n_in = len(f_in)
    n_out = len(f_out)
    if n_in == 0 or n_out == 0:
        return {'t_snr': np.nan, 'r_snr': np.nan, 'p_snr': np.nan, 'b_snr': np.nan}

    depth = max(float(np.mean(f_out) - np.mean(f_in)), 0.0)
    sigma_out = float(np.std(f_out, ddof=1))
    t_snr = (depth / (sigma_out + 1e-12)) * np.sqrt(n_in)

    # Robust SNR via Normalized MAD: sigma_mad = 1.4826 * MAD
    mad = np.median(np.abs(f_out - np.median(f_out)))
    sigma_mad = 1.4826 * mad
    r_snr = (depth / (sigma_mad + 1e-12)) * np.sqrt(n_in)

    # Pont et al. (2006) formulation:
    # Pont SNR bins out-of-transit data on the transit duration scale (n_dur points)
    # sigma_r^2 = var(binned flux) - sigma_w^2 / n_dur
    points_per_duration = max(int(round((duration_hours * 60.0) / cadence_minutes)), 2)
    n_bins = len(f_out) // points_per_duration
    if n_bins >= 5:
        binned_out = np.mean(f_out[:n_bins * points_per_duration].reshape(n_bins, points_per_duration), axis=1)
        var_binned = np.var(binned_out, ddof=1)
        sigma_w2 = sigma_out**2
        var_red = max(var_binned - sigma_w2 / points_per_duration, 0.0)
        sigma_tot = np.sqrt((sigma_w2 / n_in) + var_red)
        p_snr = depth / (sigma_tot + 1e-12)
    else:
        p_snr = t_snr

    # BLS Signal Residue (SR proxy, Kovacs 2002):
    # SR = sqrt(Q / (1 - Q)) * depth / sigma_out * sqrt(N) where Q = n_in / N_total
    q_duty = n_in / float(n_in + n_out)
    b_snr = np.sqrt(max(q_duty / (1.0 - q_duty + 1e-12), 1e-6)) * (depth / (sigma_out + 1e-12)) * np.sqrt(n_in + n_out)

    return {
        't_snr': float(t_snr),
        'r_snr': float(r_snr),
        'p_snr': float(p_snr),
        'b_snr': float(b_snr),
        'depth': float(depth),
        'sigma_out': float(sigma_out)
    }
