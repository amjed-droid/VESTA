"""
Log-SNRAS Core Algorithm and Metric Suite
=========================================
Implements the variance-stabilized Log-SNRAS metric with deterministic
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


def compute_transit_shape_residuals(f_in, time_in=None, period=None, t0=None, method='poly'):
    """
    Subtracts deterministic transit profile from in-transit points to isolate
    genuine photometric noise scatter (sigma_in,res).

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

    if method == 'poly':
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


def calculate_log_snras(f_norm, transit_mask, time=None, period=None, t0=None, snr_trad=None, apply_shape_correction=True):
    """
    Calculates the Log-SNRAS metric and dispersion contrast.

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

    Returns
    -------
    dict:
        log_snras : float
        penalty : float, ln(1 + psi)
        psi : float, dispersion contrast
        d_factor : float, 1 + penalty
        sigma_in : float
        sigma_out : float
        depth : float
        n_in : int
        n_out : int
    """
    f_in = f_norm[transit_mask]
    f_out = f_norm[~transit_mask]

    n_in = len(f_in)
    n_out = len(f_out)

    if n_in == 0 or n_out == 0:
        return {
            'log_snras': np.nan, 'penalty': np.nan, 'psi': np.nan,
            'd_factor': np.nan, 'sigma_in': np.nan, 'sigma_out': np.nan,
            'depth': np.nan, 'n_in': n_in, 'n_out': n_out
        }

    depth = max(float(np.mean(f_out) - np.mean(f_in)), 0.0)
    sigma_out = float(np.std(f_out, ddof=1))

    if apply_shape_correction and n_in >= 6:
        time_in = time[transit_mask] if time is not None else None
        sigma_in, _ = compute_transit_shape_residuals(
            f_in, time_in=time_in, period=period, t0=t0
        )
    else:
        sigma_in = float(np.std(f_in, ddof=1))

    # Dispersion contrast centered at zero
    if sigma_out <= 0 or np.isnan(sigma_out):
        psi = 0.0
    else:
        psi = abs(sigma_in - sigma_out) / sigma_out

    penalty = float(np.log(1.0 + psi))
    d_factor = 1.0 + penalty

    if snr_trad is None:
        snr_trad = (depth / (sigma_out + 1e-12)) * np.sqrt(n_in)

    log_snras = snr_trad / d_factor

    return {
        'log_snras': float(log_snras),
        'penalty': float(penalty),
        'psi': float(psi),
        'd_factor': float(d_factor),
        'sigma_in': float(sigma_in),
        'sigma_out': float(sigma_out),
        'depth': float(depth),
        'n_in': n_in,
        'n_out': n_out,
        'snr_trad': float(snr_trad)
    }


def get_calibrated_tier(penalty, n_in, n_out=1000, alpha_t1=0.05, alpha_t2=0.001):
    """
    Assigns classification tier based on calibrated F-distribution percentiles.
    
    Under H0 (homoscedastic Gaussian noise), (s_in/s_out)^2 ~ F(n_in-1, n_out-1).
    psi = |s_in/s_out - 1|. We determine critical psi values from the two-tailed F-quantiles.
    """
    if n_in < 4:
        return 'Excluded (N_in < 4)'

    df1 = n_in - 1
    df2 = max(n_out - 1, 100)

    # Upper and lower quantiles for alpha_t1 (Tier 1 threshold)
    f_hi_t1 = stats.f.ppf(1 - alpha_t1 / 2, df1, df2)
    f_lo_t1 = stats.f.ppf(alpha_t1 / 2, df1, df2)
    psi_crit_t1 = max(abs(np.sqrt(f_hi_t1) - 1.0), abs(np.sqrt(f_lo_t1) - 1.0))
    p_crit_t1 = np.log(1.0 + psi_crit_t1)

    # Upper and lower quantiles for alpha_t2 (Tier 2/3 veto threshold)
    f_hi_t2 = stats.f.ppf(1 - alpha_t2 / 2, df1, df2)
    f_lo_t2 = stats.f.ppf(alpha_t2 / 2, df1, df2)
    psi_crit_t2 = max(abs(np.sqrt(f_hi_t2) - 1.0), abs(np.sqrt(f_lo_t2) - 1.0))
    p_crit_t2 = np.log(1.0 + psi_crit_t2)

    if penalty <= p_crit_t1:
        return 'Tier 1 (Clean)'
    elif penalty <= p_crit_t2:
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
