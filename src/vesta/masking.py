"""
Uniform Transit Masking Protocols
=================================
Provides consistent, unbiased in-transit and out-of-transit data partitioning
using ephemeris-based orbital phasing and buffer exclusion.
"""

import numpy as np


def create_ephemeris_mask(time, t0, period, duration_days, buffer_factor=2.0):
    """
    Creates uniform in-transit and out-of-transit masks based on orbital ephemeris.

    Parameters
    ----------
    time : np.ndarray
        Time array in days (BJD/BTJD).
    t0 : float
        Reference epoch of mid-transit in days.
    period : float
        Orbital period in days.
    duration_days : float
        Total transit duration (T_14) in days.
    buffer_factor : float
        Multiplier of duration removed around transit to prevent contamination
        of out-of-transit baseline by ingress/egress tails (default 2.0).

    Returns
    -------
    transit_mask : np.ndarray (bool)
        True for points strictly within mid-transit +/- (T_dur / 2).
    baseline_mask : np.ndarray (bool)
        True for pure out-of-transit points (excluding transit and buffer).
    """
    time = np.asarray(time)
    
    # Calculate phase in range [-0.5, 0.5)
    phase = np.mod(time - t0, period) / period
    phase[phase > 0.5] -= 1.0
    
    half_dur_phase = (duration_days / 2.0) / period
    buffer_phase = (buffer_factor * duration_days / 2.0) / period

    transit_mask = np.abs(phase) <= half_dur_phase
    buffer_zone = (np.abs(phase) > half_dur_phase) & (np.abs(phase) <= buffer_phase)
    baseline_mask = ~transit_mask & ~buffer_zone

    return transit_mask, baseline_mask


def create_event_mask(time, t_center, duration_days, buffer_factor=2.0):
    """
    Creates uniform window mask for localized events (e.g. single transit, flare, burst).
    """
    time = np.asarray(time)
    dt = np.abs(time - t_center)
    half_dur = duration_days / 2.0
    buffer_dur = buffer_factor * half_dur

    transit_mask = dt <= half_dur
    buffer_zone = (dt > half_dur) & (dt <= buffer_dur)
    baseline_mask = ~transit_mask & ~buffer_zone

    return transit_mask, baseline_mask
