"""
Light Curve Ingestion and FITS Reader
=====================================
Safely ingests Kepler and TESS light curves using explicit FITS column names
rather than arbitrary column indices.
"""

import numpy as np
from astropy.io import fits


FLUX_CANDIDATE_COLUMNS = [
    'PDCSAP_FLUX', 'SAP_FLUX', 'DVT_FLUX', 'FLUX', 'CAL_FLUX', 'CORR_FLUX'
]

TIME_CANDIDATE_COLUMNS = [
    'TIME', 'BTJD', 'BKJD', 'BJD', 'TIME_UTC'
]


def load_fits_light_curve(file_path, preferred_flux='PDCSAP_FLUX'):
    """
    Loads time and flux arrays from a standard FITS binary table.

    Parameters
    ----------
    file_path : str
        Path to the FITS file.
    preferred_flux : str
        First-choice flux column name.

    Returns
    -------
    time : np.ndarray
    flux : np.ndarray (normalized)
    meta : dict
    """
    with fits.open(file_path) as hdul:
        # Find binary table extension (typically ext 1)
        bintable_ext = None
        for ext in hdul:
            if isinstance(ext, (fits.BinTableHDU, fits.TableHDU)):
                bintable_ext = ext
                break

        if bintable_ext is None:
            raise ValueError(f"No binary table found in {file_path}")

        data = bintable_ext.data
        col_names = [c.name.upper() for c in bintable_ext.columns]

        # Locate time column
        time_col = None
        for cand in TIME_CANDIDATE_COLUMNS:
            if cand in col_names:
                time_col = cand
                break
        if time_col is None:
            time_col = col_names[0]

        # Locate flux column
        flux_col = None
        candidates = [preferred_flux] + [c for c in FLUX_CANDIDATE_COLUMNS if c != preferred_flux]
        for cand in candidates:
            if cand in col_names:
                flux_col = cand
                break
        if flux_col is None:
            raise KeyError(f"Could not locate a valid flux column in {file_path}. Columns available: {col_names}")

        raw_time = np.array(data[time_col], dtype=float)
        raw_flux = np.array(data[flux_col], dtype=float)

        # Quality filtering and NaN removal
        valid = np.isfinite(raw_time) & np.isfinite(raw_flux) & (raw_flux > 0)
        time = raw_time[valid]
        flux = raw_flux[valid]

        if len(flux) == 0:
            raise ValueError(f"No finite positive flux points in {file_path}")

        # Normalize by median
        median_flux = np.median(flux)
        f_norm = flux / median_flux

        header = bintable_ext.header
        meta = {
            'target_name': header.get('OBJECT', header.get('TARGNAME', 'UNKNOWN')),
            'teff': header.get('TEFF', np.nan),
            'time_col': time_col,
            'flux_col': flux_col,
            'n_points': len(f_norm),
            'median_flux': median_flux
        }

    return time, f_norm, meta
