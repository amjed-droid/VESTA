import os
import json
import urllib.request
import urllib.parse
import numpy as np
import pandas as pd
from astropy.io import fits
import sys

# Add src to path
sys.path.insert(0, os.path.abspath('.'))
from src.log_snras.core import calculate_log_snras, get_calibrated_tier, compute_baseline_metrics
from src.log_snras.masking import create_ephemeris_mask


def get_ephemerides(planet_names):
    """Query NASA Exoplanet Archive for ephemerides."""
    names_str = "','".join([p.replace("'", "''") for p in planet_names])
    query = f"select pl_name,hostname,pl_orbper,pl_tranmid,pl_trandur,pl_trandep from ps where pl_name in ('{names_str}')"
    url = 'https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=' + urllib.parse.quote_plus(query) + '&format=json'
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode())
            ephems = {}
            for row in data:
                name = row['pl_name']
                if name not in ephems and row.get('pl_orbper') is not None and row.get('pl_tranmid') is not None:
                    ephems[name] = {
                        'period': float(row['pl_orbper']),
                        't0': float(row['pl_tranmid']),
                        'dur_hours': float(row['pl_trandur']) if row.get('pl_trandur') else 2.5,
                        'depth_ppm': float(row['pl_trandep']) * 10000 if row.get('pl_trandep') else 10000.0
                    }
            return ephems
    except Exception as e:
        print(f"Query error: {e}")
        return {}


def test_kepler_targets():
    print("\n" + "=" * 75)
    print("TESTING DEEP TRANSITING PLANETS FROM D:/nrebig (KEPLER DATASETS)")
    print("=" * 75)

    targets_data = json.load(open('D:/nrebig/kepler/targets.json', encoding='utf-8'))
    kepler_dirs = [d for d in os.listdir('D:/nrebig/kepler/lightcurves') if os.path.isdir(os.path.join('D:/nrebig/kepler/lightcurves', d))]
    
    planets_to_test = []
    dir_map = {}
    for d in kepler_dirs:
        rank_str = d.split('_')[0]
        try:
            rank = int(rank_str)
            target_info = next((t for t in targets_data if t['rank'] == rank), None)
            if target_info:
                planets_to_test.append(target_info['planet'])
                dir_map[target_info['planet']] = d
        except ValueError:
            pass

    print(f"Found {len(planets_to_test)} downloaded Kepler target folders.")
    ephems = get_ephemerides(planets_to_test[:15])
    print(f"Retrieved NASA ephemerides for {len(ephems)} targets.\n")

    results = []
    for planet_name, ephem in ephems.items():
        folder_name = dir_map[planet_name]
        folder_path = os.path.join('D:/nrebig/kepler/lightcurves', folder_name)
        
        # Find all fits files recursively
        fits_files = []
        for root, _, files in os.walk(folder_path):
            for f in files:
                if f.endswith('_llc.fits'): # Long cadence
                    fits_files.append(os.path.join(root, f))
        if not fits_files:
            continue

        fits_file = fits_files[0]
        try:
            with fits.open(fits_file) as hdul:
                data = hdul[1].data
                time = data['TIME']
                flux = data['PDCSAP_FLUX']
                
                valid = np.isfinite(time) & np.isfinite(flux) & (flux > 0)
                time = np.ascontiguousarray(time[valid], dtype=np.float64)
                flux = np.ascontiguousarray(flux[valid], dtype=np.float64)
                flux /= np.median(flux)

            period = ephem['period']
            dur_days = ephem['dur_hours'] / 24.0
            
            # Kepler time is BKJD (BJD - 2454833)
            t0 = ephem['t0']
            if t0 > 2450000:
                t0_bkjd = t0 - 2454833.0
            else:
                t0_bkjd = t0

            mask_in, mask_out = create_ephemeris_mask(time, t0_bkjd, period, dur_days, buffer_factor=2.0)
            n_in = int(np.sum(mask_in))
            n_out = int(np.sum(mask_out))

            if n_in < 4 or n_out < 10:
                continue

            comb = mask_in | mask_out
            sub_f = flux[comb]
            sub_m = mask_in[comb]
            sub_t = time[comb]

            # Compute WITHOUT shape correction
            res_raw = calculate_log_snras(sub_f, sub_m, time=sub_t, period=period, t0=t0_bkjd, apply_shape_correction=False)
            # Compute WITH shape correction
            res_corr = calculate_log_snras(sub_f, sub_m, time=sub_t, period=period, t0=t0_bkjd, apply_shape_correction=True)
            tier_corr = get_calibrated_tier(res_corr['penalty'], n_in, n_out)

            baselines = compute_baseline_metrics(sub_f, sub_m, duration_hours=ephem['dur_hours'])

            results.append({
                'Planet': planet_name,
                'Depth (%)': ephem['depth_ppm'] / 10000.0,
                'N_in': n_in,
                'Raw Penalty': res_raw['penalty'],
                'Corr Penalty': res_corr['penalty'],
                'Delta Penalty': res_raw['penalty'] - res_corr['penalty'],
                'Tier': tier_corr,
                'T-SNR': baselines['t_snr'],
                'L-SNRAS': res_corr['log_snras']
            })

            print(f"Target: {planet_name:18s} | Depth: {ephem['depth_ppm']/10000.:.2f}% | "
                  f"Raw P: {res_raw['penalty']:.4f} -> Corr P: {res_corr['penalty']:.4f} | "
                  f"Tier: {tier_corr:15s}")

        except Exception as err:
            # print(f"Error testing {planet_name}: {err}")
            pass

    if results:
        df = pd.DataFrame(results)
        print("\n" + "=" * 75)
        print("SUMMARY TABLE OF KEPLER TARGETS:")
        print("=" * 75)
        print(df.to_string(index=False))


if __name__ == '__main__':
    test_kepler_targets()
