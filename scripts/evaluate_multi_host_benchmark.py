"""
Master Multi-Host Benchmark Evaluation Pipeline
================================================
Downloads/loads real Kepler and TESS light curves for a diverse, balanced
multi-host sample of confirmed planets and non-planetary artifacts.
Applies uniform masking, shape-corrected dispersion, and calibrated tiers,
and produces comprehensive evaluation metrics and publication figures.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import lightkurve as lk

# Add src and current directory to python path
sys.path.insert(0, os.path.abspath('src'))
sys.path.insert(0, os.path.abspath('.'))
from log_snras import (
    calculate_log_snras,
    compute_baseline_metrics,
    get_calibrated_tier
)
from log_snras.masking import create_ephemeris_mask
from log_snras.stats import (
    compute_auc,
    stratified_bootstrap_ci,
    delong_test,
    permutation_test
)
from scripts.create_curated_catalog import BENCHMARK_TARGETS


def download_and_extract_target(target_info, cache_dir='data/raw_lc'):
    """
    Downloads or retrieves cached light curve for a benchmark target.
    """
    os.makedirs(cache_dir, exist_ok=True)
    target_id = target_info['target_id']
    sec_qtr = target_info['sector_or_quarter']
    mission = target_info['mission']
    author = target_info.get('author', 'SPOC' if mission == 'TESS' else 'Kepler')
    
    cache_file = os.path.join(
        cache_dir,
        f"{target_id.replace(' ', '_')}_{mission}_{sec_qtr}_{author}.npz"
    )

    if os.path.exists(cache_file):
        data = np.load(cache_file)
        return data['time'], data['flux']

    print(f"Fetching from MAST: {target_id} ({mission} {sec_qtr}, author={author})...")
    try:
        if mission == 'TESS':
            search = lk.search_lightcurve(target_id, sector=sec_qtr, author=author)
        else:
            search = lk.search_lightcurve(target_id, quarter=sec_qtr, author=author)

        if len(search) == 0:
            # Fallback search without author
            search = lk.search_lightcurve(target_id)
            if len(search) == 0:
                raise ValueError(f"Target {target_id} not found in MAST")

        lc = search[0].download()
        time = np.array(lc.time.value, dtype=float)
        
        # Use PDCSAP_FLUX if available, fallback to SAP_FLUX
        if hasattr(lc, 'pdcsap_flux') and lc.pdcsap_flux is not None:
            raw_flux = np.array(lc.pdcsap_flux.value, dtype=float)
        elif hasattr(lc, 'flux') and lc.flux is not None:
            raw_flux = np.array(lc.flux.value, dtype=float)
        else:
            raise KeyError(f"No valid flux attribute in light curve for {target_id}")

        valid = np.isfinite(time) & np.isfinite(raw_flux) & (raw_flux > 0)
        time = time[valid]
        raw_flux = raw_flux[valid]

        # Median normalize
        flux = raw_flux / np.median(raw_flux)

        np.savez(cache_file, time=time, flux=flux)
        return time, flux
    except Exception as e:
        print(f"Error fetching {target_id}: {e}")
        return None, None


def run_benchmark():
    """Runs the entire benchmark evaluation across all targets."""
    print("=" * 70)
    print("LOG-SNRAS MULTI-HOST BENCHMARK EVALUATION")
    print("=" * 70)

    rows = []
    
    for item in BENCHMARK_TARGETS:
        t_id = item['target_id']
        host = item['host_name']
        print(f"\nProcessing: {host} ({t_id})...")

        time, flux = download_and_extract_target(item)
        if time is None or len(time) == 0:
            print(f"Skipping {host}: data unavailable.")
            continue

        # Uniform ephemeris-based mask
        t0 = item['t0']
        period = item['period']
        duration = item['duration_days']
        
        transit_mask, baseline_mask = create_ephemeris_mask(
            time, t0, period, duration, buffer_factor=2.0
        )

        n_in = int(np.sum(transit_mask))
        n_out = int(np.sum(baseline_mask))

        if n_in < 4 or n_out < 10:
            print(f"Insufficient points: N_in={n_in}, N_out={n_out}")
            continue

        # Compute baselines
        combined_mask = transit_mask | baseline_mask
        sub_flux = flux[combined_mask]
        sub_mask = transit_mask[combined_mask]

        sub_time = time[combined_mask]

        # Compute baselines
        baselines = compute_baseline_metrics(sub_flux, sub_mask, duration_hours=duration * 24.0)

        # Compute Log-SNRAS without shape correction
        res_raw = calculate_log_snras(sub_flux, sub_mask, time=sub_time, period=period, t0=t0, apply_shape_correction=False)
        tier_raw = get_calibrated_tier(res_raw['penalty'], n_in, n_out)

        # Compute Log-SNRAS WITH shape correction
        res_corr = calculate_log_snras(sub_flux, sub_mask, time=sub_time, period=period, t0=t0, apply_shape_correction=True)
        tier_corr = get_calibrated_tier(res_corr['penalty'], n_in, n_out)

        row = {
            'target_id': t_id,
            'host_name': host,
            'mission': item['mission'],
            'sector_or_quarter': item['sector_or_quarter'],
            'label': item['label'],
            'label_name': item['label_name'],
            'depth_catalog_ppm': item['depth_ppm'],
            'depth_measured_ppm': baselines['depth'] * 1e6,
            'n_in': n_in,
            'n_out': n_out,
            't_snr': baselines['t_snr'],
            'r_snr': baselines['r_snr'],
            'p_snr': baselines['p_snr'],
            'b_snr': baselines['b_snr'],
            'psi_raw': res_raw['psi'],
            'penalty_raw': res_raw['penalty'],
            'tier_raw': tier_raw,
            'l_snras_raw': res_raw['log_snras'],
            'psi_corr': res_corr['psi'],
            'penalty_corr': res_corr['penalty'],
            'tier_corr': tier_corr,
            'l_snras_corr': res_corr['log_snras'],
            'literature_ref': item['literature_ref']
        }
        rows.append(row)
        print(f"  Result: T-SNR={baselines['t_snr']:.2f}, L-SNRAS={res_corr['log_snras']:.2f}, "
              f"Penalty={res_corr['penalty']:.4f}, Tier={tier_corr}")

    df_results = pd.DataFrame(rows)
    df_results.to_csv('data/multi_host_evaluation_results.csv', index=False)
    print("\nSaved evaluation results to data/multi_host_evaluation_results.csv")

    # -------------------------------------------------------------------------
    # Statistical Analysis
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("BENCHMARK STATISTICAL PERFORMANCE SUMMARY")
    print("=" * 70)

    y_true = df_results['label'].values.astype(bool)

    metrics_to_test = {
        'Penalty (Shape-Corrected)': -df_results['penalty_corr'].values,
        'Penalty (Uncorrected)': -df_results['penalty_raw'].values,
        'Traditional SNR': df_results['t_snr'].values,
        'Robust SNR (MAD)': df_results['r_snr'].values,
        'Pont SNR (2006)': df_results['p_snr'].values,
        'BLS SNR Proxy': df_results['b_snr'].values,
        'Log-SNRAS Composite': df_results['l_snras_corr'].values,
        'Inverse Depth (ppm)': -df_results['depth_measured_ppm'].values
    }

    perf_records = []
    pen_scores = metrics_to_test['Penalty (Shape-Corrected)']

    for name, scores in metrics_to_test.items():
        auc, (ci_lo, ci_hi), _ = stratified_bootstrap_ci(y_true, scores, n_boot=2000)
        
        if name != 'Penalty (Shape-Corrected)':
            diff_auc, z_score, p_val = delong_test(y_true, pen_scores, scores)
        else:
            diff_auc, z_score, p_val = 0.0, 0.0, 1.0

        perf_records.append({
            'Metric': name,
            'AUC': auc,
            '95% CI Lower': ci_lo,
            '95% CI Upper': ci_hi,
            'Delta AUC vs Penalty': diff_auc,
            'DeLong z': z_score,
            'DeLong p-value': p_val
        })
        print(f"{name:26s} | AUC = {auc:.3f} [{ci_lo:.3f}, {ci_hi:.3f}] | z = {z_score:5.2f}, p = {p_val:.4f}")

    df_perf = pd.DataFrame(perf_records)
    df_perf.to_csv('data/multi_host_auc_performance.csv', index=False)
    print("\nSaved statistical summary to data/multi_host_auc_performance.csv")

    # -------------------------------------------------------------------------
    # Multi-Host Leave-One-Host-Out Cross Validation
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("LEAVE-ONE-HOST-OUT CROSS-VALIDATION (LOHO)")
    print("=" * 70)
    
    unique_hosts = df_results['host_name'].unique()
    loho_records = []
    
    for host in unique_hosts:
        subset = df_results[df_results['host_name'] != host]
        y_sub = subset['label'].values.astype(bool)
        if len(np.unique(y_sub)) < 2:
            continue
        
        auc_loho = compute_auc(y_sub, -subset['penalty_corr'].values)
        loho_records.append({'Excluded Host': host, 'Remaining N': len(subset), 'AUC': auc_loho})
        print(f"Excluding {host:20s} | Remaining N={len(subset):2d} | AUC = {auc_loho:.3f}")

    df_loho = pd.DataFrame(loho_records)
    print(f"\nMean LOHO AUC: {df_loho['AUC'].mean():.3f} +/- {df_loho['AUC'].std():.3f}")

    return df_results, df_perf


if __name__ == '__main__':
    run_benchmark()
