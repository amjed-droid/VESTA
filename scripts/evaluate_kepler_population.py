"""
evaluate_kepler_population.py
==============================
خط أنابيب التقييم السكاني الشامل لـ 50 كوكباً من بعثة Kepler عبر 4 رتب عشرية
من عمق العبور (من 32,000 ppm إلى 12 ppm) إضافة إلى أهداف الشذوذ الفيزيائي.
"""

from __future__ import annotations
import sys, io, os, time, glob, warnings
import numpy as np
import pandas as pd
import lightkurve as lk
from sklearn.metrics import roc_auc_score

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
warnings.filterwarnings('ignore')

sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('src'))
from src.log_snras.core import calculate_log_snras, get_calibrated_tier, compute_baseline_metrics
from src.log_snras.masking import create_ephemeris_mask
from src.log_snras.stats import compute_auc, stratified_bootstrap_ci, delong_test

IN_CATALOG  = "data/kepler_population_catalog.csv"
OUT_RESULTS = "data/kepler_population_results.csv"
OUT_PERF    = "data/kepler_population_performance.csv"


def get_fits_files(fits_dir: str):
    """البحث عن كافة ملفات FITS في المجلد والمجلدات الفرعية مع تفضيل llc."""
    all_f = []
    for root, _, fnames in os.walk(fits_dir):
        for fn in fnames:
            if fn.endswith('_llc.fits') or fn.endswith('.fits'):
                all_f.append(os.path.join(root, fn))
    llc = [f for f in all_f if f.endswith('_llc.fits')]
    return sorted(list(set(llc if llc else all_f)))


def load_kepler_lightcurve_from_fits(fits_dir: str, max_files: int = 8):
    """تحميل ودمج ملفات FITS للهدف بدقة وسرعة."""
    files = get_fits_files(fits_dir)
    if not files:
        return None, None

    selected = files[:max_files]
    try:
        lcs = []
        for f in selected:
            try:
                lc = lk.read(f)
                lcs.append(lc)
            except Exception:
                continue
        if not lcs:
            return None, None

        if len(lcs) == 1:
            stitched = lcs[0]
        else:
            stitched = lk.LightCurveCollection(lcs).stitch()

        stitched = stitched.normalize().remove_outliers(sigma=5)
        t = np.array(stitched.time.value, dtype=float)
        f = np.array(stitched.flux.value, dtype=float)

        mask = np.isfinite(t) & np.isfinite(f) & (f > 0)
        return t[mask], f[mask]
    except Exception:
        return None, None


def load_artifact_lightcurve(cache_file: str):
    """تحميل بيانات الشذوذ الفلكي من ملف npz."""
    if not cache_file or not os.path.exists(cache_file):
        return None, None
    try:
        data = np.load(cache_file)
        t = data['time']
        f = data['flux']
        mask = np.isfinite(t) & np.isfinite(f) & (f > 0)
        return t[mask], f[mask]
    except Exception:
        return None, None


def main():
    print("=" * 75)
    print("LOG-SNRAS POPULATION-SCALE BENCHMARK EVALUATION (KEPLER N=50 + ARTIFACTS)")
    print("=" * 75)

    if not os.path.exists(IN_CATALOG):
        print(f"ERROR: {IN_CATALOG} not found.")
        sys.exit(1)

    df_cat = pd.read_csv(IN_CATALOG)
    print(f"Loaded {len(df_cat)} systems from catalog.\n")

    results = []
    t_start = time.time()

    for idx, row in df_cat.iterrows():
        rank = row['rank']
        name = row['target_name']
        ttype = row['target_type']
        lbl = int(row['label'])
        period = float(row['period'])
        t0 = float(row['t0_bkjd'])
        dur_h = float(row['duration_hours'])
        dur_d = dur_h / 24.0
        depth_ppm = float(row['depth_ppm'])
        stratum = row['depth_stratum']
        fits_dir = str(row.get('fits_dir', ''))
        cache_file = str(row.get('cache_file', ''))

        lbl_str = "PLANET" if lbl == 1 else "ARTIFACT"
        print(f"[{idx+1:2d}/{len(df_cat)}] Rank {rank:3d} | {name:<24} | {lbl_str:<8} | {stratum:<32} ... ", end='', flush=True)

        t0_sub = time.time()
        if ttype == 'confirmed_planet':
            time_arr, flux_arr = load_kepler_lightcurve_from_fits(fits_dir, max_files=8)
        else:
            time_arr, flux_arr = load_artifact_lightcurve(cache_file)

        if time_arr is None or len(time_arr) < 50:
            print(f"FAIL (no data) [{time.time()-t0_sub:.1f}s]")
            continue

        # Ephemeris masking
        try:
            transit_mask, baseline_mask = create_ephemeris_mask(
                time_arr, t0, period, dur_d, buffer_factor=1.5
            )
        except Exception:
            transit_mask = np.zeros(len(time_arr), dtype=bool)
            baseline_mask = np.ones(len(time_arr), dtype=bool)

        n_in = int(np.sum(transit_mask))
        n_out = int(np.sum(baseline_mask))

        if n_in < 4 or n_out < 10:
            print(f"SKIP (n_in={n_in}) [{time.time()-t0_sub:.1f}s]")
            continue

        # Compute Log-SNRAS
        try:
            res_dict = calculate_log_snras(
                flux_arr, transit_mask,
                time=time_arr,
                period=period if period > 0 else None,
                t0=t0 if t0 > 0 else None,
                apply_shape_correction=True
            )
            penalty = float(res_dict['penalty'])
            psi = float(res_dict['psi'])
            log_snras = float(res_dict['log_snras'])
            measured_depth = float(res_dict['depth'])

            tier_str = get_calibrated_tier(penalty, n_in=max(n_in, 5), n_out=max(n_out, 50))
            tier_num = 1 if 'Tier 1' in tier_str else (2 if 'Tier 2' in tier_str else 3)

            baselines = compute_baseline_metrics(
                flux_arr, transit_mask, duration_hours=dur_h
            )
            t_snr = float(baselines['t_snr'])
            r_snr = float(baselines['r_snr'])
            p_snr = float(baselines['p_snr'])
            b_snr = float(baselines['b_snr'])

            elapsed = time.time() - t0_sub
            print(f"OK | P={penalty:.4f} | {tier_str} | SNR_T={t_snr:.1f} | L_SNRAS={log_snras:.1f} [{elapsed:.1f}s]")

            results.append({
                'rank': rank,
                'target_name': name,
                'star_name': row['star_name'],
                'target_type': ttype,
                'label': lbl,
                'period': period,
                't0_bkjd': t0,
                'duration_hours': dur_h,
                'literature_depth_ppm': depth_ppm,
                'measured_depth': measured_depth,
                'depth_stratum': stratum,
                'n_in': n_in,
                'n_out': n_out,
                'n_total_points': len(time_arr),
                'penalty_corr': penalty,
                'psi': psi,
                'log_snras': log_snras,
                'tier': tier_str,
                'tier_num': tier_num,
                't_snr': t_snr,
                'r_snr': r_snr,
                'p_snr': p_snr,
                'b_snr': b_snr,
                'status': 'ok'
            })
        except Exception as e:
            print(f"ERR: {str(e)[:30]} [{time.time()-t0_sub:.1f}s]")

    # Save results
    res_df = pd.DataFrame(results)
    res_df.to_csv(OUT_RESULTS, index=False)
    total_time = time.time() - t_start

    print("\n" + "=" * 75)
    print(f"EVALUATION COMPLETE: {len(res_df)} / {len(df_cat)} targets processed in {total_time:.1f}s")
    print(f"Output saved to: {OUT_RESULTS}")
    print("=" * 75)

    ok_planets = res_df[res_df['label'] == 1]
    ok_artifacts = res_df[res_df['label'] == 0]
    n_p = len(ok_planets)
    n_a = len(ok_artifacts)

    print(f"\nPopulation Statistics (Planets: {n_p}, Artifacts: {n_a}):")
    print(f"  Planets mean Penalty:   {ok_planets['penalty_corr'].mean():.4f} +/- {ok_planets['penalty_corr'].std():.4f}")
    print(f"  Artifacts mean Penalty: {ok_artifacts['penalty_corr'].mean():.4f} +/- {ok_artifacts['penalty_corr'].std():.4f}")

    print("\nPlanets Performance by Depth Stratum:")
    for strat, grp in ok_planets.groupby('depth_stratum'):
        clean_pct = (grp['tier_num'] == 1).mean() * 100
        print(f"  {strat:<36}: N={len(grp):2d} | Mean P={grp['penalty_corr'].mean():.4f} | Tier 1 (Clean): {clean_pct:.1f}%")

    veto_pct = (ok_artifacts['tier_num'] >= 2).mean() * 100
    print(f"\nArtifacts Rejection Rate: {veto_pct:.1f}% placed in Review/Veto Tiers")

    # AUC & Multi-metric benchmarking
    y_true = res_df['label'].values
    pen_score = -res_df['penalty_corr'].values
    pen_auc = roc_auc_score(y_true, pen_score)

    metrics_to_eval = [
        ('Penalty (Shape-Corrected)', pen_score),
        ('Log-SNRAS Composite',       res_df['log_snras'].values),
        ('Traditional SNR',           res_df['t_snr'].values),
        ('Robust SNR (MAD)',          res_df['r_snr'].values),
        ('Pont SNR (2006)',           res_df['p_snr'].values),
        ('BLS SNR Proxy',             res_df['b_snr'].values),
        ('Inverse Depth (1/delta)',   1.0 / (res_df['measured_depth'].values + 1e-12))
    ]

    perf_records = []
    print("\n" + "-" * 88)
    print(f"{'Metric':<28} | {'AUC':<6} | {'95% CI (Bootstrap)':<22} | {'DeLong vs Pen (p)':<18}")
    print("-" * 88)

    for name, scores in metrics_to_eval:
        auc_val = roc_auc_score(y_true, scores)
        _, (ci_lo, ci_hi), _ = stratified_bootstrap_ci(y_true, scores, n_boot=2000, alpha=0.05)

        if name == 'Penalty (Shape-Corrected)':
            dz, dp = 0.0, 1.0
            delta_auc = 0.0
        else:
            try:
                dz, dp = delong_test(y_true, pen_score, scores)
                delta_auc = pen_auc - auc_val
            except Exception:
                dz, dp, delta_auc = np.nan, np.nan, np.nan

        print(f"{name:<28} | {auc_val:6.4f} | [{ci_lo:6.4f}, {ci_hi:6.4f}]       | z={dz:4.2f}, p={dp:.4f}")

        perf_records.append({
            'Metric': name,
            'AUC': auc_val,
            'CI_Lower': ci_lo,
            'CI_Upper': ci_hi,
            'Delta_AUC_vs_Penalty': delta_auc,
            'DeLong_z': dz,
            'DeLong_p': dp
        })

    perf_df = pd.DataFrame(perf_records)
    perf_df.to_csv(OUT_PERF, index=False)
    print("-" * 88)
    print(f"Performance summary saved to: {OUT_PERF}")


if __name__ == '__main__':
    main()
