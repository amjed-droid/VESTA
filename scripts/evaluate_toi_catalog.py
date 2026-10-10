"""
evaluate_toi_catalog.py — تشغيل VESTA على كتالوج TESS TOI الكامل
======================================================================
يقرأ: data/toi_catalog_labeled.csv
يُنزّل ويُخزن منحنيات الضوء مع حد أقصى للمهلة (Timeout 30s)
يصحح الـ epoch من BJD إلى BTJD (BJD - 2457000)
يحسب: VESTA, penalty_corr, Tiers, T-SNR, R-SNR, P-SNR, B-SNR
يكتب: data/toi_evaluation_results.csv
"""

from __future__ import annotations
import sys, io, os, time, warnings, argparse, concurrent.futures
import numpy as np
import pandas as pd

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
warnings.filterwarnings('ignore')

# Add project root and src to path
sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('src'))
from src.vesta.core import calculate_vesta, get_calibrated_tier, compute_baseline_metrics
from src.vesta.masking import create_ephemeris_mask

IN_CSV    = "data/toi_catalog_labeled.csv"
OUT_CSV   = "data/toi_evaluation_results.csv"
CACHE_DIR = "data/raw_lc/toi"


def _fetch_from_mast(tic_id: int):
    """جلب منحنى الضوء من MAST عبر lightkurve."""
    import lightkurve as lk
    for cadence in ['2min', '10min', '30min']:
        try:
            search = lk.search_lightcurve(
                f"TIC {tic_id}", mission="TESS",
                author="SPOC", cadence=cadence
            )
            if len(search) == 0:
                continue
            lcc = search[0].download(flux_column='pdcsap_flux')
            if lcc is None:
                continue
            lc = lcc.normalize().remove_outliers(sigma=5)
            time_arr = np.array(lc.time.value, dtype=float)
            flux_arr = np.array(lc.flux.value, dtype=float)
            if hasattr(lc, 'flux_err') and lc.flux_err is not None:
                ferr_arr = np.array(lc.flux_err.value, dtype=float)
            else:
                ferr_arr = np.full_like(flux_arr, float(np.nanstd(flux_arr) * 0.01))

            mask = np.isfinite(time_arr) & np.isfinite(flux_arr) & np.isfinite(ferr_arr) & (flux_arr > 0)
            t_c, f_c, e_c = time_arr[mask], flux_arr[mask], ferr_arr[mask]
            if len(t_c) >= 50:
                return t_c, f_c, e_c, cadence
        except Exception:
            continue
    return None, None, None, None


def get_tess_lc(tic_id: int, timeout_sec: int = 30):
    """جلب منحنى الضوء مع التخزين المؤقت وحد زمني صارم."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    cache_path = os.path.join(CACHE_DIR, f"TIC_{tic_id}.npz")

    if os.path.exists(cache_path):
        try:
            d = np.load(cache_path)
            return d['time'], d['flux'], d['flux_err'], str(d['cadence'])
        except Exception:
            pass

    # Fetch with timeout
    with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
        future = executor.submit(_fetch_from_mast, tic_id)
        try:
            t_c, f_c, e_c, cadence = future.result(timeout=timeout_sec)
            if t_c is not None and len(t_c) >= 50:
                np.savez_compressed(
                    cache_path,
                    time=t_c, flux=f_c, flux_err=e_c, cadence=cadence
                )
                return t_c, f_c, e_c, cadence
        except concurrent.futures.TimeoutError:
            return None, None, None, 'timeout'
        except Exception:
            return None, None, None, None

    return None, None, None, None


def evaluate_one(row: pd.Series) -> dict:
    """تقييم TOI واحد وحساب كافة المقاييس الإحصائية بدقة."""
    tic_id = int(row.get('tic_id', 0))
    result = {
        'toi':             row.get('toi', '?'),
        'tic_id':          tic_id,
        'toi_disposition': row.get('toi_disposition', '?'),
        'label':           int(row.get('label', -1)),
        'penalty_corr':    np.nan,
        'tier':            'None',
        'tier_num':        -1,
        'vesta_score':     np.nan,
        'psi':             np.nan,
        't_snr':           np.nan,
        'r_snr':           np.nan,
        'p_snr':           np.nan,
        'b_snr':           np.nan,
        'cadence':         None,
        'n_transit_pts':   0,
        'n_oot_pts':       0,
        'status':          'pending',
    }

    if tic_id == 0:
        result['status'] = 'no_tic'
        return result

    period = float(row.get('pl_orbper',   0) or 0)
    t0_raw = float(row.get('pl_tranmid',  0) or 0)
    dur_h  = float(row.get('pl_trandurh', 0) or 0)
    dur_d  = dur_h / 24.0

    # Convert BJD to BTJD if necessary
    t0 = (t0_raw - 2457000.0) if t0_raw > 2400000.0 else t0_raw

    time_arr, flux_arr, ferr_arr, cad = get_tess_lc(tic_id, timeout_sec=25)
    if cad == 'timeout':
        result['status'] = 'timeout'
        return result
    if time_arr is None or len(time_arr) < 50:
        result['status'] = 'no_lc'
        return result

    result['cadence'] = cad

    flux_med = np.nanmedian(flux_arr)
    if flux_med <= 0 or np.isnan(flux_med):
        result['status'] = 'bad_flux'
        return result
    f_norm = flux_arr / flux_med

    # Transit and baseline masks
    try:
        if period > 0 and t0 > 0 and dur_d > 0:
            transit_mask, baseline_mask = create_ephemeris_mask(
                time_arr, t0, period, dur_d, buffer_factor=1.5
            )
        else:
            transit_mask  = np.zeros(len(time_arr), dtype=bool)
            baseline_mask = np.ones(len(time_arr),  dtype=bool)
    except Exception:
        transit_mask  = np.zeros(len(time_arr), dtype=bool)
        baseline_mask = np.ones(len(time_arr),  dtype=bool)

    n_in  = int(np.sum(transit_mask))
    n_out = int(np.sum(baseline_mask))
    result['n_transit_pts'] = n_in
    result['n_oot_pts']     = n_out

    if n_in < 4 or n_out < 10:
        result['status'] = 'insufficient_pts'
        return result

    # Compute VESTA
    try:
        res_dict = calculate_vesta(
            f_norm, transit_mask,
            time=time_arr,
            period=period if period > 0 else None,
            t0=t0 if t0 > 0 else None,
            apply_shape_correction=True,
        )
        penalty  = res_dict['penalty']
        tier_str = get_calibrated_tier(penalty, n_in=max(n_in, 5), n_out=max(n_out, 50))
        tier_num = 1 if 'Tier 1' in tier_str else (2 if 'Tier 2' in tier_str else 3)

        result['penalty_corr'] = float(penalty)
        result['tier']         = str(tier_str)
        result['tier_num']     = int(tier_num)
        result['vesta_score']  = float(res_dict['vesta_score'])
        result['psi']          = float(res_dict['psi'])

        # Compute comparator baselines
        baselines = compute_baseline_metrics(
            f_norm, transit_mask, duration_hours=dur_h if dur_h > 0 else 2.5
        )
        result['t_snr'] = float(baselines.get('t_snr', np.nan))
        result['r_snr'] = float(baselines.get('r_snr', np.nan))
        result['p_snr'] = float(baselines.get('p_snr', np.nan))
        result['b_snr'] = float(baselines.get('b_snr', np.nan))

        result['status'] = 'ok'
    except Exception as e:
        result['status'] = f'err:{str(e)[:40]}'

    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--max',    type=int, default=50,   help='Max TOIs to process (default 50, balanced CP/FP)')
    parser.add_argument('--resume', action='store_true',    help='Resume from existing output')
    args = parser.parse_args()

    if not os.path.exists(IN_CSV):
        print(f"ERROR: {IN_CSV} not found. Run download_toi_catalog.py first.")
        sys.exit(1)

    df = pd.read_csv(IN_CSV)
    print(f"Loaded {len(df)} TOIs from {IN_CSV}")

    # Resume mode
    processed_tics = set()
    existing_rows = []
    if args.resume and os.path.exists(OUT_CSV):
        done = pd.read_csv(OUT_CSV)
        # Keep valid rows
        done_ok = done[done['status'].isin(['ok', 'no_lc', 'insufficient_pts', 'timeout'])]
        processed_tics = set(done_ok['tic_id'].astype(int).tolist())
        existing_rows = done_ok.to_dict('records')
        print(f"Resuming: {len(processed_tics)} already processed")
        df = df[~df['tic_id'].astype(int).isin(processed_tics)].reset_index(drop=True)

    if args.max > 0 and len(df) > args.max:
        n_each = args.max // 2
        cp_df = df[df['label'] == 1].head(n_each)
        fp_df = df[df['label'] == 0].head(n_each)
        df = pd.concat([cp_df, fp_df]).sample(frac=1.0, random_state=42).reset_index(drop=True)

    print(f"Processing {len(df)} TOIs (CP: {(df['label']==1).sum()}, FP: {(df['label']==0).sum()})...\n")

    results = list(existing_rows)
    n_ok = sum(1 for r in results if r.get('status') == 'ok')
    n_fail = len(results) - n_ok

    for i, row in df.iterrows():
        tic = int(row.get('tic_id', 0))
        toi = row.get('toi', '?')
        lbl = int(row.get('label', -1))
        lbl_str = 'CP' if lbl == 1 else 'FP'

        print(f"[{i+1:3d}/{len(df)}] TIC {tic:10d}  TOI {str(toi):<8} {lbl_str}  ", end='', flush=True)

        t0_start = time.time()
        res = evaluate_one(row)
        elapsed = time.time() - t0_start

        p = res.get('penalty_corr', np.nan)
        tier = res.get('tier', '?')
        status = res.get('status', '?')

        if status == 'ok':
            print(f"P={p:.4f}  {tier}  [{elapsed:.1f}s]")
            n_ok += 1
        else:
            print(f"SKIP [{status}]  [{elapsed:.1f}s]")
            n_fail += 1

        results.append(res)

        # Checkpoint every 5 TOIs
        if (i + 1) % 5 == 0 or (i + 1) == len(df):
            pd.DataFrame(results).to_csv(OUT_CSV, index=False)

    final = pd.DataFrame(results)
    final.to_csv(OUT_CSV, index=False)

    print(f"\n{'='*60}")
    print(f"  DONE: {n_ok} OK, {n_fail} skipped/failed, Total={len(final)}")
    print(f"  Output: {OUT_CSV}")

    ok = final[final['status'] == 'ok']
    if len(ok) > 0:
        n_cp = (ok['label'] == 1).sum()
        n_fp = (ok['label'] == 0).sum()
        print(f"\nResults on successful targets (N={len(ok)}: CP={n_cp}, FP={n_fp}):")
        print(f"  CP mean penalty: {ok[ok['label']==1]['penalty_corr'].mean():.4f}")
        print(f"  FP mean penalty: {ok[ok['label']==0]['penalty_corr'].mean():.4f}")

        if n_cp >= 2 and n_fp >= 2:
            from sklearn.metrics import roc_auc_score
            for col, name in [
                ('penalty_corr', 'VESTA Penalty (Inv)'),
                ('vesta_score',  'VESTA Metric'),
                ('t_snr',        'Traditional SNR'),
                ('r_snr',        'Robust SNR (MAD)'),
                ('p_snr',        'Pont SNR'),
                ('b_snr',        'BLS SNR proxy'),
            ]:
                if col in ok.columns:
                    vals = ok[col].dropna()
                    valid_idx = vals.index
                    y_true = ok.loc[valid_idx, 'label']
                    if col == 'penalty_corr':
                        score = -ok.loc[valid_idx, col]
                    else:
                        score = ok.loc[valid_idx, col]
                    if len(np.unique(y_true)) > 1:
                        try:
                            auc = roc_auc_score(y_true, score)
                            print(f"  AUC {name:<26}: {auc:.4f}")
                        except Exception:
                            pass


if __name__ == "__main__":
    main()
