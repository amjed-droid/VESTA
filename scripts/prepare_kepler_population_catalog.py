"""
prepare_kepler_population_catalog.py
====================================
يجمع المعاملات المدارية الدقيقة (koi_period, koi_time0bk, koi_duration, koi_depth)
لكافة الكواكب الـ 50 من أرشيف ناسا ويطابقها مع مجلدات FITS المحلية في D:/nrebig/kepler.
ويضيف أهداف الشذوذ الفيزيائي الـ 7 المعتمدة في الأدبيات (المخزنة محلياً في data/raw_lc).
"""

import os, sys, glob
sys.path.insert(0, os.path.abspath('.'))
import numpy as np
import pandas as pd
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
from scripts.create_curated_catalog import BENCHMARK_TARGETS

BASE_DIR  = "D:/nrebig/kepler"
LOCAL_CSV = os.path.join(BASE_DIR, "targets_en.csv")
OUT_CSV   = "data/kepler_population_catalog.csv"

def build_catalog():
    print("Loading local 50-planet list...")
    df_local = pd.read_csv(LOCAL_CSV)
    print(f"Loaded {len(df_local)} planets from {LOCAL_CSV}")

    # Query NASA KOI cumulative table
    print("Querying NASA Exoplanet Archive KOI cumulative table...")
    koi_tab = NasaExoplanetArchive.query_criteria(
        table='cumulative',
        select='kepid,kepoi_name,kepler_name,koi_period,koi_time0bk,koi_duration,koi_depth'
    ).to_pandas()
    print(f"Retrieved {len(koi_tab)} KOIs")

    # Query PS table as backup
    print("Querying NASA PS table for any missing targets...")
    ps_tab = NasaExoplanetArchive.query_criteria(
        table='pscomppars',
        select='pl_name,hostname,pl_orbper,pl_tranmid,pl_trandur,pl_trandep'
    ).to_pandas()

    target_dirs = os.listdir(os.path.join(BASE_DIR, "lightcurves"))
    records = []

    for _, row in df_local.iterrows():
        rank = int(row['rank'])
        planet = row['planet'].strip()
        star = row['star'].strip()
        depth_pct = float(row['transit_depth_percent'])
        depth_ppm = depth_pct * 10000.0

        prefix = f"{rank:02d}_"
        match_dirs = [d for d in target_dirs if d.startswith(prefix)]
        if not match_dirs:
            continue
        target_dir = os.path.join(BASE_DIR, "lightcurves", match_dirs[0])

        koi_match = koi_tab[koi_tab['kepler_name'] == planet]
        if len(koi_match) > 0:
            k_row = koi_match.iloc[0]
            kepid = int(k_row['kepid'])
            t0 = float(k_row['koi_time0bk'])
            period = float(k_row['koi_period'])
            dur_hours = float(k_row['koi_duration'])
            koi_name = str(k_row['kepoi_name'])
        else:
            ps_match = ps_tab[ps_tab['pl_name'] == planet]
            if len(ps_match) > 0:
                p_row = ps_match.iloc[0]
                kepid = 0
                period = float(p_row['pl_orbper'])
                dur_hours = float(p_row['pl_trandur'])
                t0_raw = float(p_row['pl_tranmid'])
                t0 = (t0_raw - 2454833.0) if t0_raw > 2400000.0 else t0_raw
                koi_name = planet
            else:
                continue

        if depth_ppm >= 10000.0:
            stratum = 'Ultra-Deep (>10,000 ppm)'
        elif depth_ppm >= 1000.0:
            stratum = 'Moderate (1,000-10,000 ppm)'
        elif depth_ppm >= 100.0:
            stratum = 'Shallow Super-Earth (100-1,000 ppm)'
        else:
            stratum = 'Ultra-Shallow Sub-Earth (<100 ppm)'

        records.append({
            'target_type': 'confirmed_planet',
            'label': 1,
            'rank': rank,
            'target_name': planet,
            'star_name': star,
            'kepid': kepid,
            'koi_name': koi_name,
            'period': period,
            't0_bkjd': t0,
            'duration_hours': dur_hours,
            'depth_ppm': depth_ppm,
            'depth_stratum': stratum,
            'fits_dir': target_dir,
            'cache_file': ''
        })

    # Add 7 verified literature artifacts from BENCHMARK_TARGETS
    artifact_id = 101
    for b in BENCHMARK_TARGETS:
        if b['label'] == 0:
            t_id = b['target_id']
            mission = b['mission']
            sec_qtr = b['sector_or_quarter']
            author = b.get('author', 'SPOC' if mission == 'TESS' else 'Kepler')
            cache_file = os.path.join(
                'data/raw_lc',
                f"{t_id.replace(' ', '_')}_{mission}_{sec_qtr}_{author}.npz"
            )
            records.append({
                'target_type': 'artifact',
                'label': 0,
                'rank': artifact_id,
                'target_name': b['host_name'],
                'star_name': b['target_id'],
                'kepid': 0,
                'koi_name': b['label_name'],
                'period': float(b['period']),
                't0_bkjd': float(b['t0']),
                'duration_hours': float(b['duration_days']) * 24.0,
                'depth_ppm': float(b['depth_ppm']),
                'depth_stratum': 'Astrophysical Artifact',
                'fits_dir': '',
                'cache_file': cache_file
            })
            artifact_id += 1

    out_df = pd.DataFrame(records)
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    out_df.to_csv(OUT_CSV, index=False)
    print(f"\nWrote {len(out_df)} total systems to {OUT_CSV}")
    print(out_df['depth_stratum'].value_counts())

if __name__ == '__main__':
    build_catalog()
