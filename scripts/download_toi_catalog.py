"""
download_toi_catalog.py — تنزيل كتالوج TESS TOI من NASA / ExoFOP
"""

from __future__ import annotations
import sys, io, warnings, requests, pandas as pd, numpy as np, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
warnings.filterwarnings('ignore')

N_CP_TARGET = 200
N_FP_TARGET = 150
OUT_CSV = "data/toi_catalog_labeled.csv"

HEADERS = {
    "User-Agent": "Mozilla/5.0 (compatible; AcademicResearch/1.0)",
    "Accept": "text/csv,text/plain,*/*",
}


def try_astroquery():
    from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
    print("  [astroquery] Querying NasaExoplanetArchive...")
    tbl = NasaExoplanetArchive.query_criteria(
        table="toi",
        select="toi,tic_id,toi_disposition,pl_orbper,pl_trandep,pl_trandurh,pl_tranmid,st_tmag,ra,dec,sectors",
        where="toi_disposition in ('CP','FP') and pl_trandep is not null",
    )
    return tbl.to_pandas()


def try_exofop():
    url = "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv"
    print(f"  [ExoFOP] {url}")
    resp = requests.get(url, headers=HEADERS, timeout=60)
    resp.raise_for_status()
    from io import StringIO
    return pd.read_csv(StringIO(resp.text), comment='#')


def try_nexsci_api():
    url = (
        "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
        "?QUERY=SELECT+toi,tic_id,toi_disposition,pl_orbper,pl_trandep,"
        "pl_trandurh,pl_tranmid,st_tmag,ra,dec"
        "+FROM+toi"
        "+WHERE+toi_disposition+IN+('CP','FP')"
        "+AND+pl_trandep+IS+NOT+NULL"
        "&FORMAT=csv"
        "&LANG=ADQL"
    )
    print(f"  [TAP API] {url[:60]}...")
    resp = requests.get(url, headers=HEADERS, timeout=120)
    resp.raise_for_status()
    from io import StringIO
    return pd.read_csv(StringIO(resp.text), comment='#')


def download_catalog():
    print("=" * 64)
    print("  NASA TOI CATALOG DOWNLOAD")
    print("=" * 64)
    methods = [
        ("astroquery", try_astroquery),
        ("ExoFOP",     try_exofop),
        ("TAP API",    try_nexsci_api),
    ]
    for name, fn in methods:
        try:
            print(f"\n[Trying] {name}...")
            df = fn()
            if df is not None and len(df) > 10:
                print(f"  [OK] {len(df)} rows via {name}")
                return df
        except Exception as e:
            print(f"  [FAIL] {e}")
    raise RuntimeError("All download methods failed")


def process(df: pd.DataFrame) -> pd.DataFrame:
    print(f"\nRaw rows: {len(df)}")
    df.columns = [c.lower().strip() for c in df.columns]

    # ExoFOP uses 'tfopwg disposition', NASA uses 'toi_disposition'
    disp_col = next((c for c in ['toi_disposition', 'tfopwg disposition',
                                  'disposition', 'tfopwg_disp']
                     if c in df.columns), None)
    if not disp_col:
        raise ValueError(f"No disposition column. Have: {list(df.columns)}")

    # TIC ID column: ExoFOP uses 'tic id', NASA uses 'tic_id'
    tic_col = next((c for c in ['tic_id', 'tic id'] if c in df.columns), None)
    if tic_col:
        df = df.rename(columns={tic_col: 'tic_id'})

    # Orbital params: ExoFOP uses 'period (days)', 'duration (hours)', 'epoch (bjd)', 'depth (ppm)'
    col_map = {
        'period (days)':    'pl_orbper',
        'duration (hours)': 'pl_trandurh',
        'epoch (bjd)':      'pl_tranmid',
        'depth (ppm)':      'pl_trandep',
        'tess mag':         'st_tmag',
        'toi':              'toi',
    }
    for old, new in col_map.items():
        if old in df.columns and new not in df.columns:
            df = df.rename(columns={old: new})

    # Filter by disposition (CP = confirmed planet, FP = false positive)
    # ExoFOP values: 'CP', 'FP', 'PC', 'APC'
    valid = df[df[disp_col].isin(['CP', 'FP'])].copy()
    valid['toi_disposition'] = valid[disp_col]
    valid['label'] = (valid[disp_col] == 'CP').astype(int)

    # Filter: need depth and duration
    for col in ['pl_trandep', 'pl_trandurh', 'pl_orbper']:
        if col in valid.columns:
            valid = valid[valid[col].notna()]

    # Tmag filter: bright enough for good LC
    if 'st_tmag' in valid.columns:
        valid['st_tmag'] = pd.to_numeric(valid['st_tmag'], errors='coerce')
        valid = valid[valid['st_tmag'] < 14.0]

    n_cp = (valid['label']==1).sum()
    n_fp = (valid['label']==0).sum()
    print(f"  After filters — CP: {n_cp}  |  FP: {n_fp}")

    if n_cp == 0 or n_fp == 0:
        raise ValueError(f"Not enough samples after filtering: CP={n_cp}, FP={n_fp}")

    cp_sel = valid[valid['label']==1].sample(n=min(N_CP_TARGET, n_cp), random_state=42)
    fp_sel = valid[valid['label']==0].sample(n=min(N_FP_TARGET, n_fp), random_state=42)
    out = pd.concat([cp_sel, fp_sel]).sample(frac=1, random_state=42).reset_index(drop=True)
    print(f"Selected: CP={len(cp_sel)}, FP={len(fp_sel)}, Total={len(out)}")
    return out



def main():
    os.makedirs('data', exist_ok=True)
    df = download_catalog()
    df = process(df)
    df.to_csv(OUT_CSV, index=False)
    print(f"\n[SAVED] {len(df)} rows → {OUT_CSV}")
    print(df['label'].value_counts().rename({1:'CP', 0:'FP'}))


if __name__ == "__main__":
    main()
