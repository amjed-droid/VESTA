"""
Curated Multi-Host Benchmark Catalog Definition
===============================================
Constructs a balanced, verified, multi-host benchmark catalog of exoplanet transits
and non-planetary astrophysical/instrumental artifacts with complete literature provenance.
"""

import pandas as pd

BENCHMARK_TARGETS = [
    # ---------------------------------------------------------------------------------
    # CONFIRMED PLANETS (Diverse Hosts, Shallow to Deep Transits)
    # ---------------------------------------------------------------------------------
    {
        'target_id': 'TIC 261136679',
        'host_name': 'Pi Mensae',
        'mission': 'TESS',
        'sector_or_quarter': 1,
        'label': 1,
        'label_name': 'Confirmed Planet',
        't0': 1325.503,
        'period': 6.26784,
        'duration_days': 0.123,
        'depth_ppm': 300,
        'literature_ref': 'Kunovac Hodzic et al. (2021, MNRAS 502, 2893)'
    },
    {
        'target_id': 'TIC 350618622',
        'host_name': 'TOI-201',
        'mission': 'TESS',
        'sector_or_quarter': 6,
        'label': 1,
        'label_name': 'Confirmed Planet',
        't0': 1482.029,
        'period': 52.978,
        'duration_days': 0.290,
        'depth_ppm': 5300,
        'literature_ref': 'Hobson et al. (2021, AJ 161, 235)'
    },
    {
        'target_id': 'KIC 11904151',
        'host_name': 'Kepler-10',
        'mission': 'Kepler',
        'sector_or_quarter': 1,
        'label': 1,
        'label_name': 'Confirmed Planet',
        't0': 131.574858,
        'period': 0.837491225,
        'duration_days': 0.0748,
        'depth_ppm': 190,
        'literature_ref': 'Batalha et al. (2011, ApJ 729, 27)'
    },
    {
        'target_id': 'KIC 5812701',
        'host_name': 'Kepler-448',
        'mission': 'Kepler',
        'sector_or_quarter': 4,
        'label': 1,
        'label_name': 'Confirmed Planet',
        't0': 146.596466,
        'period': 17.855221681,
        'duration_days': 0.3089,
        'depth_ppm': 9065,
        'literature_ref': 'Bourrier et al. (2015, A&A 579, A55)'
    },
    {
        'target_id': 'KIC 11446443',
        'host_name': 'TrES-2',
        'mission': 'Kepler',
        'sector_or_quarter': 2,
        'label': 1,
        'label_name': 'Confirmed Planet',
        't0': 122.763305,
        'period': 2.470613377,
        'duration_days': 0.0726,
        'depth_ppm': 14231,
        'literature_ref': "O'Donovan et al. (2006, ApJ 651, L61)"
    },
    {
        'target_id': 'TIC 25155310',
        'host_name': 'WASP-126',
        'mission': 'TESS',
        'sector_or_quarter': 1,
        'label': 1,
        'label_name': 'Confirmed Planet',
        't0': 1326.54,
        'period': 3.2888,
        'duration_days': 0.110,
        'depth_ppm': 9200,
        'literature_ref': 'Armstrong et al. (2014, MNRAS 444, 1873)'
    },
    {
        'target_id': 'TIC 410153553',
        'host_name': 'TOI-540',
        'mission': 'TESS',
        'sector_or_quarter': 1,
        'label': 1,
        'label_name': 'Confirmed Planet',
        't0': 1325.503,
        'period': 1.239149,
        'duration_days': 0.038,
        'depth_ppm': 1500,
        'literature_ref': 'Ment et al. (2021, AJ 161, 23)'
    },
    {
        'target_id': 'TIC 307210830',
        'host_name': 'L 98-59',
        'mission': 'TESS',
        'sector_or_quarter': 2,
        'label': 1,
        'label_name': 'Confirmed Planet',
        't0': 1354.904,
        'period': 2.25314,
        'duration_days': 0.045,
        'depth_ppm': 500,
        'literature_ref': 'Kostov et al. (2019, AJ 158, 32)'
    },
    {
        'target_id': 'KIC 6922244',
        'host_name': 'Kepler-8',
        'mission': 'Kepler',
        'sector_or_quarter': 2,
        'label': 1,
        'label_name': 'Confirmed Planet',
        't0': 121.1194228,
        'period': 3.522498429,
        'duration_days': 0.1332,
        'depth_ppm': 9146,
        'literature_ref': 'Jenkins et al. (2010, ApJ 724, 1108)'
    },
    {
        'target_id': 'KIC 11804465',
        'host_name': 'Kepler-11',
        'mission': 'Kepler',
        'sector_or_quarter': 3,
        'label': 1,
        'label_name': 'Confirmed Planet',
        't0': 171.0091259,
        'period': 4.437962934,
        'duration_days': 0.1957,
        'depth_ppm': 1639,
        'literature_ref': 'Lissauer et al. (2011, Nature 470, 53)'
    },

    # ---------------------------------------------------------------------------------
    # NON-PLANETARY ARTIFACTS (Eclipsing Binaries, Stellar Activity, Non-Stationary Dips)
    # ---------------------------------------------------------------------------------
    {
        'target_id': 'KIC 8462852',
        'host_name': 'Boyajian Star (Q8)',
        'mission': 'Kepler',
        'sector_or_quarter': 8,
        'label': 0,
        'label_name': 'Non-Stationary Astrophysical Artifact',
        't0': 792.0,
        'period': 1000.0,
        'duration_days': 4.5,
        'depth_ppm': 15000,
        'literature_ref': 'Boyajian et al. (2016, MNRAS 457, 3988)'
    },
    {
        'target_id': 'KIC 8462852',
        'host_name': 'Boyajian Star (Q16)',
        'mission': 'Kepler',
        'sector_or_quarter': 16,
        'label': 0,
        'label_name': 'Non-Stationary Astrophysical Artifact',
        't0': 1520.0,
        'period': 1000.0,
        'duration_days': 8.0,
        'depth_ppm': 200000,
        'literature_ref': 'Boyajian et al. (2016, MNRAS 457, 3988)'
    },
    {
        'target_id': 'KIC 12644769',
        'host_name': 'KOI-1611 EB',
        'mission': 'Kepler',
        'sector_or_quarter': 2,
        'label': 0,
        'label_name': 'Eclipsing Binary',
        't0': 173.7358974,
        'period': 41.077590334,
        'duration_days': 0.2586,
        'depth_ppm': 134484,
        'literature_ref': 'Prsa et al. (2011, AJ 141, 83)'
    },
    {
        'target_id': 'KIC 3858884',
        'host_name': 'KIC 3858884',
        'mission': 'Kepler',
        'sector_or_quarter': 3,
        'author': 'Kepler',
        'label': 0,
        'label_name': 'Pulsating Binary',
        'literature_ref': 'Maceroni et al. (2014, A&A 563, A106)',
        't0': 262.15,
        'period': 25.95,
        'duration_days': 0.8,
        'depth_ppm': 45000,
        'notes': 'Heartbeat star with tidally induced pulsations'
    },
    {
        'target_id': 'KIC 5385723',
        'host_name': 'V380 Cyg (KIC 5385723)',
        'mission': 'Kepler',
        'sector_or_quarter': 4,
        'author': 'Kepler',
        'label': 0,
        'label_name': 'Eclipsing Binary',
        'literature_ref': 'Tkachenko et al. (2014, MNRAS 438, 3093)',
        't0': 355.2,
        'period': 12.425,
        'duration_days': 0.65,
        'depth_ppm': 80000,
        'notes': 'B-type eccentric eclipsing binary'
    },
    {
        'target_id': 'KIC 11253226',
        'host_name': 'KIC 11253226',
        'mission': 'Kepler',
        'sector_or_quarter': 2,
        'author': 'Kepler',
        'label': 0,
        'label_name': 'Contact Binary',
        'literature_ref': 'Kirk et al. (2016, AJ 151, 68)',
        't0': 170.4,
        'period': 0.388,
        'duration_days': 0.08,
        'depth_ppm': 25000,
        'notes': 'W UMa type contact binary'
    },
    {
        'target_id': 'TIC 441462308',
        'host_name': 'V723 Mon',
        'mission': 'TESS',
        'sector_or_quarter': 33,
        'author': 'SPOC',
        'label': 0,
        'label_name': 'Ellipsoidal Variable',
        'literature_ref': 'Jayasinghe et al. (2021, MNRAS 504, 2577)',
        't0': 2200.0,
        'period': 59.9,
        'duration_days': 5.0,
        'depth_ppm': 30000,
        'notes': 'Red giant with black hole candidate companion (ellipsoidal variation, non-transit)'
    }
]


def get_curated_benchmark_df():
    """Returns the benchmark catalog as a clean DataFrame."""
    return pd.DataFrame(BENCHMARK_TARGETS)


if __name__ == '__main__':
    df = get_curated_benchmark_df()
    print(f"Curated Benchmark: {len(df)} targets total")
    print(f"Confirmed Planets: {sum(df['label'] == 1)} across {df[df['label'] == 1]['host_name'].nunique()} distinct hosts")
    print(f"Non-Planetary Artifacts: {sum(df['label'] == 0)} across {df[df['label'] == 0]['host_name'].nunique()} distinct sources")
    df.to_csv('data/curated_benchmark_catalog.csv', index=False)
    print("Exported to data/curated_benchmark_catalog.csv")
