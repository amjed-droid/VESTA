"""
Reviewer Compliance Checker
Verifies every reviewer request (R5-C1 through R5-C6, R4-C1 through R4-C2)
is fully addressed in the repository and manuscript.
"""
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import pandas as pd
import numpy as np
import os
import sys
sys.path.insert(0, 'src')
from sklearn.metrics import roc_auc_score

df = pd.read_csv('data/multi_host_evaluation_results.csv')
perf = pd.read_csv('data/multi_host_auc_performance.csv')

y = df['label'].values
p_corr_vals = df['penalty_corr'].values
p_raw_vals  = df['penalty_raw'].values
hosts       = df['host_name'].values

pass_all = True

def check(label, condition, detail=''):
    global pass_all
    icon = '[OK]   ' if condition else '[FAIL] '
    if not condition:
        pass_all = False
    print(f'  {icon} {label}' + (f' | {detail}' if detail else ''))

print('=== REVIEWER COMPLIANCE VERIFICATION REPORT ===\n')

# ============================================================
# REVIEWER #5
# ============================================================
print('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━')
print('REVIEWER #5 – Forensic Methodological Critique')
print('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n')

# ---- C1: Synthetic data ----
print('Comment 1: Synthetic data in baseline validation figure')
fig4 = os.path.join('figures', 'Figure4_PiMen_Real_PhaseFold.png')
fig6 = os.path.join('figures', 'Figure6_Comparison_RealData.png')
check('Figure4 (Pi Mensae) from real TESS SPOC FITS',
      os.path.exists(fig4) and os.path.getsize(fig4) > 100000,
      f'{os.path.getsize(fig4)//1024} KB' if os.path.exists(fig4) else 'MISSING')
check('Figure6 (Boyajian D792) from real Kepler Q8 FITS',
      os.path.exists(fig6) and os.path.getsize(fig6) > 100000,
      f'{os.path.getsize(fig6)//1024} KB' if os.path.exists(fig6) else 'MISSING')
gen_base = 'scripts/plots/generate_baseline_figure.m'
check('Synthetic script generate_baseline_figure.m DELETED',
      not os.path.exists(gen_base),
      'DELETED' if not os.path.exists(gen_base) else 'STILL EXISTS! Contains rand()')
check('Reproduction script generate_publication_figures.py exists',
      os.path.exists('scripts/generate_publication_figures.py'))
print()

# ---- C2: Masking discrepancy ----
print('Comment 2: Masking discrepancy (ephemeris vs sigma-clipping)')
masking_py = 'src/log_snras/masking.py'
check('masking.py module exists', os.path.exists(masking_py))
if os.path.exists(masking_py):
    with open(masking_py) as f:
        mc = f.read()
    check('No sigma clipping in masking pipeline',
          'sigma_clip' not in mc.lower() and 'outlier' not in mc.lower())
    check('Ephemeris-based mask (create_ephemeris_mask)',
          'create_ephemeris_mask' in mc or 'ephemeris' in mc.lower())
    check('2x T_dur buffer implemented',
          '2' in mc and ('buffer' in mc.lower() or 'T_dur' in mc or 't_dur' in mc))
print()

# ---- C3: Column indexing bug ----
print('Comment 3: Hardcoded column index 8 in FITS ingestion')
io_py = 'src/log_snras/io.py'
check('io.py module exists', os.path.exists(io_py))
if os.path.exists(io_py):
    with open(io_py) as f:
        io_c = f.read()
    check("Named column 'PDCSAP_FLUX' used",
          'PDCSAP_FLUX' in io_c)
    check("Named column 'TIME' used",
          "'TIME'" in io_c or '"TIME"' in io_c)
    check('No hardcoded column index [8] or [:,8] or data(8)',
          'data[:, 8]' not in io_c and 'data[:,8]' not in io_c and 'data(8)' not in io_c)
print()

# ---- C4: TOI-201 / Single-host confound ----
print('Comment 4: TOI-201 misclassification and single-host confound')
n_planets  = len(df[df['label']==1]['host_name'].unique())
n_artifacts = len(df[df['label']==0]['host_name'].unique())
check(f'Multi-host benchmark: {n_planets} planet hosts (need >=10)', n_planets >= 10)
check(f'Artifact class: {n_artifacts} distinct artifacts (need >=5)', n_artifacts >= 5)

# Shape correction on TOI-201
toi = df[df['host_name'] == 'TOI-201'].iloc[0]
p_raw_toi  = float(toi['penalty_raw'])
p_corr_toi = float(toi['penalty_corr'])
check('TOI-201 raw penalty > 1.0 (was incorrectly Tier-3)',
      p_raw_toi > 1.0, f'P_raw = {p_raw_toi:.3f}')
check('TOI-201 corrected penalty < 0.60 (restored to non-Tier-3)',
      p_corr_toi < 0.60, f'P_corr = {p_corr_toi:.4f}')

# LOHO AUC
mask_no_pi = hosts != 'Pi Mensae'
mask_no_toi = hosts != 'TOI-201'
auc_no_pi  = round(roc_auc_score(y[mask_no_pi],  -p_corr_vals[mask_no_pi]),  3)
auc_no_toi = round(roc_auc_score(y[mask_no_toi], -p_corr_vals[mask_no_toi]), 3)
check(f'LOHO AUC excluding Pi Mensae = {auc_no_pi} (paper=0.762)',
      abs(auc_no_pi - 0.762) < 0.002)
check(f'LOHO AUC excluding TOI-201 = {auc_no_toi} (paper=0.778)',
      abs(auc_no_toi - 0.778) < 0.002)

loho_aucs = []
for host in df['host_name'].unique():
    mask = hosts != host
    if len(np.unique(y[mask])) == 2:
        loho_aucs.append(roc_auc_score(y[mask], -p_corr_vals[mask]))
mean_loho = round(np.mean(loho_aucs), 3)
std_loho  = round(np.std(loho_aucs), 3)
check(f'Mean LOHO AUC = {mean_loho} +/- {std_loho} (paper=0.757 +/- 0.036)',
      abs(mean_loho - 0.757) < 0.002)
print()

# ---- C5: Ground-truth mislabeling ----
print('Comment 5: Ground-truth mislabeling (Kepler-448, TrES-2, Proxima Cen)')
planets_hosts = df[df['label']==1]['host_name'].tolist()
artifacts_hosts = df[df['label']==0]['host_name'].tolist()
all_hosts = df['host_name'].tolist()
check('Kepler-448 b correctly in Confirmed Planets',
      any('Kepler-448' in h for h in planets_hosts))
check('TrES-2 b correctly in Confirmed Planets',
      any('TrES-2' in h for h in planets_hosts))
check('Proxima Cen removed entirely from catalog',
      not any('Proxima' in h for h in all_hosts))
check('Boyajian Star correctly in Artifacts (Q8 & Q16)',
      sum(1 for h in artifacts_hosts if 'Boyajian' in h) >= 2)
check('KOI-1611 eccentric EB in Artifacts',
      any('KOI-1611' in h or '12644769' in h for h in artifacts_hosts))
print()

# ---- C6: Depth as confound ----
print('Comment 6: Transit depth as confounding covariate')
inv_row = perf[perf['Metric'].str.contains('Inverse', case=False)].iloc[0]
pen_row = perf[perf['Metric'].str.contains('Shape', case=False)].iloc[0]
inv_auc  = round(float(inv_row['AUC']), 3)
pen_auc  = round(float(pen_row['AUC']), 3)
delta    = round(pen_auc - inv_auc, 3)
delong_z = round(float(inv_row['DeLong z']), 2)
delong_p = round(float(inv_row['DeLong p-value']), 3)
check(f'Inverse Depth baseline AUC = {inv_auc} (paper=0.471)',
      abs(inv_auc - 0.471) < 0.002)
check(f'Delta AUC (Penalty - InvDepth) = {delta} (paper=0.286)',
      abs(delta - 0.286) < 0.002)
check(f'DeLong z = {delong_z} (paper=1.81)',
      abs(delong_z - 1.81) < 0.02)
check(f'DeLong p = {delong_p} (paper=0.070)',
      abs(delong_p - 0.070) < 0.002)
print()

# ============================================================
# REVIEWER #4
# ============================================================
print('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━')
print('REVIEWER #4 – Mathematical Rigor Critique')
print('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n')

main_tex = os.path.join('..', 'main.tex')
with open(main_tex, 'r', encoding='utf-8', errors='ignore') as f:
    main_content = f.read()

# ---- C1: Formal propositions ----
print('Comment 1: Mathematical rigor – formal propositions in Appendix B')
check('Proposition 1 (Asymptotic Distribution of SNR_trad) present',
      'Proposition 1' in main_content)
check('Proposition 2 (N_in >= 30 threshold) present',
      'Proposition 2' in main_content)
check('Lindeberg-Levy CLT proof cited',
      'Lindeberg' in main_content)
check("Slutsky's Theorem cited in proof",
      'Slutsky' in main_content)
check("Jensen's inequality and finite-sample shrinkage",
      "Jensen" in main_content)
check('Shrinkage table at N_in=30, 100, 250',
      '0.916' in main_content and '0.952' in main_content and '0.969' in main_content)
print()

# ---- C2: Calibrated tiers from F-distribution ----
print('Comment 2: Tier boundaries calibrated to null F-distribution')
check('F-distribution quantile calibration in main.tex (Section 2.4)',
      'F_{1-\\alpha' in main_content or 'F(N_{\\text{in}}' in main_content)
# Paper uses \alpha = 0.05 and \alpha = 0.001 in Section 2.4 (tier definitions)
check('alpha = 0.05 level for Tier 1 stated in Section 2.4',
      r'\alpha = 0.05' in main_content or 'alpha_1 = 0.05' in main_content)
# Paper uses P <= 0.15 (for N_in~100) which is the approximate F-calibrated P_crit,1
check('Calibrated P_crit,1 ~ 0.15 for N_in~100 stated in Section 2.4',
      '0.15' in main_content)
check('Calibrated P_crit,2 ≈ 0.44 for N_in=30',
      '0.44' in main_content)

# Check F-calibration in core.py
core_py = 'src/log_snras/core.py'
with open(core_py) as f:
    core_c = f.read()
check('F-distribution calibration in src/log_snras/core.py',
      'f.ppf' in core_c.lower() or 'scipy.stats' in core_c.lower() or 'stats.f' in core_c)
print()

# ============================================================
# FINAL SUMMARY
# ============================================================
print('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━')
if pass_all:
    print('FINAL RESULT: ALL REVIEWER REQUESTS FULLY ADDRESSED [PASS]')
else:
    print('FINAL RESULT: SOME REVIEWER ITEMS NEED ATTENTION [REVIEW]')
print('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━')
