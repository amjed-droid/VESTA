"""
Verification script: checks ALL numerical claims in main.tex against reproduce_all.py output.
"""
import pandas as pd
import numpy as np
import os
from sklearn.metrics import roc_auc_score

df_eval = pd.read_csv('data/multi_host_evaluation_results.csv')
perf = pd.read_csv('data/multi_host_auc_performance.csv')

print('=== VERIFICATION REPORT: reproduce_all.py vs main.tex ===\n')

all_ok = True

# -----------------------------------------------------------------------
# TABLE 1 (tab:performance_summary) -- AUC values from manuscript
# -----------------------------------------------------------------------
# Exact mapping: CSV Metric name -> (expected_auc, expected_ci_lo, expected_ci_hi)
expected_aucs = {
    'Penalty (Shape-Corrected)': (0.757, 0.486, 0.971),
    'Penalty (Uncorrected)':     (0.600, 0.286, 0.871),
    'Log-SNRAS Composite':       (0.700, 0.414, 0.943),
    'Traditional SNR':           (0.629, 0.321, 0.914),
    'Robust SNR (MAD)':          (0.614, 0.314, 0.886),
    'BLS SNR Proxy':             (0.614, 0.314, 0.886),
    'Pont SNR (2006)':           (0.600, 0.271, 0.886),
    'Inverse Depth (ppm)':       (0.471, 0.143, 0.857),
}

print(f'TABLE 1 (Performance Summary) AUC Comparison:')
print(f'  {"Metric":<32} | {"Paper":>6} | {"Script":>6} | {"Paper CI":<15} | {"Script CI":<15} | Status')
print('  ' + '-'*95)

for _, row in perf.iterrows():
    name = str(row['Metric'])
    auc = round(float(row['AUC']), 3)
    ci_lo = round(float(row['95% CI Lower']), 3)
    ci_hi = round(float(row['95% CI Upper']), 3)

    # Exact match by name
    matched = False
    for exp_name, (exp_auc, exp_lo, exp_hi) in expected_aucs.items():
        if exp_name.lower() in name.lower() or name.lower() in exp_name.lower():
            auc_ok  = abs(auc  - exp_auc) < 0.002
            lo_ok   = abs(ci_lo - exp_lo) < 0.002
            hi_ok   = abs(ci_hi - exp_hi) < 0.002
            ok = auc_ok and lo_ok and hi_ok
            status = '[OK]' if ok else '[MISMATCH!]'
            if not ok:
                all_ok = False
            paper_ci = f'[{exp_lo:.3f},{exp_hi:.3f}]'
            script_ci = f'[{ci_lo:.3f},{ci_hi:.3f}]'
            print(f'  {name:<32} | {exp_auc:>6.3f} | {auc:>6.3f} | {paper_ci:<15} | {script_ci:<15} | {status}')
            matched = True
            break
    if not matched:
        print(f'  {name:<32} | {"?":>6} | {auc:>6.3f} | {"?":15} | [{"[%.3f,%.3f]" % (ci_lo,ci_hi)}:<15] | [unmatched]')

print()

# -----------------------------------------------------------------------
# LOHO Cross-Validation
# -----------------------------------------------------------------------
print('LOHO Cross-Validation (Section 3.3):')
y = df_eval['label'].values
p_corr = df_eval['penalty_corr'].values
host_names = df_eval['host_name'].values

loho_claims = [
    ('Pi Mensae', 0.762),
    ('TOI-201', 0.778),
]
loho_aucs = []
for host in df_eval['host_name'].unique():
    mask = host_names != host
    if len(np.unique(y[mask])) == 2:
        auc_fold = roc_auc_score(y[mask], -p_corr[mask])
        loho_aucs.append(auc_fold)

for host_claim, expected_auc in loho_claims:
    mask = host_names != host_claim
    if len(np.unique(y[mask])) < 2:
        print(f'  Excluding {host_claim}: [SKIP - insufficient classes]')
        continue
    produced = round(roc_auc_score(y[mask], -p_corr[mask]), 3)
    ok = abs(produced - expected_auc) < 0.002
    status = '[OK]' if ok else '[MISMATCH!]'
    if not ok:
        all_ok = False
    print(f'  Excluding {host_claim:<15}: Paper={expected_auc:.3f} | Script={produced:.3f} {status}')

mean_loho = round(np.mean(loho_aucs), 3)
std_loho  = round(np.std(loho_aucs), 3)
mean_ok   = abs(mean_loho - 0.757) < 0.002
std_ok    = abs(std_loho  - 0.036) < 0.002
status = '[OK]' if (mean_ok and std_ok) else '[MISMATCH!]'
if not (mean_ok and std_ok):
    all_ok = False
print(f'  Mean LOHO AUC         : Paper=0.757 +/-0.036 | Script={mean_loho:.3f} +/-{std_loho:.3f} {status}')
print()

# -----------------------------------------------------------------------
# Individual Target Penalty Values from manuscript body (Section 4)
# -----------------------------------------------------------------------
print('Individual Target Penalties (Section 4, manuscript body):')
claims_by_host = {
    'TOI-201':   {'p_raw': 1.640, 'p_corr': 0.498, 'tol': 0.05},
    'Kepler-8':  {'p_raw': 0.951, 'p_corr': 0.004, 'tol': 0.05},
    'Kepler-11': {'p_raw': 1.552, 'p_corr': 0.031, 'tol': 0.10},
    'Kepler-448':{'p_raw': 1.707, 'p_corr': 0.189, 'tol': 0.05},
}

for _, row in df_eval.iterrows():
    host = str(row['host_name'])
    p_raw  = float(row['penalty_raw'])
    p_corr_val = float(row['penalty_corr'])

    for key, claim in claims_by_host.items():
        if key in host:
            tol = claim['tol']
            raw_ok  = 'N/A' if claim['p_raw'] is None else abs(p_raw - claim['p_raw']) < tol
            corr_ok = abs(p_corr_val - claim['p_corr']) < tol
            ok = (raw_ok is True or raw_ok == 'N/A') and corr_ok
            status = '[OK]' if ok else '[MISMATCH!]'
            if not ok:
                all_ok = False
            print(f'  {key:<12}: P_raw={p_raw:.3f} (paper={claim["p_raw"]}) | '
                  f'P_corr={p_corr_val:.4f} (paper={claim["p_corr"]}) {status}')
            break

    if 'Boyajian' in host:
        range_ok = 1.87 <= p_corr_val <= 2.00
        status = '[OK]' if range_ok else '[MISMATCH!]'
        if not range_ok:
            all_ok = False
        print(f'  Boyajian     : P_corr={p_corr_val:.4f} (paper claims 1.87-1.99) {status}')

print()

# -----------------------------------------------------------------------
# DeLong Test Key Claim (tab:delong_penalty)
# -----------------------------------------------------------------------
print('DeLong Test Table (tab:delong_penalty):')
# From perf CSV: row for Inverse Depth should have z=1.81, p=0.070
for _, row in perf.iterrows():
    if 'Inverse' in str(row['Metric']):
        z  = round(float(row['DeLong z']), 2)
        p  = round(float(row['DeLong p-value']), 3)
        da = round(float(row['Delta AUC vs Penalty']), 3)
        z_ok  = abs(z  - 1.81) < 0.01
        p_ok  = abs(p  - 0.070) < 0.002
        da_ok = abs(da - 0.286) < 0.002
        ok = z_ok and p_ok and da_ok
        status = '[OK]' if ok else '[MISMATCH!]'
        if not ok:
            all_ok = False
        print(f'  Penalty vs Inverse Depth: DeltaAUC={da:.3f} (paper=0.286), z={z:.2f} (1.81), p={p:.3f} (0.070) {status}')

print()

# -----------------------------------------------------------------------
# Figure files
# -----------------------------------------------------------------------
print('FIGURE FILES (Figures 1-6 in manuscript):')
expected_figs = [
    ('Figure1_ResponseSurface.png',      200),
    ('Figure2_EmpiricalDiagnostic.png',  200),
    ('Figure3_ROC_MultiHost.png',        200),
    ('Figure4_PiMen_Real_PhaseFold.png', 300),
    ('Figure5_TOI201_ShapeDiagnostic.png', 300),
    ('Figure6_Comparison_RealData.png',  300),
]
for figname, min_kb in expected_figs:
    path = os.path.join('figures', figname)
    if os.path.exists(path):
        size_kb = os.path.getsize(path) / 1024
        ok = size_kb >= min_kb
        status = '[OK]' if ok else '[TOO SMALL?]'
        if not ok:
            all_ok = False
        print(f'  {figname}: {size_kb:.0f} KB {status}')
    else:
        print(f'  {figname}: MISSING [FAIL]')
        all_ok = False

print()
print('=' * 70)
if all_ok:
    print('  OVERALL: ALL CHECKS PASSED - Reproduction matches the manuscript!')
else:
    print('  OVERALL: SOME MISMATCHES DETECTED - See items above.')
print('=' * 70)
