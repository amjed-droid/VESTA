"""
reproduce_table1_roc.py — Reproduces Table 1 & Figure 3 (ROC) from the manuscript.

Replaces: reproduce_table3_roc.m

Log-SNRAS: A Computationally Efficient Variance-Stabilized Metric
for Vetting Heteroscedastic Light Curves

This script loads the curated multi-host benchmark evaluation results,
computes the AUC for all 8 metrics (including shape-corrected Log-SNRAS),
computes stratified bootstrap 95% confidence intervals (B=2000, seed=42),
performs exact DeLong tests, and generates the publication ROC curve.

Usage
-----
    python reproduce_table1_roc.py
    # or as part of reproduce_all.py (called automatically)
"""

from __future__ import annotations

import sys
import os
import io
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Add repo root to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scripts.utils.delong import delong_test, compute_auc
from sklearn.metrics import roc_curve


def reproduce_table1_roc(
        csv_path: str = 'data/multi_host_evaluation_results.csv',
        fig_out:  str = 'figures/Figure3_ROC_MultiHost.png',
        csv_out:  str = 'data/multi_host_auc_performance.csv',
        B: int = 2000,
        seed: int = 42,
) -> pd.DataFrame:
    """
    Full reproduction of Table 1 and Figure 3.

    Parameters
    ----------
    csv_path : path to evaluation results CSV
    fig_out  : output path for ROC figure (PNG)
    csv_out  : output path for performance summary (CSV)
    B        : bootstrap resamples (default 2000)
    seed     : random seed (default 42, must match manuscript)

    Returns
    -------
    summary DataFrame
    """
    print('=' * 68)
    print('  REPRODUCING BENCHMARK CLASSIFICATION PERFORMANCE (TABLE 1 & ROC)')
    print('=' * 68)

    # ------------------------------------------------------------------
    # 1. Load evaluation data
    # ------------------------------------------------------------------
    if not os.path.isfile(csv_path):
        raise FileNotFoundError(
            f'{csv_path} not found.\n'
            'Run "python scripts/evaluate_multi_host_benchmark.py" first.'
        )

    df = pd.read_csv(csv_path)
    print(f'\nLoaded {len(df)} targets from: {csv_path}\n')

    y = df['label'].values.astype(int)

    # ------------------------------------------------------------------
    # 2. Define metrics
    # ------------------------------------------------------------------
    metrics = [
        ('Penalty (Shape-Corrected)', -df['penalty_corr'].values),
        ('Penalty (Uncorrected)',     -df['penalty_raw'].values),
        ('Traditional SNR',           df['t_snr'].values),
        ('Robust SNR (MAD)',          df['r_snr'].values),
        ('Pont SNR (2006)',           df['p_snr'].values),
        ('BLS SNR Proxy',             df['b_snr'].values),
        ('Log-SNRAS Composite',       df['l_snras_corr'].values),
        ('Inverse Depth (ppm)',      -df['depth_measured_ppm'].values),
    ]

    # ------------------------------------------------------------------
    # 3. Point AUCs (Wilcoxon / U-statistic)
    # ------------------------------------------------------------------
    point_aucs = np.array([compute_auc(y, s) for _, s in metrics])

    # ------------------------------------------------------------------
    # 4. Stratified bootstrap 95% CIs  (B=2000, seed=42)
    # ------------------------------------------------------------------
    rng = np.random.default_rng(seed)
    pos_idx = np.where(y == 1)[0]
    neg_idx = np.where(y == 0)[0]

    boot_aucs = np.zeros((B, len(metrics)))
    for b in range(B):
        b_pos = rng.choice(pos_idx, size=len(pos_idx), replace=True)
        b_neg = rng.choice(neg_idx, size=len(neg_idx), replace=True)
        b_idx = np.concatenate([b_pos, b_neg])
        y_b   = y[b_idx]
        for m, (_, s) in enumerate(metrics):
            boot_aucs[b, m] = compute_auc(y_b, s[b_idx])

    ci_lo = np.percentile(boot_aucs, 2.5,  axis=0)
    ci_hi = np.percentile(boot_aucs, 97.5, axis=0)

    # ------------------------------------------------------------------
    # 5. DeLong test: shape-corrected penalty vs. each baseline
    # ------------------------------------------------------------------
    ref_scores = metrics[0][1]
    delong_z   = np.zeros(len(metrics))
    delong_p   = np.ones(len(metrics))
    for m, (_, s) in enumerate(metrics):
        if m > 0:
            z, p = delong_test(y, ref_scores, s)
            delong_z[m], delong_p[m] = z, p

    # ------------------------------------------------------------------
    # 6. Print Table 1
    # ------------------------------------------------------------------
    print(f'{"Method":<32} | {"AUC":>6} | {"95% CI":^18} | {"Delta AUC":>9} | {"z":>6} | {"p":>7}')
    print('-' * 85)
    for m, (name, _) in enumerate(metrics):
        delta = point_aucs[0] - point_aucs[m] if m > 0 else 0.0
        print(
            f'{name:<32} | {point_aucs[m]:>6.3f} | '
            f'[{ci_lo[m]:.3f}, {ci_hi[m]:.3f}] | '
            f'{delta:>+9.3f} | '
            f'{delong_z[m]:>6.2f} | '
            f'{delong_p[m]:>7.4f}'
        )
    print('=' * 85)

    # ------------------------------------------------------------------
    # 7. Build summary DataFrame  (exact format used in reproduce_all.py)
    # ------------------------------------------------------------------
    delta_aucs = point_aucs[0] - point_aucs
    delta_aucs[0] = 0.0

    summary = pd.DataFrame({
        'Metric':               [n for n, _ in metrics],
        'AUC':                  point_aucs,
        '95% CI Lower':         ci_lo,
        '95% CI Upper':         ci_hi,
        'Delta AUC vs Penalty': delta_aucs,
        'DeLong z':             delong_z,
        'DeLong p-value':       delong_p,
    })
    summary.to_csv(csv_out, index=False)
    print(f'\nSaved statistical summary to: {csv_out}')

    # ------------------------------------------------------------------
    # 8. Generate ROC figure (Figure 3)
    # ------------------------------------------------------------------
    _plot_roc(y, metrics, point_aucs, fig_out)
    return summary


def _plot_roc(y, metrics, point_aucs, outpath):
    """Publication-quality multi-metric ROC curve."""
    style = {
        0: ('#003399', 2.5, '-',  'Penalty (Shape-Corrected)'),
        1: ('#999999', 1.4, '--', 'Penalty (Uncorrected)'),
        2: ('#CC3300', 1.4, '-',  'Traditional SNR'),
        6: ('#009933', 1.8, '-',  'Log-SNRAS Composite'),
        7: ('#FF8800', 1.4, '--', 'Inverse Depth'),
    }

    fig, ax = plt.subplots(figsize=(8, 6.5))
    ax.plot([0, 1], [0, 1], 'k--', lw=1.2, label='Random Chance (AUC = 0.500)')

    for m_idx, (color, lw, ls, _) in style.items():
        name, s = metrics[m_idx]
        fpr, tpr, _ = roc_curve(y, s)
        ax.plot(fpr, tpr, color=color, lw=lw, ls=ls,
                label=f'{name} (AUC = {point_aucs[m_idx]:.3f})')

    ax.set_xlabel('False Positive Rate (1 − Specificity)', fontsize=12, fontweight='bold')
    ax.set_ylabel('True Positive Rate (Sensitivity)',      fontsize=12, fontweight='bold')
    ax.set_title('Multi-Host Benchmark ROC Analysis (N = 17)',
                 fontsize=13, fontweight='bold')
    ax.legend(loc='lower right', fontsize=9)
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1.02])
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved ROC figure to: {outpath}')


if __name__ == '__main__':
    reproduce_table1_roc()
