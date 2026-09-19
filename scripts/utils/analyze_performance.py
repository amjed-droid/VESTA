"""
analyze_performance.py — Full statistical validation pipeline for Log-SNRAS.

Replaces scripts/utils/analyze_performance.m

Implements the complete benchmarking procedure underlying Table 1 and
Section 3 of the manuscript "Log-SNRAS: A Computationally Efficient
Variance-Stabilized Metric for Vetting Heteroscedastic Light Curves."

Methodology
-----------
1. Point-estimate AUC for all 8 metrics.
2. Stratified bootstrap (B=2000, seed=42) 95% confidence intervals.
3. Jackknife stability analysis on the shape-corrected penalty.
4. Exact DeLong test: shape-corrected penalty vs. every baseline.
5. Stratified permutation test (N=10 000 shuffles).

Usage
-----
    from scripts.utils.analyze_performance import run_full_analysis
    results = run_full_analysis(df)

    # or from CLI:
    python scripts/utils/analyze_performance.py
"""

from __future__ import annotations

import sys
import io
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# -- local
sys.path.insert(0, '.')
from scripts.utils.delong import delong_test, compute_auc


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run_full_analysis(df: pd.DataFrame,
                      B: int = 2000,
                      n_perm: int = 10_000,
                      seed: int = 42,
                      save_roc_fig: str | None = 'figures/Figure3_ROC_MultiHost.png',
                      verbose: bool = True) -> pd.DataFrame:
    """
    Run the full statistical validation pipeline.

    Parameters
    ----------
    df   : DataFrame with columns:
               label, penalty_corr, penalty_raw, t_snr, r_snr, p_snr,
               b_snr, l_snras_corr, depth_measured_ppm
    B    : number of bootstrap resamples (default 2000)
    n_perm : number of permutations for permutation test (default 10 000)
    seed : random seed for exact reproducibility (default 42)
    save_roc_fig : path to save the ROC figure (None = don't save)
    verbose : print progress and results

    Returns
    -------
    summary : DataFrame with columns
              [Metric, AUC, 95% CI Lower, 95% CI Upper,
               Delta AUC vs Penalty, DeLong z, DeLong p-value]
    """
    rng = np.random.default_rng(seed)

    y = df['label'].values.astype(int)

    # ------------------------------------------------------------------
    # Metric definitions  (name, scores, direction)
    # direction = +1 → higher score = more likely planet
    # direction = -1 → lower score  = more likely planet (penalty-type)
    # ------------------------------------------------------------------
    metrics = [
        ('Penalty (Shape-Corrected)', -df['penalty_corr'].values),
        ('Penalty (Uncorrected)',     -df['penalty_raw'].values),
        ('Traditional SNR',           df['t_snr'].values),
        ('Robust SNR (MAD)',          df['r_snr'].values),
        ('Pont SNR (2006)',           df['p_snr'].values),
        ('BLS SNR Proxy',             df['b_snr'].values),
        ('Log-SNRAS Composite',       df['l_snras_corr'].values),
        ('Inverse Depth (ppm)',       -df['depth_measured_ppm'].values),
    ]

    # Reference = shape-corrected penalty (index 0)
    ref_scores = metrics[0][1]
    ref_auc    = compute_auc(y, ref_scores)

    # ------------------------------------------------------------------
    # 1. Point AUCs
    # ------------------------------------------------------------------
    point_aucs = np.array([compute_auc(y, s) for _, s in metrics])

    # ------------------------------------------------------------------
    # 2. Stratified bootstrap 95% CI  (B = 2000)
    # ------------------------------------------------------------------
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
    # 3. Jackknife stability on shape-corrected penalty
    # ------------------------------------------------------------------
    n = len(y)
    jack_aucs = np.array([
        compute_auc(np.delete(y, i), np.delete(ref_scores, i))
        for i in range(n)
    ])

    # ------------------------------------------------------------------
    # 4. DeLong test: penalty vs. each other metric
    # ------------------------------------------------------------------
    delong_z = np.zeros(len(metrics))
    delong_p = np.ones(len(metrics))
    for m, (_, s) in enumerate(metrics):
        if m == 0:
            delong_z[m], delong_p[m] = 0.0, 1.0
        else:
            z, p = delong_test(y, ref_scores, s)
            delong_z[m], delong_p[m] = z, p

    # ------------------------------------------------------------------
    # 5. Permutation test (10 000 shuffles)
    # ------------------------------------------------------------------
    perm_aucs = np.array([
        compute_auc(rng.permutation(y), ref_scores)
        for _ in range(n_perm)
    ])
    perm_p = float(np.mean(perm_aucs >= ref_auc))

    # ------------------------------------------------------------------
    # Print results
    # ------------------------------------------------------------------
    if verbose:
        _print_results(metrics, point_aucs, ci_lo, ci_hi,
                       delong_z, delong_p,
                       jack_aucs, perm_p, ref_auc)

    # ------------------------------------------------------------------
    # 6. Generate ROC figure
    # ------------------------------------------------------------------
    if save_roc_fig:
        _plot_roc(y, metrics, point_aucs, save_roc_fig)

    # ------------------------------------------------------------------
    # Build summary DataFrame
    # ------------------------------------------------------------------
    delta_aucs = ref_auc - point_aucs
    delta_aucs[0] = 0.0

    summary = pd.DataFrame({
        'Metric':                [name for name, _ in metrics],
        'AUC':                   point_aucs,
        '95% CI Lower':          ci_lo,
        '95% CI Upper':          ci_hi,
        'Delta AUC vs Penalty':  delta_aucs,
        'DeLong z':              delong_z,
        'DeLong p-value':        delong_p,
    })
    return summary


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _print_results(metrics, point_aucs, ci_lo, ci_hi,
                   delong_z, delong_p,
                   jack_aucs, perm_p, ref_auc):
    """Print a formatted results table to stdout."""
    out = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sep = '=' * 68

    print(sep)
    print('  BENCHMARK STATISTICAL PERFORMANCE SUMMARY')
    print(sep)
    fmt = '{:<30} | AUC = {:.3f} [{:.3f}, {:.3f}] | z = {:5.2f}, p = {:.4f}'
    for m, (name, _) in enumerate(metrics):
        print(fmt.format(name, point_aucs[m], ci_lo[m], ci_hi[m],
                         delong_z[m], delong_p[m]))

    print()
    print('JACKKNIFE STABILITY (Shape-Corrected Penalty)')
    print(f'  Mean = {jack_aucs.mean():.4f} +/- {jack_aucs.std():.6f}')

    print()
    print('PERMUTATION TEST (10 000 shuffles)')
    print(f'  Empirical p-value = {perm_p:.4f}')
    print(sep)


def _plot_roc(y, metrics, point_aucs, outpath):
    """Generate and save the multi-metric ROC curve (Figure 3)."""
    highlight = [0, 1, 2, 6, 7]   # indices to draw explicitly
    colors = {
        0: ('#003399', 2.5),   # Shape-Corrected Penalty — dark blue, thick
        1: ('#888888', 1.4),   # Uncorrected Penalty
        2: ('#CC3300', 1.4),   # Traditional SNR
        6: ('#009933', 1.8),   # Log-SNRAS Composite — green
        7: ('#FF8800', 1.4),   # Inverse Depth
    }

    fig, ax = plt.subplots(figsize=(8, 6.5))
    ax.plot([0, 1], [0, 1], 'k--', lw=1.2, label='Random Chance (AUC = 0.500)')

    for m in highlight:
        name, s = metrics[m]
        fpr, tpr, _ = roc_curve(y, s)
        color, lw = colors[m]
        label = f'{name} (AUC = {point_aucs[m]:.3f})'
        ax.plot(fpr, tpr, color=color, lw=lw, label=label)

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


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    import os
    csv_path = 'data/multi_host_evaluation_results.csv'
    if not os.path.isfile(csv_path):
        print(f'ERROR: {csv_path} not found.')
        print('Run "python scripts/evaluate_multi_host_benchmark.py" first.')
        sys.exit(1)

    df = pd.read_csv(csv_path)
    summary = run_full_analysis(df, B=2000, n_perm=10_000, seed=42)

    out_csv = 'data/multi_host_auc_performance.csv'
    summary.to_csv(out_csv, index=False)
    print(f'\nSaved statistical summary to: {out_csv}')
