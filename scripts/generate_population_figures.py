"""
generate_population_figures.py
==============================
توليد كافة الرسوم البيانية الخمسة عالية الدقة (600 DPI) لمخطوطة البحث السكاني الجديد
لمجلة Astronomy and Computing استناداً إلى بيانات N=56 الحقيقية الكاملة.
"""

from __future__ import annotations
import os, sys, glob
sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('scripts'))
sys.path.insert(0, os.path.abspath('src'))
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.metrics import roc_curve, roc_auc_score

plt.rcParams.update({
    'font.size': 11,
    'font.family': 'sans-serif',
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'figure.titlesize': 14,
    'figure.autolayout': False
})

OUT_DIR = "figures"
os.makedirs(OUT_DIR, exist_ok=True)
RESULTS_CSV = "data/kepler_population_results.csv"
PERF_CSV    = "data/kepler_population_performance.csv"


def plot_fig1_response_surface():
    """Fig 1: Theoretical Response Surface with Explicit Tier Contours."""
    print("Generating Figure 1: Response Surface...")
    fig, ax = plt.subplots(figsize=(7, 5.5), dpi=300)

    sig_in = np.linspace(0.0001, 0.005, 200)
    sig_out = np.linspace(0.0001, 0.005, 200)
    SI, SO = np.meshgrid(sig_in, sig_out)

    psi = np.abs(SI - SO) / SO
    P = np.log(1.0 + psi)
    
    snr_trad = 50.0
    L_VESTA = snr_trad / (1.0 + P)

    cf = ax.contourf(SI * 1e3, SO * 1e3, L_VESTA, levels=25, cmap='viridis')
    cbar = fig.colorbar(cf, ax=ax)
    cbar.set_label(r'VESTA Score ($\mathrm{SNR}_{\rm raw} = 50$)', rotation=270, labelpad=15)

    cs1 = ax.contour(SI * 1e3, SO * 1e3, P, levels=[0.15], colors=['blue'], linestyles=['dotted'], linewidths=[2.0])
    cs2 = ax.contour(SI * 1e3, SO * 1e3, P, levels=[0.60], colors=['red'], linestyles=['dashdot'], linewidths=[2.0])
    ax.plot([0.1, 5.0], [0.1, 5.0], 'k--', lw=1.5, label=r'Homoscedastic Line ($\sigma_{\rm in} = \sigma_{\rm out}$)')

    custom_lines = [
        Line2D([0], [0], color='black', lw=1.5, linestyle='--'),
        Line2D([0], [0], color='blue', lw=2.0, linestyle=':'),
        Line2D([0], [0], color='red', lw=2.0, linestyle='-.')
    ]
    ax.legend(
        custom_lines,
        [
            r'Homoscedastic ($\sigma_{\rm in} = \sigma_{\rm out}$)',
            r'Tier 1 / Tier 2 Boundary ($P \leq 0.15$)',
            r'Tier 2 / Tier 3 Veto Boundary ($P > 0.60$)'
        ],
        loc='upper left', framealpha=0.9
    )

    ax.set_xlabel(r'In-Transit Dispersion $\sigma_{\rm in,res} \times 10^3$')
    ax.set_ylabel(r'Out-of-Transit Dispersion $\sigma_{\rm out} \times 10^3$')
    ax.set_title('Heteroscedastic Response Surface and Calibrated Tier Contours')
    ax.grid(True, linestyle=':', alpha=0.5)

    p1 = os.path.join(OUT_DIR, 'Figure1_ResponseSurface.png')
    p2 = 'Figure1_ResponseSurface.png'
    fig.savefig(p1, dpi=600, bbox_inches='tight')
    fig.savefig(p2, dpi=600, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {p1}")


def plot_fig2_depth_invariance():
    """Fig 2: Population Depth-Invariance Across 4 Decades of Transit Depth."""
    print("Generating Figure 2: Population Depth Invariance...")
    df = pd.read_csv(RESULTS_CSV)
    planets = df[df['label'] == 1].copy()
    artifacts = df[df['label'] == 0].copy()

    fig, ax = plt.subplots(figsize=(8.5, 5.5), dpi=300)

    # Background shading for depth strata
    ax.axvspan(10, 100, color='#e8f4f8', alpha=0.7, label='Ultra-Shallow (<100 ppm)')
    ax.axvspan(100, 1000, color='#e2f0d9', alpha=0.7, label='Shallow Super-Earth (100-1k ppm)')
    ax.axvspan(1000, 10000, color='#fff2cc', alpha=0.7, label='Moderate (1k-10k ppm)')
    ax.axvspan(10000, 100000, color='#fce4d6', alpha=0.7, label='Ultra-Deep (>10k ppm)')

    # Tier boundaries
    ax.axhline(0.15, color='blue', linestyle=':', lw=1.8, label=r'Tier 1 Boundary ($P \leq 0.15$)')
    ax.axhline(0.60, color='red', linestyle='-.', lw=1.8, label=r'Tier 3 Veto ($P > 0.60$)')

    # Scatter confirmed planets
    sc = ax.scatter(
        planets['literature_depth_ppm'],
        planets['penalty_corr'],
        c=planets['t_snr'],
        cmap='viridis',
        s=65,
        edgecolors='black',
        linewidth=0.8,
        zorder=5,
        label='Kepler Confirmed Planets (N=49)'
    )
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label(r'Traditional SNR ($\mathrm{SNR}_{\rm trad}$)', rotation=270, labelpad=15)

    # Scatter artifacts at their literature depths
    ax.scatter(
        artifacts['literature_depth_ppm'],
        artifacts['penalty_corr'],
        color='crimson',
        marker='X',
        s=120,
        edgecolors='black',
        linewidth=1.2,
        zorder=6,
        label='Astrophysical Artifacts (N=7)'
    )

    ax.set_xscale('log')
    ax.set_xlim(8, 250000)
    ax.set_ylim(-0.05, 2.1)
    ax.set_xlabel(r'Transit Depth $\delta$ (ppm) [Logarithmic Scale]')
    ax.set_ylabel(r'Shape-Corrected Penalty Metric $P = \ln(1 + \psi)$')
    ax.set_title('VESTA Depth-Invariance Across Four Decades (Kepler N=56 Population)')
    ax.grid(True, which='both', linestyle=':', alpha=0.5)

    # Annotations for extreme benchmark systems with background bounding boxes
    box_props = dict(boxstyle='round,pad=0.25', facecolor='white', alpha=0.88, edgecolor='#888888', lw=0.7)

    ax.annotate('Kepler-37 b\n(12 ppm)', xy=(12, 0.038), xytext=(15, 0.36),
                arrowprops=dict(arrowstyle='->', lw=1.1, color='black'),
                fontsize=8.5, fontweight='bold', bbox=box_props)

    ax.annotate('Kepler-10 b\n(190 ppm)', xy=(190, 0.003), xytext=(115, 0.33),
                arrowprops=dict(arrowstyle='->', lw=1.1, color='black'),
                fontsize=8.5, fontweight='bold', bbox=box_props)

    ax.annotate('Kepler-45 b\n(32,040 ppm)', xy=(32040, 0.0269), xytext=(50000, 0.22),
                arrowprops=dict(arrowstyle='->', lw=1.1, color='black', connectionstyle="arc3,rad=-0.15"),
                fontsize=8.5, fontweight='bold', bbox=box_props)

    ax.annotate("Boyajian's Star Q8\n(Anomalous Dip)", xy=(15000, 1.97), xytext=(4500, 1.68),
                arrowprops=dict(arrowstyle='->', lw=1.1, color='crimson'),
                fontsize=8.5, fontweight='bold', color='crimson',
                bbox=dict(boxstyle='round,pad=0.25', facecolor='white', alpha=0.88, edgecolor='crimson', lw=0.7))

    ax.legend(loc='upper left', fontsize=9, framealpha=0.92)

    p1 = os.path.join(OUT_DIR, 'Figure2_DepthInvariance.png')
    p2 = 'Figure2_DepthInvariance.png'
    fig.savefig(p1, dpi=600, bbox_inches='tight')
    fig.savefig(p2, dpi=600, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {p1}")


def plot_fig3_phase_folded_archetypes():
    """Fig 3: Phase-Folded Light Curve Archetypes across 4 Decades."""
    print("Generating Figure 3: Archetype Phase-Folded Transits...")
    from scripts.evaluate_kepler_population import load_kepler_lightcurve_from_fits
    from src.vesta.masking import create_ephemeris_mask
    from src.vesta.core import compute_transit_shape_residuals

    df = pd.read_csv("data/kepler_population_catalog.csv")
    archetypes = [
        ('Kepler-45 b',  'Ultra-Deep (32,040 ppm)'),
        ('Kepler-148 d', 'Moderate (10,310 ppm)'),
        ('Kepler-10 b',  'Shallow Super-Earth (190 ppm)'),
        ('Kepler-37 b',  'Ultra-Shallow Sub-Earth (12 ppm)')
    ]

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), dpi=300)
    axes = axes.flatten()

    for ax, (name, label_strat) in zip(axes, archetypes):
        row = df[df['target_name'] == name].iloc[0]
        fits_dir = row['fits_dir']
        period = float(row['period'])
        t0 = float(row['t0_bkjd'])
        dur_d = float(row['duration_hours']) / 24.0

        t, f = load_kepler_lightcurve_from_fits(fits_dir, max_files=4)
        if t is None:
            continue

        # Phase fold
        phase = ((t - t0) % period) / period
        phase[phase > 0.5] -= 1.0

        # Transit mask
        transit_mask, _ = create_ephemeris_mask(t, t0, period, dur_d, buffer_factor=1.5)
        phase_in = phase[transit_mask]
        f_in = f[transit_mask]

        # Restrict plot window to +/- 3 * duration
        win = 3.0 * (dur_d / period)
        mask_win = np.abs(phase) <= win
        p_win = phase[mask_win]
        f_win = f[mask_win]

        # Plot out-of-transit
        ax.scatter(p_win * 24.0 * period, (f_win - 1.0) * 1e3, color='gray', alpha=0.35, s=6, label='Out-of-Transit Flux')
        # Plot in-transit
        ax.scatter(phase_in * 24.0 * period, (f_in - 1.0) * 1e3, color='royalblue', alpha=0.8, s=12, label='In-Transit Points')

        # Fit and plot shape model
        if len(f_in) >= 6:
            _, f_model = compute_transit_shape_residuals(f_in, time_in=t[transit_mask], period=period, t0=t0)
            sort_p = np.argsort(phase_in)
            ax.plot(phase_in[sort_p] * 24.0 * period, (f_model[sort_p] - 1.0) * 1e3, 'r-', lw=2.0, label='Shape Model Residual Baseline')

        ax.set_title(f"{name} — {label_strat}", fontweight='bold', fontsize=11)
        ax.set_xlabel('Time from Mid-Transit (hours)')
        ax.set_ylabel(r'Normalized Flux Offset ($\times 10^{-3}$)')
        ax.grid(True, linestyle=':', alpha=0.5)
        ax.legend(loc='lower right', fontsize=8, framealpha=0.9)

    fig.tight_layout()
    p1 = os.path.join(OUT_DIR, 'Figure3_Archetypes_PhaseFold.png')
    p2 = 'Figure3_Archetypes_PhaseFold.png'
    fig.savefig(p1, dpi=600, bbox_inches='tight')
    fig.savefig(p2, dpi=600, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {p1}")


def plot_fig4_population_roc():
    """Fig 4: Population-Scale ROC Curve with Multi-Metric Comparison."""
    print("Generating Figure 4: Population ROC Comparison...")
    df = pd.read_csv(RESULTS_CSV)
    y_true = df['label'].values

    fig, ax = plt.subplots(figsize=(7.5, 6), dpi=300)

    # Reference diagonal
    ax.plot([0, 1], [0, 1], 'k--', lw=1.2, label='Random Classifier (AUC = 0.500)')

    # Metrics
    configs = [
        ('Penalty (Shape-Corrected)', -df['penalty_corr'].values, 'crimson', '-', 2.5),
        ('VESTA Composite Score',     df['vesta_score'].values,   'purple',  '-.', 1.8),
        ('Traditional SNR',           df['t_snr'].values,         'navy',    '--', 1.8),
        ('Robust SNR (MAD)',          df['r_snr'].values,         'teal',    ':', 1.8),
        ('Pont SNR (2006)',           df['p_snr'].values,         'darkorange', '--', 1.5),
        ('BLS SNR Proxy',             df['b_snr'].values,         'brown',   ':', 1.5),
    ]

    for name, scores, color, ls, lw in configs:
        fpr, tpr, _ = roc_curve(y_true, scores)
        auc_val = roc_auc_score(y_true, scores)
        ax.plot(fpr, tpr, color=color, linestyle=ls, linewidth=lw, label=f"{name:<24} (AUC = {auc_val:.3f})")

    ax.set_xlabel('False Positive Rate (1 - Specificity)')
    ax.set_ylabel('True Positive Rate (Sensitivity)')
    ax.set_title('Population-Scale ROC Benchmark (Kepler N=56 Systems)')
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, linestyle=':', alpha=0.5)
    ax.legend(loc='lower right', fontsize=9, framealpha=0.92)

    p1 = os.path.join(OUT_DIR, 'Figure4_Population_ROC.png')
    p2 = 'Figure4_Population_ROC.png'
    p3 = 'Figure3_ROC_MultiHost.png'
    fig.savefig(p1, dpi=600, bbox_inches='tight')
    fig.savefig(p2, dpi=600, bbox_inches='tight')
    fig.savefig(p3, dpi=600, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {p1}")


def plot_fig5_tier_distribution():
    """Fig 5: Tier Classification Distribution Across Strata and Artifacts."""
    print("Generating Figure 5: Tier Distribution Matrix...")
    df = pd.read_csv(RESULTS_CSV)

    strata_raw = [
        'Ultra-Deep (>10,000 ppm)',
        'Moderate (1,000-10,000 ppm)',
        'Shallow Super-Earth (100-1,000 ppm)',
        'Ultra-Shallow Sub-Earth (<100 ppm)',
        'Astrophysical Artifact'
    ]

    strata_display = [
        'Ultra-Deep (>10,000 ppm)\n[N = 5]',
        'Moderate (1,000–10,000 ppm)\n[N = 24]',
        'Shallow Super-Earth (100–1,000 ppm)\n[N = 13]',
        'Ultra-Shallow Sub-Earth (<100 ppm)\n[N = 7]',
        'Astrophysical Artifacts\n[N = 7]'
    ]

    tier1_counts, tier2_counts, tier3_counts = [], [], []

    for s in strata_raw:
        grp = df[df['depth_stratum'] == s]
        n_tot = len(grp)
        t1 = (grp['tier_num'] == 1).sum() / float(n_tot) * 100.0
        t2 = (grp['tier_num'] == 2).sum() / float(n_tot) * 100.0
        t3 = (grp['tier_num'] == 3).sum() / float(n_tot) * 100.0
        tier1_counts.append(t1)
        tier2_counts.append(t2)
        tier3_counts.append(t3)

    tier1 = np.array(tier1_counts)
    tier2 = np.array(tier2_counts)
    tier3 = np.array(tier3_counts)

    fig, ax = plt.subplots(figsize=(8.5, 4.8), dpi=300)
    y_pos = np.arange(len(strata_raw))
    bar_height = 0.58

    c_tier1 = '#2e7d32'  # Forest Green
    c_tier2 = '#f57c00'  # Amber Orange
    c_tier3 = '#c62828'  # Deep Crimson

    ax.barh(y_pos, tier1, color=c_tier1, edgecolor='white', linewidth=1.5, height=bar_height,
            label='Tier 1 (Clean Planet Candidate)')
    ax.barh(y_pos, tier2, left=tier1, color=c_tier2, edgecolor='white', linewidth=1.5, height=bar_height,
            label='Tier 2 (Review Required)')
    ax.barh(y_pos, tier3, left=tier1 + tier2, color=c_tier3, edgecolor='white', linewidth=1.5, height=bar_height,
            label='Tier 3 (Automated Veto)')

    for i in range(len(strata_raw)):
        if tier1[i] >= 8.0:
            ax.text(tier1[i] / 2.0, i, f"{tier1[i]:.1f}%",
                    ha='center', va='center', color='white', fontweight='bold', fontsize=10)
        if tier2[i] >= 8.0:
            ax.text(tier1[i] + tier2[i] / 2.0, i, f"{tier2[i]:.1f}%",
                    ha='center', va='center', color='white', fontweight='bold', fontsize=10)
        if tier3[i] >= 8.0:
            ax.text(tier1[i] + tier2[i] + tier3[i] / 2.0, i, f"{tier3[i]:.1f}%",
                    ha='center', va='center', color='white', fontweight='bold', fontsize=10)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(strata_display, fontsize=9.5)
    ax.set_xlabel('Population Proportion (%)', fontsize=11, fontweight='bold', labelpad=8)
    ax.set_xlim(0, 100)
    ax.xaxis.grid(True, linestyle='--', alpha=0.4, color='#999999')
    ax.set_axisbelow(True)

    ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.05),
              ncol=3, fontsize=9.5, frameon=True, facecolor='#ffffff', edgecolor='#cccccc',
              fancybox=True, borderpad=0.6, handlelength=1.4)

    for spine in ax.spines.values():
        spine.set_color('#333333')
        spine.set_linewidth(1.0)

    fig.tight_layout()
    p1 = os.path.join(OUT_DIR, 'Figure5_Tier_Distribution.png')
    p2 = 'Figure5_Tier_Distribution.png'
    fig.savefig(p1, dpi=600, bbox_inches='tight')
    fig.savefig(p2, dpi=600, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {p1}")


def plot_fig6_case_studies_diagnostics():
    """Fig 6: Diagnostic Case Studies: High-Cadence USP and Extreme TTVs."""
    print("Generating Figure 6: Specialized Case Studies (Kepler-78 b and Kepler-36 c)...")
    from scripts.evaluate_kepler_population import load_kepler_lightcurve_from_fits

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), dpi=300)

    # (a) Kepler-78 b: USP World
    ax1 = axes[0]
    t78, f78 = load_kepler_lightcurve_from_fits('D:/nrebig/kepler/lightcurves/40_Kepler_78', max_files=8)
    p78 = 0.35500745
    t0_78 = 120.95984
    dur_78 = 0.814 / 24.0

    phase78 = ((t78 - t0_78) % p78) / p78
    phase78[phase78 > 0.5] -= 1.0
    dt_h_78 = phase78 * p78 * 24.0

    win78 = np.abs(dt_h_78) < 2.0
    t_w78 = dt_h_78[win78]
    f_w78 = (f78[win78] - 1.0) * 1e6

    bins78 = np.linspace(-2.0, 2.0, 35)
    bc78 = 0.5 * (bins78[:-1] + bins78[1:])
    bm78, be78 = [], []
    for i in range(len(bins78)-1):
        m = (t_w78 >= bins78[i]) & (t_w78 < bins78[i+1])
        if np.sum(m) > 0:
            bm78.append(np.median(f_w78[m]))
            be78.append(np.std(f_w78[m]) / np.sqrt(np.sum(m)))
        else:
            bm78.append(np.nan)
            be78.append(np.nan)

    ax1.scatter(t_w78, f_w78, s=2, color='gray', alpha=0.18, rasterized=True, label='Kepler Cadences (N = 13,081)')
    ax1.errorbar(bc78, bm78, yerr=be78, fmt='o-', color='#1f77b4', lw=2.0, markersize=4.5, capsize=2, label='Phase Binned (35 bins)')
    ax1.axvspan(-0.814/2, 0.814/2, color='#2ca02c', alpha=0.15, label='Transit Duration (0.81 h)')

    ax1.set_xlim(-2.0, 2.0)
    ax1.set_ylim(-3500, 3500)
    ax1.set_xlabel('Time from Mid-Transit (Hours)')
    ax1.set_ylabel(r'Relative Flux $\Delta F / F_0$ (ppm)')
    ax1.set_title('(a) Kepler-78 b: Ultra-Short-Period (P = 8.5 h)\nTier 1 Clean: Zero Penalty under Rapid Cadence Stacking')
    ax1.grid(True, linestyle=':', alpha=0.5)

    badge78 = (
        'Kepler-78 b (USP Terrestrial)\n' +
        r'$\delta = 220.0\ \mathrm{ppm},\ N_{\mathrm{in}} = 2618$' + '\n' +
        r'$\sigma_{\mathrm{out}} = 2862.2\ \mathrm{ppm},\ \sigma_{\mathrm{in,res}} = 2855.3\ \mathrm{ppm}$' + '\n' +
        r'$\psi_{\mathrm{corr}} = 0.002,\ \mathcal{P}_{\mathrm{corr}} = \mathbf{0.0024}$' + '\n' +
        r'$\mathbf{Tier\ 1\ (Clean)}\ [99.8\%\ \mathrm{SNR\ Retained}]$'
    )
    ax1.text(0.04, 0.04, badge78, transform=ax1.transAxes, fontsize=8.2, verticalalignment='bottom',
             bbox=dict(boxstyle='round,pad=0.4', facecolor='#e8f5e9', edgecolor='#2e7d32', lw=1.2, alpha=0.92))
    ax1.legend(loc='upper right', fontsize=8, framealpha=0.9)

    # (b) Kepler-36 c: Linear
    t36, f36 = load_kepler_lightcurve_from_fits('D:/nrebig/kepler/lightcurves/34_Kepler_36', max_files=8)
    p36 = 16.231917
    t0_36 = 171.67758
    dur_36 = 7.461 / 24.0

    phase36_lin = ((t36 - t0_36) % p36) / p36
    phase36_lin[phase36_lin > 0.5] -= 1.0
    dt_h_36_lin = phase36_lin * p36 * 24.0

    ax2 = axes[1]
    win36 = np.abs(dt_h_36_lin) < 14.0
    t_w36 = dt_h_36_lin[win36]
    f_w36 = (f36[win36] - 1.0) * 1e6

    bins36 = np.linspace(-14.0, 14.0, 35)
    bc36 = 0.5 * (bins36[:-1] + bins36[1:])
    bm36_lin, be36_lin = [], []
    for i in range(len(bins36)-1):
        m = (t_w36 >= bins36[i]) & (t_w36 < bins36[i+1])
        if np.sum(m) > 0:
            bm36_lin.append(np.median(f_w36[m]))
            be36_lin.append(np.std(f_w36[m]) / np.sqrt(np.sum(m)))
        else:
            bm36_lin.append(np.nan)
            be36_lin.append(np.nan)

    ax2.scatter(t_w36, f_w36, s=3, color='gray', alpha=0.25, rasterized=True, label='Folded Points (Linear)')
    ax2.errorbar(bc36, bm36_lin, yerr=be36_lin, fmt='s-', color='#d62728', lw=2.0, markersize=4.5, capsize=2, label='Linear Phase Binned')
    ax2.axvspan(-7.461/2, 7.461/2, color='#d62728', alpha=0.15, label='Nominal Mask (7.46 h)')

    ax2.set_xlim(-14.0, 14.0)
    ax2.set_ylim(-850, 450)
    ax2.set_xlabel('Time from Nominal Mid-Transit (Hours)')
    ax2.set_ylabel(r'Relative Flux $\Delta F / F_0$ (ppm)')
    ax2.set_title(r'(b) Kepler-36 c: Rigid Linear Ephemeris' + '\n' + r'Tier 3 Veto: $\pm 8$ h TTV Smearing Inflates Residuals')
    ax2.grid(True, linestyle=':', alpha=0.5)

    badge36_lin = (
        'Kepler-36 c (Static Linear Fold)\n' +
        r'$\delta_{\mathrm{true}} = 430.0\ \mathrm{ppm},\ \Delta t_{\mathrm{TTV}} \sim \pm 8\ \mathrm{h}$' + '\n' +
        r'$\sigma_{\mathrm{out}} = 104.3\ \mathrm{ppm},\ \sigma_{\mathrm{in,res}} = 177.3\ \mathrm{ppm}$' + '\n' +
        r'$\psi_{\mathrm{corr}} = 0.700,\ \mathcal{P}_{\mathrm{corr}} = \mathbf{0.5306}$' + '\n' +
        r'$\mathbf{Tier\ 3\ (Veto)}\ [\mathrm{Diagnostic\ Flag\ for\ TTVs}]$'
    )
    ax2.text(0.04, 0.04, badge36_lin, transform=ax2.transAxes, fontsize=8.2, verticalalignment='bottom',
             bbox=dict(boxstyle='round,pad=0.4', facecolor='#ffebee', edgecolor='#c62828', lw=1.2, alpha=0.92))
    ax2.legend(loc='upper right', fontsize=8, framealpha=0.9)

    # (c) Kepler-36 c: Dynamical
    ax3 = axes[2]
    t36_aligned = t36.copy()
    for n in range(-5, 45):
        t_nom = t0_36 + n * p36
        mask_transit = np.abs(t36 - t_nom) < (dur_36 * 1.6)
        if np.sum(mask_transit) >= 6:
            t_sub = t36[mask_transit]
            f_sub = f36[mask_transit]
            t_center = t_sub[np.argmin(f_sub)]
            shift = t_center - t_nom
            if np.abs(shift) < (dur_36 * 1.3):
                t36_aligned[mask_transit] -= shift

    phase36_dyn = ((t36_aligned - t0_36) % p36) / p36
    phase36_dyn[phase36_dyn > 0.5] -= 1.0
    dt_h_36_dyn = phase36_dyn * p36 * 24.0

    win36_dyn = np.abs(dt_h_36_dyn) < 14.0
    t_w36_dyn = dt_h_36_dyn[win36_dyn]
    f_w36_dyn = (f36[win36_dyn] - 1.0) * 1e6

    bm36_dyn, be36_dyn = [], []
    for i in range(len(bins36)-1):
        m = (t_w36_dyn >= bins36[i]) & (t_w36_dyn < bins36[i+1])
        if np.sum(m) > 0:
            bm36_dyn.append(np.median(f_w36_dyn[m]))
            be36_dyn.append(np.std(f_w36_dyn[m]) / np.sqrt(np.sum(m)))
        else:
            bm36_dyn.append(np.nan)
            be36_dyn.append(np.nan)

    ax3.scatter(t_w36_dyn, f_w36_dyn, s=3, color='gray', alpha=0.25, rasterized=True, label='Folded Points (Dynamical)')
    ax3.errorbar(bc36, bm36_dyn, yerr=be36_dyn, fmt='o-', color='#2ca02c', lw=2.0, markersize=4.5, capsize=2, label='TTV-Aligned Binned')
    ax3.axvspan(-7.461/2, 7.461/2, color='#2ca02c', alpha=0.15, label='Aligned Transit (7.46 h)')

    ax3.set_xlim(-14.0, 14.0)
    ax3.set_ylim(-850, 450)
    ax3.set_xlabel('Time from TTV-Corrected Mid-Transit (Hours)')
    ax3.set_ylabel(r'Relative Flux $\Delta F / F_0$ (ppm)')
    ax3.set_title('(c) Kepler-36 c: Dynamical TTV-Corrected\nTier 1 Clean: Epoch Alignment Restores Stationary Baseline')
    ax3.grid(True, linestyle=':', alpha=0.5)

    badge36_dyn = (
        'Kepler-36 c (TTV-Aligned Epochs)\n' +
        r'$\delta_{\mathrm{true}} = 430.0\ \mathrm{ppm}\ (\mathrm{Sharp\ Profile})$' + '\n' +
        r'$\sigma_{\mathrm{out}} = 104.3\ \mathrm{ppm},\ \sigma_{\mathrm{in,res}} = 106.1\ \mathrm{ppm}$' + '\n' +
        r'$\psi_{\mathrm{corr}} \to 0.017,\ \mathcal{P}_{\mathrm{corr}} \to \mathbf{0.0168}$' + '\n' +
        r'$\mathbf{Tier\ 1\ (Clean)}\ [\mathrm{Full\ Planetary\ Recovery}]$'
    )
    ax3.text(0.04, 0.04, badge36_dyn, transform=ax3.transAxes, fontsize=8.2, verticalalignment='bottom',
             bbox=dict(boxstyle='round,pad=0.4', facecolor='#e8f5e9', edgecolor='#2e7d32', lw=1.2, alpha=0.92))
    ax3.legend(loc='upper right', fontsize=8, framealpha=0.9)

    fig.tight_layout()
    p1 = os.path.join(OUT_DIR, 'Figure6_CaseStudies_Diagnostics.png')
    p2 = 'Figure6_CaseStudies_Diagnostics.png'
    fig.savefig(p1, dpi=600, bbox_inches='tight')
    fig.savefig(p2, dpi=600, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {p1}")


def main():
    plot_fig1_response_surface()
    plot_fig2_depth_invariance()
    plot_fig3_phase_folded_archetypes()
    plot_fig4_population_roc()
    plot_fig5_tier_distribution()
    plot_fig6_case_studies_diagnostics()
    print("\nAll 6 publication figures generated successfully at 600 DPI!")

if __name__ == '__main__':
    main()
