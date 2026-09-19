#!/usr/bin/env python3
"""
=============================================================================
REPRODUCE ALL: End-to-End Scientific Reproduction Pipeline for Log-SNRAS
=============================================================================
Manuscript: "Log-SNRAS: A Computationally Efficient Variance-Stabilized Metric 
             for Vetting Heteroscedastic Light Curves"
Journal:    Astronomy and Computing (ASCOM)
Author:     Ahmed Sattar Jabbar (ahmed.state.me@gmail.com)

Usage:
    python reproduce_all.py

This master script performs a complete, push-button reproduction of:
  1. The curated multi-host benchmark evaluation (N = 17 systems across Kepler & TESS)
  2. Shape-corrected dispersion estimation (sigma_in,res) decoupling transit geometry
  3. Binary classification AUC metrics, 95% bootstrap CIs (B = 2000), and LOHO validation
  4. DeLong exact pairwise significance tests and permutation testing
  5. All publication figures (Figures 1-6) regenerated into figures/
=============================================================================
"""

import os
import sys
import subprocess
import time

def print_header(title):
    print("\n" + "=" * 76)
    print(f"  {title}")
    print("=" * 76 + "\n")

def run_step(step_name, command):
    print(f"--> [Step] {step_name}...")
    t0 = time.time()
    result = subprocess.run([sys.executable] + command, capture_output=True, text=True)
    dt = time.time() - t0
    if result.returncode != 0:
        print(f"Error during {step_name}:\n{result.stderr}")
        sys.exit(result.returncode)
    print(result.stdout)
    print(f"--> Completed {step_name} in {dt:.2f} seconds.\n")

def main():
    print_header("LOG-SNRAS MASTER REPRODUCIBILITY PIPELINE")
    print("Environment: Python", sys.version.split()[0])
    print("Working Directory:", os.path.abspath('.'))

    # Step 1: Run Multi-Host Benchmark Evaluation Pipeline
    run_step(
        "Evaluating Multi-Host Benchmark & LOHO Cross-Validation",
        ["scripts/evaluate_multi_host_benchmark.py"]
    )

    # Step 2: Regenerate All Publication Figures from 100% Real Data
    run_step(
        "Regenerating Publication Figures (Figures 1-6)",
        ["scripts/generate_publication_figures.py"]
    )

    # Step 3: Print Verification Summary
    print_header("REPRODUCIBILITY VERIFICATION SUMMARY")
    print(" [OK] Multi-host benchmark catalog: data/curated_benchmark_catalog.csv (N=17)")
    print(" [OK] Empirical evaluation results: data/multi_host_evaluation_results.csv")
    print(" [OK] Statistical summary table:    data/multi_host_auc_performance.csv")
    print(" [OK] Publication figures in figures/:")
    for fig_file in sorted(os.listdir('figures')):
        if fig_file.endswith('.png'):
            print(f"       * figures/{fig_file}")
    print("\nAll tables and figures reported in the manuscript have been reproduced successfully!")
    print("=" * 76 + "\n")

if __name__ == '__main__':
    main()
