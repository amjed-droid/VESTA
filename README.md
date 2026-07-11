# Log-SNRAS
**A Computationally Efficient Variance-Stabilized Metric for Vetting Heteroscedastic Light Curves**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Language](https://img.shields.io/badge/MATLAB-R2025b%2B-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![Release](https://img.shields.io/badge/Release-v1.1--revised-success.svg)](#)

## Overview
**Log-SNRAS** (Logarithmic Signal-to-Noise Ratio with Adjusted Statistics) is a computationally efficient **post-detection vetting filter** designed to automate the evaluation of exoplanet candidates in space-based photometry (e.g., *Kepler*, *TESS*).

Unlike traditional search algorithms (like BLS or TPS) or standard SNR metrics which assume homoscedastic (uniform) noise, Log-SNRAS operates *after* detection. It introduces a physically motivated logarithmic penalty term based on the **dispersion contrast ($\psi$)** between in-transit and out-of-transit flux. This allows it to penalize instrumental artifacts and non-stationary noise events (heteroscedasticity) that often mimic planetary signals.

**Formula:**
$$\text{Log-SNRAS} = \frac{\delta}{\sigma_{out}} \times \frac{\sqrt{N_{in}}}{1 + \ln(1 + \psi)}$$

## Known Limitation: Host-Star Confound
Our validation benchmark's Confirmed class is currently dominated by a single host star (TIC 261136679 / Pi Mensae). A leave-one-host-out analysis shows the reported discriminative performance has **not yet been established across independent host stars** (AUC drops to 0.524 on a host-diverse subset, N=22). We further traced one such misclassification (TOI-201) to a mechanistic cause: for deep transits, the transit's own geometric shape inflates the in-transit dispersion independently of genuine noise. See Section 3.5 of the manuscript, and `scripts/plots/generate_figure7_toi201_diagnostic.m`, for the full analysis and diagnostic figure. We disclose this prominently because it defines the current operational boundary of the method, and identifies a concrete direction (shape-corrected variance estimation) for future work.

## Repository Structure
To facilitate usage and reproducibility, the repository is organized as follows:
* `src/` : Contains the core standalone function (`calculate_log_snras.m`).
* `data/` : Contains the curated evaluation catalog (`evaluation_dataset_v2.csv`).
* `docs/` : Contains documentation and validation files.
* `scripts/` : Organized subdirectories containing MATLAB scripts:
  * `pipeline/` : Entrypoint orchestration (`main_pipeline.m`) and file ingestion pipeline (`process_files.m`).
  * `tables/` : Scripts (`generate_table1.m` through `generate_table6.m`) to generate each respective table's data and LaTeX/plain-text output.
  * `plots/` : Scripts to generate the ROC curves, response surfaces, and diagnostic figures, including `generate_figure7_toi201_diagnostic.m` (host-star confound diagnostic, Section 3.5).
  * `utils/` : Shared configurations, helpers, and statistical functions (DeLong exact tests, bootstrap).

## Installation
The code requires no compilation. You only need a standard MATLAB environment.
1. Clone the repository:
```bash
   git clone https://github.com/amjed-droid/Log-SNRAS.git
```
2. Add the `src/` and `scripts/` directories to your MATLAB path:
```matlab
   addpath(genpath('Log-SNRAS/src'));
   addpath(genpath('Log-SNRAS/scripts'));
```

## Usage
The core function `calculate_log_snras.m` operates in $\mathcal{O}(N)$ linear time.
```matlab
% 1. Load your normalized light curve array and Boolean mask (True = In-Transit)
flux = ...; 
mask = ...; 

% 2. Provide the traditional SNR value computed by your pipeline (e.g., BLS SNR)
snr_trad = ...; 

% 3. Run Log-SNRAS
[score, penalty_factor, psi] = calculate_log_snras(flux, snr_trad, mask);

fprintf('Log-SNRAS Score: %.2f\n', score);
fprintf('Denominator Scaling Factor (D): %.2f\n', penalty_factor);
```

## Reproducibility Statement
A major strength of this work is its complete reproducibility. This repository contains all necessary scripts, exact random seeds, and data to fully reproduce the results, tables, and figures presented in our foundational manuscript.

1. **Environment**: Developed and tested using MATLAB R2025b.
2. **Required Toolboxes**: Statistics and Machine Learning Toolbox; Curve Fitting Toolbox (required only for `generate_figure7_toi201_diagnostic.m`, which uses LOESS smoothing).
3. **Random Seed**: A fixed seed (`rng(42)`) is explicitly set in the scripts for the stratified bootstrap resampling ($B=2000$) to guarantee exact reproduction.
4. **Release Tag**: The results correspond to the `v1.1-revised` tag.

## Figure & Table Generation (Located in `/scripts`)
All figures and tables presented in the manuscript can be reproduced using the scripts under the `/scripts` directory:
1. **Tables (under `/scripts/tables`)**:
   * `generate_table1.m` through `generate_table6.m` to generate each respective table's data and LaTeX/plain-text output.
2. **Figures (under `/scripts/plots`)**:
   * `plot_results.m` and `generate_baseline_figure.m` to generate the empirical ROC curves and comparison figures.
   * `plot_3d_surface.m` to regenerate the 3D response surface mapping the theoretical penalty behavior.
   * `generate_scatter_plot.m` to create the metric stability mapping.
   * `generate_figure7_toi201_diagnostic.m` to reproduce the host-star confound diagnostic for TOI-201 (Section 3.5), isolating the transit's geometric shape from genuine photometric noise.

## Data Preparation & Pre-processing
Light curves evaluated in the empirical benchmark were downloaded from MAST, utilizing the PDCSAP flux column. The in-transit window is defined as the interval $\pm T_{dur}/2$ centered on the mid-transit time: for TIC 261136679 (Pi Mensae), the published ephemeris of Huang et al. (2018) is used; for all other targets, the BLS-derived mid-transit time is used. The out-of-transit window comprises the remaining points after removing a $2 \times T_{dur}$ buffer on each side.

* **Note on NaNs**: NaN values in the dataset correspond to segments with insufficient in-transit points ($N_{in} < 30$) and were excluded from binary AUC evaluation to ensure a statistically reliable estimate of $\sigma_{in}$ (see Appendix B of the manuscript).
* **Note on ROC curves**: All ROC curves (including T-SNR, R-SNR, P-SNR, B-SNR, and Log-SNRAS) are plotted using the real empirical data in the provided catalog `evaluation_dataset_v2.csv`.
* **Note on ground truth**: Confirmed/Artifact labels are derived independently of the metric's own output — from the published ephemeris (TIC 261136679) or from independent literature identity (all other targets) — never from the penalty, $\psi$, or Tier assignment itself. See Section 3.1 of the manuscript and the `Label_Source` column in `evaluation_dataset_v2.csv` for full provenance.

## Citation
If you use Log-SNRAS in your research or pipeline, please cite our foundational paper:

```text
Jabbar, A. S. (2026). Log-SNRAS: A Computationally Efficient Variance-Stabilized Metric for Vetting Heteroscedastic Light Curves. Astronomy and Computing. (Under Review)
```

## License
This project is licensed under the MIT License - see the LICENSE.txt file for details.
