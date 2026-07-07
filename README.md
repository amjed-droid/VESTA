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

## Repository Structure
To facilitate usage and reproducibility, the repository is organized as follows:
* `src/` : Contains the core standalone function (`calculate_log_snras.m`).
* `scripts/` : Contains MATLAB scripts to reproduce figures and tables.
* `data/` : Contains the curated evaluation catalog (`Supplementary_Table_C4.csv`).

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
2. **Required Toolboxes**: Statistics and Machine Learning Toolbox.
3. **Random Seed**: A fixed seed (`rng(42)`) is explicitly set in the scripts for the stratified bootstrap resampling ($B=2000$) to guarantee exact reproduction.
4. **Release Tag**: The results correspond to the `v1.1-revised` tag.

## Figure Generation & Scripts (Located in `/scripts`)
1. **`reproduce_table3_roc.m`**: Loads the pre-processed dataset, executes a stratified bootstrap ($B=2000$) to calculate Mean AUC and 95% CIs (Table 3), and generates the high-resolution empirical ROC curve (Figure 2).
2. **`generate_table2_stats.m`**: Extracts empirical class balance and median penalty values (Table 2).
3. **`plot_3d_surface.m`**: Regenerates the 3D response surface mapping theoretical penalty behavior (Figure 1).

## Data Preparation & Pre-processing
Light curves evaluated in the empirical benchmark were downloaded from MAST, utilizing the PDCSAP flux column. The in-transit window is defined as the interval $\pm T_{dur}/2$ centered on the BLS mid-transit time, while the out-of-transit window utilizes the remaining points after removing a $2 \times T_{dur}$ buffer on each side.

* **Note on NaNs**: NaN values in the dataset correspond to segments with insufficient in-transit points ($N_{in} < 10$) and were excluded to prevent undefined behavior.
* **Note on ROC curves**: All ROC curves (including T-SNR, R-SNR, P-SNR, B-SNR, and Log-SNRAS) are plotted using the real empirical data in the provided catalog `Supplementary_Table_C4.csv`.

## Citation
If you use Log-SNRAS in your research or pipeline, please cite our foundational paper:

```text
Jabbar, A. S. (2026). Log-SNRAS: A Computationally Efficient Variance-Stabilized Metric for Vetting Heteroscedastic Light Curves. Astronomy and Computing. (Under Review)
```

## License
This project is licensed under the MIT License - see the LICENSE.txt file for details.
