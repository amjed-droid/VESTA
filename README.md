# Log-SNRAS
**A Computationally Efficient Variance-Stabilized Metric for Vetting Heteroscedastic Light Curves**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2022b%2B-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![Release](https://img.shields.io/badge/Release-v2.0--overhaul-success.svg)](#)

## Overview
**Log-SNRAS** (Logarithmic Signal-to-Noise Ratio with Adjusted Statistics) is a computationally efficient, non-parametric **post-detection vetting filter** designed to automate the classification of exoplanetary transit candidates from space-based missions (*Kepler*, *K2*, *TESS*).

Standard SNR metrics and search algorithms (e.g., BLS, TPS) evaluate signal significance assuming stationary, homoscedastic noise. In real time-series photometry, instrumental anomalies, stellar variability, and non-stationary artifacts frequently induce false-positive detections. Log-SNRAS penalizes heteroscedastic candidates by evaluating the **relative dispersion contrast ($\psi$)** between in-transit and out-of-transit flux, scaling detection significance through a physically-motivated logarithmic penalty:

$$\text{Log-SNRAS} = \frac{\delta}{\hat{\sigma}_{\text{out}}} \times \frac{\sqrt{N_{\text{in}}}}{1 + \ln(1 + \psi_{\text{corr}})}$$

---

## Key Breakthrough: Shape-Corrected Dispersion ($\sigma_{\text{in, res}}$)
In earlier single-host formulations, deep planetary transits were susceptible to false variance penalties because deterministic geometric transit curvature (ingress, flat bottom, egress) inflated the raw in-transit standard deviation $\sigma_{\text{in}}$.

Log-SNRAS resolves this by decoupling geometric transit curvature from stochastic photometric noise:
$$\sigma_{\text{in, res}} = \sqrt{\frac{1}{N_{\text{in}} - 1} \sum_{i \in \text{in}} (f_i - \hat{f}_i)^2}$$
where $\hat{f}$ is a non-parametric transit profile fit in orbital phase space. 

This correction completely eliminates false vetoes on deep planetary transits:
* **TOI-201** ($R_p = 1.01 R_J$, depth $\sim 5{,}200$ ppm): Penalty drops from $\mathcal{P}_{\text{raw}} = 1.640$ (was Tier 3 Veto) to $\mathcal{P}_{\text{corr}} = 0.498$ (**Restored to Clean Planetary Transit**).
* **Kepler-8 b** (depth $\sim 10{,}000$ ppm): Penalty drops from $0.951$ to **$0.004$** (**Tier 1 Clean**).
* **Kepler-11 d** (depth $\sim 10{,}000$ ppm): Penalty drops from $1.552$ to **$0.031$** (**Tier 1 Clean**).
* **Deep Transits from TESS & Kepler** (tested on $\delta$ up to $56\%$): Reliably classified as **Tier 1 Clean**.

---

## Curated Multi-Host Benchmark ($N=17$) and LOHO Validation
The method is validated across a balanced, literature-verified multi-host benchmark:
* **10 Confirmed Exoplanets across 10 Distinct Hosts**: Pi Mensae c, TOI-201 b, Kepler-10 b, Kepler-448 b, TrES-2 b, WASP-126 b, TOI-540 b, L 98-59 b, Kepler-8 b, and Kepler-11 d.
* **7 Literature-Proven Non-Planetary Artifacts**: Boyajian's Star (Quarters 8 & 16 anomalous dips), KOI-1611 eccentric eclipsing binary, Kepler eccentric EB, contact binaries (W UMa KIC 11253226), and ellipsoidal variables.

### Empirical Classification Performance:
| Method | AUC | 95% Bootstrap CI | Mean LOHO AUC |
| :--- | :---: | :---: | :---: |
| **Shape-Corrected Penalty ($-\mathcal{P}_{\text{corr}}$)** | **0.757** | **[0.486, 0.971]** | **0.757 $\pm$ 0.036** |
| Log-SNRAS Composite | 0.700 | [0.414, 0.943] | 0.700 $\pm$ 0.035 |
| Traditional SNR (T-SNR) | 0.629 | [0.321, 0.914] | 0.629 $\pm$ 0.032 |
| Robust SNR (MAD) | 0.614 | [0.314, 0.886] | 0.614 $\pm$ 0.033 |
| BLS SNR Proxy | 0.614 | [0.314, 0.886] | 0.614 $\pm$ 0.033 |
| Pont SNR (2006) | 0.600 | [0.271, 0.886] | 0.600 $\pm$ 0.034 |
| Inverse Depth ($1/\delta$) | 0.471 | [0.143, 0.857] | 0.471 $\pm$ 0.045 |

**Leave-One-Host-Out (LOHO) Cross-Validation**:
* Excluding Pi Mensae yields **AUC = 0.762**.
* Excluding TOI-201 yields **AUC = 0.778**.
* Mean LOHO AUC = **$0.757 \pm 0.036$**, confirming that performance is host-invariant.

---

## Repository Structure
```
Log-SNRAS/
├── src/
│   ├── log_snras/                 # Modular Python package
│   │   ├── __init__.py
│   │   ├── core.py               # Log-SNRAS math, shape correction, tiers
│   │   ├── masking.py            # Uniform ephemeris-based window definition
│   │   └── io.py                 # FITS named-column reader (PDCSAP_FLUX)
│   └── calculate_log_snras.m     # Standalone vectorized MATLAB function
├── data/
│   ├── curated_benchmark_catalog.csv    # 17 multi-host benchmark targets
│   ├── multi_host_evaluation_results.csv # Empirical metric evaluation
│   ├── multi_host_auc_performance.csv   # Statistical AUC comparison
│   └── evaluation_dataset_v2.csv        # Archival dataset (superseded)
├── scripts/
│   ├── evaluate_multi_host_benchmark.py # Automated end-to-end evaluation
│   ├── generate_publication_figures.py  # 100% real archival FITS figures
│   ├── test_nrebig.py                   # Empirical validation on TESS targets
│   └── test_nrebig_kepler.py            # Empirical validation on Kepler targets
└── tests/
    └── test_core.py                     # Unit tests for shape correction and tiers
```

---

## Quickstart (Python)
```bash
# 1. Install dependencies
pip install numpy scipy astropy pandas matplotlib scikit-learn

# 2. Run automated multi-host evaluation
python scripts/evaluate_multi_host_benchmark.py
```

### Python API Example:
```python
from log_snras.core import compute_log_snras, compute_shape_corrected_dispersion

# time, flux, and in_mask (True = In-Transit)
# Compute shape-corrected in-transit dispersion
sig_in_res, fitted_profile = compute_shape_corrected_dispersion(time, flux, in_mask)

# Compute Log-SNRAS and vetting tier
result = compute_log_snras(flux, in_mask, out_mask, snr_trad=45.0, sig_in_res=sig_in_res)
print(f"Log-SNRAS Score: {result.log_snras:.2f}")
print(f"Vetting Tier: {result.tier} (Penalty = {result.penalty:.4f})")
```

---

## Quickstart (MATLAB)
```matlab
% Add src to path
addpath('src');

% Calculate Log-SNRAS
[score, penalty, psi, tier] = calculate_log_snras(flux, in_mask, out_mask, snr_trad);
fprintf('Tier: %s | Log-SNRAS: %.2f\n', tier, score);
```

---

## Reproducibility Statement
Every table and figure in the manuscript is generated from 100% real NASA archival Kepler and TESS light curves using explicit named columns (`PDCSAP_FLUX`), uniform orbital ephemeris masking, and deterministic seeds. Run:
```bash
python scripts/generate_publication_figures.py
```
to regenerate all publication-grade figures (600 DPI) into `figures/`.

---

## Citation
```bibtex
@article{jabbar2026logsnras,
  author  = {Jabbar, Ahmed Sattar},
  title   = {Log-SNRAS: A Computationally Efficient Variance-Stabilized Metric for Vetting Heteroscedastic Light Curves},
  journal = {Astronomy and Computing},
  year    = {2026}
}
```

## License
MIT License. See [LICENSE](LICENSE) for details.
