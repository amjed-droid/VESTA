# Benchmark Datasets and Evaluation Catalogs

This directory contains the data catalogs used to benchmark and evaluate **Log-SNRAS** against standard Signal-to-Noise Ratio (SNR) baselines.

---

## 1. Curated Multi-Host Benchmark Catalog (Current Standard)
### `curated_benchmark_catalog.csv`
The primary benchmark catalog introduced in the revised manuscript to resolve the single-host confound and demonstrate host-invariant performance. It comprises **$N = 17$ distinct astrophysical systems**:
* **10 Confirmed Exoplanets across 10 Distinct Hosts**:
  * Pi Mensae c (TIC 261136679, S1)
  * TOI-201 b (TIC 350618622, S4)
  * Kepler-10 b (KIC 11904151, Q1)
  * Kepler-448 b (KIC 5812701, Q1) — restored from literature
  * TrES-2 b (KIC 11446443, Q1) — restored from literature
  * WASP-126 b (TIC 25155310, S1)
  * TOI-540 b (TIC 200322593, S6)
  * L 98-59 b (TIC 307210830, S2)
  * Kepler-8 b (KIC 6922244, Q2)
  * Kepler-11 d (KIC 6541920, Q3)
* **7 Literature-Proven Non-Planetary Artifacts**:
  * Boyajian's Star (KIC 8462852, Q8 D792 deep anomalous dip)
  * Boyajian's Star (KIC 8462852, Q16 D1500 complex dip)
  * KOI-1611 (KIC 4544587, eccentric eclipsing binary)
  * Kepler eccentric EB (KIC 12644769)
  * W UMa contact binary (KIC 11253226)
  * Kepler ellipsoidal variable (KIC 8016214)
  * Kepler non-stationary EB (KIC 10001167)

*Every target's ground-truth classification is verified against the NASA Exoplanet Archive and peer-reviewed literature.*

### `multi_host_evaluation_results.csv`
Contains the calculated metrics across all 17 benchmark targets, including:
* Ingested from NASA MAST using explicit named column access (`PDCSAP_FLUX`).
* Uniform orbital ephemeris transit masking with a $2\times T_{\text{dur}}$ buffer.
* Shape-corrected in-transit dispersion ($\sigma_{\text{in, res}}$).
* Comparative metrics: Traditional SNR, Robust MAD SNR, Pont SNR, BLS SNR proxy, and Inverse Depth ($1/\delta$).
* Classified tiers (`Tier 1 Clean`, `Tier 2 Ambiguous`, `Tier 3 Veto`).

### `multi_host_auc_performance.csv`
Summary table reporting point AUCs, 95% stratified bootstrap confidence intervals ($B = 2000$), and Leave-One-Host-Out (LOHO) cross-validation scores.

---

## 2. Archival Catalog (Legacy)
### `evaluation_dataset_v2.csv`
* **Status**: Deprecated / Retained for Historical Provenance.
* **Description**: This file represents the evaluation table from earlier revision rounds (R1–R4). It contains 180 light curve segments, of which 76 confirmed transit segments were derived from a single quiet host star (TIC 261136679 / Pi Mensae). As noted by Reviewer #5, this single-host concentration produced an apparent AUC of 0.960 that dropped when tested across diverse hosts.
* **Resolution**: Replaced by `curated_benchmark_catalog.csv` ($N = 17$) and the shape-corrected metric, which achieves a host-invariant AUC of **0.757 $\pm$ 0.036**.

---

## Reproduction
To evaluate the multi-host benchmark and regenerate the results:
```bash
python scripts/evaluate_multi_host_benchmark.py
```
