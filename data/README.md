# VESTA Benchmark Datasets and Evaluation Catalogs

This directory contains the observational catalogs, empirical results, and diagnostic tables used to evaluate **VESTA** (**V**ariance-**E**valuated **S**tability for **T**ransit **A**nalysis) across four orders of magnitude in transit depth.

---

## 1. Primary Population-Scale Benchmark Catalog ($N = 56$)

### `kepler_population_catalog.csv`
The core archival catalog establishing population-scale validation across four full decades of transit depth ($12$~ppm to $32{,}040$~ppm), matched with official NASA Mikulski Archive for Space Telescopes (MAST) *Kepler* FITS products.
* **49 Confirmed *Kepler* Exoplanets** across 4 observational depth strata:
  1. **Ultra-Deep Transits (>10,000 ppm, $N=5$):** Giant exoplanets such as Kepler-45 b ($32{,}040$~ppm), Kepler-468 b ($23{,}530$~ppm), Kepler-548 b ($16{,}040$~ppm), Kepler-762 b ($13{,}890$~ppm), and Kepler-148 d ($10{,}310$~ppm).
  2. **Moderate Transits (1,000–10,000 ppm, $N=24$):** Sub-Jupiter and Neptune-sized worlds, including Kepler-702 b, Kepler-39 b, Kepler-43 b, Kepler-46 b, Kepler-117 c, and Kepler-42 b.
  3. **Shallow Super-Earths (100–1,000 ppm, $N=13$):** Small rocky and volatile-rich planets, including Kepler-79 c ($790$~ppm), Kepler-11 b ($300$~ppm), Kepler-78 b ($220$~ppm), and Kepler-10 b ($190$~ppm).
  4. **Ultra-Shallow Sub-Earths (<100 ppm, $N=7$):** Terrestrial candidates at the photon-noise limit, including Kepler-62 c ($70$~ppm), Kepler-444 e ($51$~ppm), and Kepler-37 b ($12$~ppm).
* **7 Certified Non-Planetary Astrophysical Variables and Artifacts**:
  * Boyajian's Star (KIC 8462852, Quarter 8 D792 deep anomalous dip)
  * Boyajian's Star (KIC 8462852, Quarter 16 D1500 complex dip)
  * KOI-1611 EB (KIC 4544587, eccentric eclipsing binary)
  * KIC 3858884 (pulsating eccentric binary)
  * V380 Cyg (KIC 5385723, eccentric heartbeat binary)
  * KIC 11253226 (W UMa contact binary)
  * V723 Mon (ellipsoidal variable)

Every target's physical classification and ephemerides ($T_0, P, T_{\text{dur}}$) are verified against the NASA Exoplanet Archive cumulative table and peer-reviewed literature.

### `kepler_population_results.csv`
Contains the complete evaluation metrics across all 56 systems processed through the VESTA pipeline:
* Multi-quarter FITS light curve stitching from `D:\nrebig\kepler\lightcurves\`.
* Ephemeris-based transit window definition with a $1.5\times T_{\text{dur}}$ buffer.
* Non-parametric transit shape detrending ($\sigma_{\text{in, res}}$) and residual dispersion calculation.
* Computed stability penalty ($\mathcal{P}_{\text{corr}}$), dispersion contrast ($\psi$), and composite VESTA detection scores.
* Comparator baseline metrics: Traditional matched-filter SNR ($\text{SNR}_{\text{trad}}$), Robust MAD SNR ($\text{SNR}_{\text{MAD}}$), Pont Red-Noise SNR ($\text{SNR}_{\text{Pont}}$), BLS SNR Proxy ($\text{SNR}_{\text{BLS}}$), and Inverse Depth ($1/\delta$).
* Calibrated statistical tier assignments (`Tier 1 Clean`, `Tier 2 Review`, `Tier 3 Veto`).

### `kepler_population_performance.csv`
Statistical performance summary reporting:
* Area Under the ROC Curve (AUC) for all metrics.
* Stratified 95% bootstrap confidence intervals ($B = 2000$).
* DeLong paired statistical significance tests comparing VESTA against each baseline.

### `latex_catalog_table_clean.tex`
LaTeX-formatted reference catalog table presenting all 56 benchmark systems, designed for direct inclusion in journal manuscripts.

---

## 2. Legacy Catalogs (Retained for Historical Provenance)

### `curated_benchmark_catalog.csv`
* **Description**: Preliminary multi-host proof-of-concept catalog ($N = 17$) developed during earlier revisions (R3–R4).
* **Status**: Superseded by the population-scale 56-system benchmark (`kepler_population_catalog.csv`).

### `evaluation_dataset_v2.csv`
* **Description**: Archival 180-segment single-host catalog (Pi Mensae c) from original submission (R0–R2).
* **Status**: Deprecated due to host-star concentration bias.

---

## 3. Reproduction Pipeline

To evaluate the full population benchmark and reproduce all catalog results:
```bash
python scripts/evaluate_kepler_population.py
```
Execution takes approximately $34$~seconds on a standard CPU workstation for all 56 multi-quarter systems.
