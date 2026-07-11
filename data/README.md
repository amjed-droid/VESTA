# Supplementary Data: Complete Evaluation Catalog

This directory contains the unified and pre-processed evaluation dataset used to validate the performance of the proposed **Log-SNRAS** algorithm against traditional and robust signal-to-noise ratio (SNR) baselines.

## File Contents

### `evaluation_dataset_v2.csv`
A comma-separated values (CSV) file containing the calculated Signal-to-Noise Ratio (SNR) metrics, spatial dispersion coefficients, morphological properties, and vetting outcomes for all $N = 180$ light curve segments in the benchmark catalog.

---

## Column Descriptions

1. **`Filename`**: The original filename of the light curve segment (FITS or CSV).
2. **`SourceType`**: The data source type (e.g., `TESS_fits` or `BBS_csv`).
3. **`Label`**: Ground truth diagnostic classification label (`1` = confirmed transit, `0` = false positive/literature artifact, `NaN` = complex/excluded).
4. **`Label_Source`**: The source of the ground truth label (e.g., `literature_confirmed`, `literature_artifact`, `ephemeris_confirmed`, `complex_excluded`, `unknown_excluded`).
5. **`N_total`**: Total number of data points (measurements) in the light curve.
6. **`N_in`**: Number of data points falling inside the expected transit window.
7. **`N_out`**: Number of data points falling outside the transit window.
8. **`Depth_ppm`**: Transit depth in parts per million (ppm).
9. **`T_SNR`**: Traditional Signal-to-Noise Ratio (T-SNR).
10. **`R_SNR`**: Robust Signal-to-Noise Ratio (R-SNR) using Median Absolute Deviation (MAD).
11. **`P_SNR`**: Pont et al. red-noise corrected SNR proxy.
12. **`B_SNR`**: Box Least Squares (BLS) power proxy SNR.
13. **`L_SNRAS`**: Proposed Log-SNRAS (local-variance penalized metric).
14. **`Psi`**: Spatially-confined dispersion contrast ($\psi$) at the transit boundary.
15. **`Penalty_pct`**: Adaptive variance penalty percentage applied to the traditional SNR.
16. **`Suppression_pct`**: The signal suppression percentage under correction.
17. **`Tier`**: Classified vetting tier outcome (`Tier 1`, `Tier 2`, or `Tier 3`).
18. **`Transit_Expected`**: Binary indicator (0 or 1) showing whether a transit event is expected.
19. **`T_start_BTJD`**: Observation start time in Barycentric TESS Julian Date (BTJD).
20. **`T_end_BTJD`**: Observation end time in BTJD.

---

## Reproducing Table 3 and ROC Curves

The validation script `reproduce_table3_roc.m` in the repository's main directory is designed to run directly using this dataset. It loads the dataset, performs stratified bootstrap analysis ($B = 2000$), computes confidence intervals, prints performance summaries, and generates the final ROC curve.

### How to Run:
1. Ensure `reproduce_table3_roc.m` and `Data/evaluation_dataset_v2.csv` are in your MATLAB working path.
2. In the script `reproduce_table3_roc.m`, make sure the file path points to the new dataset:
   ```matlab
   csv_file = 'Data/evaluation_dataset_v2.csv';
   ```
3. Run the script in the MATLAB Command Window:
   ```matlab
   reproduce_table3_roc
   ```
