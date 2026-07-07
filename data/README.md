# Supplementary Data: Complete Benchmark Catalog (Table C4)

This directory contains the pre-processed diagnostic metrics and classification outcomes for the complete curated catalog of $N = 151$ target light curves. This dataset serves as the foundational data for validating the performance of **Log-SNRAS** against traditional and robust baselines.

## File Contents

### `Supplementary_Table_C4.csv`
A comma-separated values (CSV) file containing the calculated Signal-to-Noise Ratio (SNR) metrics and vetting outcomes for all $151$ segments in the benchmark.

### Column Descriptions:
1. **`ID`**: Cleaned target name formatted for LaTeX compatibility.
2. **`T_SNR`**: Traditional Signal-to-Noise Ratio (T-SNR).
3. **`R_SNR`**: Robust Signal-to-Noise Ratio (R-SNR) using Median Absolute Deviation (MAD).
4. **`P_SNR`**: Pont et al. red-noise corrected SNR proxy.
5. **`B_SNR`**: Box Least Squares (BLS) power proxy SNR.
6. **`L_SNRAS`**: Log-SNRAS (Proposed local-variance penalized metric).
7. **`Psi`**: Spatially-confined dispersion contrast ($\psi$) at the transit boundary.
8. **`Penalty`**: Adaptive variance penalty percentage applied to the traditional SNR.
9. **`Outcome`**: Vetting tier classification (`Confirmed`, `Flagged`, or `Vetoed`).
10. **`Rank`**: Numerical rank mapping the outcome classification (`1` = Confirmed, `2` = Flagged, `3` = Vetoed).
11. **`Runtime`**: Execution time in seconds for the Log-SNRAS algorithm.
12. **`RawFilename`**: Original FITS filename of the light curve.

---

## Reproducing Table 3 and ROC Curves

The script `reproduce_table3_roc.m` is designed to run directly using this dataset. It will load `Supplementary_Table_C4.csv`, perform the stratified bootstrap analysis ($B = 2000$), print the performance metrics, and generate the final ROC curve.

### How to run:
1. Ensure `reproduce_table3_roc.m` and `data/Supplementary_Table_C4.csv` are in your MATLAB working path.
2. Run the script:
   ```matlab
   reproduce_table3_roc
   ```
