%% MAIN_PIPELINE Master orchestration script for Log-SNRAS manuscript.
%
%   MAIN_PIPELINE executes the complete end-to-end processing, analysis,
%   and figure-generation pipeline for the manuscript "Log-SNRAS: A
%   Computationally Efficient Variance-Stabilized Metric for Vetting
%   Heteroscedastic Light Curves."
%
%   EXECUTION FLOW
%   --------------
%   1. Configuration: Loads all ephemeris parameters, folder paths, tier
%      thresholds, and ground-truth keyword lists via CONFIG().
%   2. File Processing: Parses TESS (.fits) and BBS (.csv/.xlsx) light
%      curves, applies ephemeris/literature ground-truth labeling, computes
%      all candidate metrics (traditional/robust/PONT/BLS SNR and Log-SNRAS),
%      and exports the unified evaluation table to evaluation_dataset_v2.csv.
%   3. Dataset Statistics: Computes the median Psi and D for Confirmed and
%      Artifact classes and writes Table 1 to LaTeX.
%   4. Performance Analysis: Validates the discriminative power of each
%      metric via AUC, stratified bootstrap (B=2000), Jackknife, DeLong
%      exact tests, and permutation tests (10,000 shuffles), writing
%      Table 3 to LaTeX.
%   5. Visualization: Generates and exports the ROC comparison figure and
%      the bootstrap distribution histogram of the penalty-based AUC.
%
%   DEPENDENCIES
%   ------------
%   The following files must be present in the same directory:
%       config.m, helpers.m, process_files.m, generate_table1.m,
%       analyze_performance.m, plot_results.m, calculate_log_snras.m
%   and subfolders ./TESS/ and ./BBS/ containing the raw light-curve data.
%
%   OUTPUTS
%   -------
%   Evaluation catalog:       evaluation_dataset_v2.csv
%   Tables (LaTeX):           Table1_DatasetStats.tex, Table3_Performance.tex
%   Figures (600 DPI PNG):    ROC_EphemerisGroundTruth.png,
%                             Bootstrap_PenaltyDistribution.png
%
%   AUTHOR
%   ------
%   Ahmed Sattar Jabbar.
clear; clc; close all;

% 1. Load configurations
cfg = config();

% 2. Process files and generate evaluation_dataset_v2.csv
results = process_files(cfg);

% 3. Generate Table 1 (LaTeX and Command Window output)
generate_table1(results);

% 4. Analyze performance, compute Bootstrap, and statistical tests
[y_true_v, scores_inv, raw_auc, boot_pen, auc_penalty, n_valid, raw_scores_cell] = analyze_performance(results);

% 5. Plot ROC and Bootstrap distribution figures
plot_results(y_true_v, scores_inv, raw_auc, boot_pen, auc_penalty, n_valid, raw_scores_cell);

fprintf('\n=== DONE. All tables and figures regenerated from clean ground truth. ===\n');
