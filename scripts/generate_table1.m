%% ==========================================================
%  Standalone script to compute Table 1 statistics 
%  ==========================================================
%  This script reads the CSV file created by rocess_files.m
%  and prints the dataset statistics to the Command Window.
%  ==========================================================
clear; clc;

% 1. Load the evaluation dataset
csv_file = 'evaluation_dataset_v2.csv';
if ~isfile(csv_file)
    error('File %s not found. Please ensure it exists in the current directory.', csv_file);
end
results = readtable(csv_file);

% 2. Compute statistics (exactly as in master_pipeline.m)
idx_confirmed = results.Label == 1;
idx_artifact  = results.Label == 0;
n_confirmed = sum(idx_confirmed);
n_artifact  = sum(idx_artifact);
n_excluded  = sum(isnan(results.Label));

med_psi_confirmed = median(results.Psi(idx_confirmed), 'omitnan');
med_D_confirmed   = median(1 + results.Penalty_pct(idx_confirmed)/100, 'omitnan');
med_psi_artifact  = median(results.Psi(idx_artifact), 'omitnan');
med_D_artifact    = median(1 + results.Penalty_pct(idx_artifact)/100, 'omitnan');

% 3. Print results directly to the Command Window
fprintf('\n=====================================================\n');
fprintf('              TABLE 1: DATASET STATISTICS\n');
fprintf('=====================================================\n');
fprintf('Class                  N      Median Psi      Median D\n');
fprintf('-----------------------------------------------------\n');
fprintf('Confirmed              %-5d   %-12.3f   %-12.3f\n', n_confirmed, med_psi_confirmed, med_D_confirmed);
fprintf('Artifact               %-5d   %-12.3f   %-12.3f\n', n_artifact, med_psi_artifact, med_D_artifact);
fprintf('-----------------------------------------------------\n');
fprintf('Excluded (complex/ambiguous) : N = %d\n', n_excluded);
fprintf('=====================================================\n');
