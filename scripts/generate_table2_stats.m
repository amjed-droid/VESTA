%% Script to extract empirical statistics for Table 2
%  (Refactored to use Supplementary_Table_C4.csv for full reproducibility)
clc; clear;

% 1. Load the evaluation catalog (already included in the repository)
data = readtable('../data/Supplementary_Table_C4.csv');
data.Properties.VariableNames = strtrim(data.Properties.VariableNames);

% 2. Filter out segments excluded from final analysis (N_in < 10)
%    The README states that NaN in L_SNRAS means insufficient in-transit points.
%    We remove those rows to match the manuscript's evaluation.
valid_idx = ~isnan(data.L_SNRAS);
data_clean = data(valid_idx, :);

fprintf('Total rows in CSV: %d\n', height(data));
fprintf('Rows after removing NaNs (N_in < 10): %d\n', height(data_clean));

% 3. Label each segment as:
%      1 = Confirmed Planet
%      0 = Heteroscedastic Artifact
%    Logic replicated from the original FITS-based script.
nRows = height(data_clean);
all_labels = zeros(nRows, 1);
tic_count = 0;   % Counter for host star TIC 261136679

for i = 1:nRows
    name = lower(data_clean.Filename{i});  % case-insensitive matching
    
    if contains(name, 'planet')   || contains(name, 'confirmed') || ...
       contains(name, '261136679') || contains(name, '1625')
        all_labels(i) = 1;        % planet
        if contains(name, '261136679')
            tic_count = tic_count + 1;
        end
    else
        all_labels(i) = 0;        % artifact
    end
end

% 4. Extract median Psi and Penalty for each class
idx_planets   = (all_labels == 1);
idx_artifacts = (all_labels == 0);

num_planets   = sum(idx_planets);
num_artifacts = sum(idx_artifacts);
tic_percentage = round((tic_count / num_planets) * 100);

med_psi_planets   = median(data_clean.Psi(idx_planets));
med_pen_planets   = median(data_clean.Penalty_Val(idx_planets));

med_psi_artifacts   = median(data_clean.Psi(idx_artifacts));
med_pen_artifacts   = median(data_clean.Penalty_Val(idx_artifacts));

% 5. Print results ready to be copied into LaTeX / manuscript
fprintf('\n--- Statistics for Table 2 ---\n');
fprintf('Confirmed Planets: N = %d | Median Psi = %.3f | Median Penalty = %.3f\n', ...
        num_planets, med_psi_planets, med_pen_planets);
fprintf('Heteroscedastic Artifacts: N = %d | Median Psi = %.3f | Median Penalty = %.3f\n', ...
        num_artifacts, med_psi_artifacts, med_pen_artifacts);
fprintf('Note: %d%% of the Confirmed segments originate from TIC 261136679.\n', tic_percentage);
fprintf('=========================================================\n');