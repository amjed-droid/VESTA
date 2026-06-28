% =========================================================================
% Script: Stratified Bootstrap AUC and Confidence Intervals Validation
% Author: Ahmed S. Jabbar
% Description: Performs B=2000 bootstrap iterations to calculate the mean
% Area Under the Curve (AUC) and 95% Confidence Intervals for Log-SNRAS
% and baseline metrics.
% =========================================================================
clear; clc; close all;

disp('Loading dataset...');
% Ensure the Supplementary_Table_C4 dataset is in the current working directory
data = readtable('../data/Supplementary_Table_C4.csv');

% Clean column names
data.Properties.VariableNames = strtrim(data.Properties.VariableNames);

% 1. Create Ground Truth Labels
% 1 = Confirmed Planet, 0 = Artifact
n_rows = height(data);
y_true = nan(n_rows, 1);

for i = 1:n_rows
    f_name = upper(data.Filename{i});
    
    % Identify confirmed planetary transits
    if contains(f_name, '261136679') || contains(f_name, 'TOI_201') || ...
       contains(f_name, 'KEPLER_1625') || contains(f_name, 'TESS_S')
        y_true(i) = 1; 
        
    % Identify heteroscedastic artifacts / false positives
    elseif contains(f_name, '12644769') || contains(f_name, '5812701') || ...
           contains(f_name, 'V723_MON') || contains(f_name, 'TIC14444029')
        y_true(i) = 0; 
    end
end

% Filter out unlabeled rows
valid_idx = ~isnan(y_true);
y_true_clean = y_true(valid_idx);
data_clean = data(valid_idx, :);

% 2. Setup Evaluation Parameters
score_columns = {'T_SNR', 'R_SNR', 'L_SNRAS'}; 
metric_names = {'Traditional SNR', 'Robust SNR', 'Log-SNRAS'};

B = 2000; % Number of bootstrap iterations
disp(['Starting Stratified Bootstrap with B = ', num2str(B)]);
disp('--------------------------------------------------');

% Determine indices for stratified sampling
pos_idx = find(y_true_clean == 1);
neg_idx = find(y_true_clean == 0);

% Set fixed random seed to guarantee exact reproducibility
rng(42); 

% 3. Calculate AUC and 95% CI for each evaluated metric
for k = 1:length(score_columns)
    
    current_col = score_columns{k};
    scores = data_clean.(current_col);
    scores(isnan(scores)) = 0; % Handle missing values
    
    % Invert scores for AUC calculation logic
    scores_to_test = -scores; 
    
    auc_results = zeros(B, 1);
    
    % Execute bootstrap iterations
    for b = 1:B
        % Stratified sampling with replacement
        sample_pos = datasample(pos_idx, length(pos_idx));
        sample_neg = datasample(neg_idx, length(neg_idx));
        boot_idx = [sample_pos; sample_neg];
        
        % Compute AUC for the current bootstrap sample
        [~, ~, ~, auc_results(b)] = perfcurve(y_true_clean(boot_idx), scores_to_test(boot_idx), 1);
    end
    
    % Compute summary statistics
    final_mean_auc = mean(auc_results);
    ci_95 = prctile(auc_results, [2.5, 97.5]);
    
    % Display metric evaluation results
    fprintf('%s:\n', metric_names{k});
    fprintf('Mean AUC = %.3f | 95%% CI = [%.3f, %.3f]\n\n', final_mean_auc, ci_95(1), ci_95(2));
end

disp('Bootstrap evaluation complete.');
