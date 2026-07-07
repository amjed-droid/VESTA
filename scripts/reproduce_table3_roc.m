% =========================================================================
% Script: Stratified Bootstrap AUC and Confidence Intervals Validation
% Author: Ahmed S. Jabbar
% Description: Performs B=2000 bootstrap iterations to calculate the mean
% Area Under the Curve (AUC) and 95% Confidence Intervals for Log-SNRAS
% and baseline metrics.
% =========================================================================
clear; clc; close all;

% --- 1. Load the pre-processed benchmark dataset ---
csv_file = 'Supplementary_Table_C4.csv';
if ~exist(csv_file, 'file')
    error('Could not find %s. Please run generate_table_c4.m first to produce the CSV file.', csv_file);
end

data_csv = readtable(csv_file);
fprintf('Loaded %s with %d targets.\n', csv_file, height(data_csv));

% --- 2. Determine Classification Ground Truth (Confirmed vs Vetoed) ---
ranks = data_csv.Rank;
y_true = nan(size(ranks));
y_true(ranks == 1) = 1; % Confirmed (Rank 1)
y_true(ranks == 3) = 0; % Vetoed (Rank 3)
valid = ~isnan(y_true);
y_true_v = y_true(valid); % Column vector of size 68 x 1

m_names = {'T-SNR', 'R-SNR', 'P-SNR', 'B-SNR', 'Log-SNRAS'};
csv_cols = {'T_SNR', 'R_SNR', 'P_SNR', 'B_SNR', 'L_SNRAS'};
perf_results = zeros(5, 5); 
threshold = 7.1;

fprintf('\nCalculating baseline classification metrics (Youden-optimal threshold = %.1f)...\n', threshold);
for m = 1:5
    current_col = csv_cols{m};
    scores = data_csv.(current_col)(valid);
    y_pred = (scores >= threshold);
    
    tp = sum(y_pred == 1 & y_true_v == 1);
    fp = sum(y_pred == 1 & y_true_v == 0);
    fn = sum(y_pred == 0 & y_true_v == 1);
    tn = sum(y_pred == 0 & y_true_v == 0);
    
    acc = (tp + tn) / length(y_true_v);
    prec = tp / (tp + fp + eps);
    rec = tp / (tp + fn + eps);
    f1 = 2 * (prec * rec) / (prec + rec + eps);
    [~,~,~,auc] = perfcurve(y_true_v, scores, 1);
    perf_results(m, :) = [acc, prec, rec, f1, auc];
end

PerformanceTable = table(m_names', perf_results(:,1), perf_results(:,2), ...
    perf_results(:,3), perf_results(:,4), perf_results(:,5), ...
    'VariableNames', {'Method', 'Accuracy', 'Precision', 'Recall', 'F1_Score', 'AUC'});
disp('--- Baseline Performance Metrics ---');
disp(PerformanceTable);

% --- 3. Stratified Bootstrap (B=2000) for Table 3 CI and Delta AUC ---
fprintf('Running Stratified Bootstrap (B = 2000) to compute 95%% CIs...\n');
B = 2000;
rng(42); % Seed for reproducibility
pos_idx = find(y_true_v == 1);
neg_idx = find(y_true_v == 0);

boot_auc = zeros(B, 5);
for b = 1:B
    sample_pos = datasample(pos_idx, length(pos_idx));
    sample_neg = datasample(neg_idx, length(neg_idx));
    boot_idx = [sample_pos; sample_neg];
    
    for m = 1:5
        current_col = csv_cols{m};
        scores_all = data_csv.(current_col)(valid);
        scores_boot = scores_all(boot_idx);
        
        [~, ~, ~, boot_auc(b, m)] = perfcurve(y_true_v(boot_idx), scores_boot, 1);
    end
end

% Calculate Confidence Intervals and print diagnostic checks
ci_auc = cell(5, 1);
disp('--- Bootstrap Diagnostic Checks ---');
for m = 1:5
    point_est = perf_results(m, 5);
    boot_vals = boot_auc(:, m);
    boot_mean = mean(boot_vals);
    boot_median = median(boot_vals);
    ci_95 = prctile(boot_vals, [2.5, 97.5]);
    ci_auc{m} = sprintf('[%.3f, %.3f]', ci_95(1), ci_95(2));
    
    pct_above = 100 * sum(boot_vals > point_est) / B;
    fprintf('%s: Point estimate=%.3f | Bootstrap mean=%.3f | Bootstrap median=%.3f | CI=[%.3f, %.3f] | %.1f%% of iterations above point estimate\n', ...
        m_names{m}, point_est, boot_mean, boot_median, ci_95(1), ci_95(2), pct_above);
end
disp('----------------------------------');

% Calculate Delta AUC vs Ours (Log-SNRAS is index 5)
delta_auc_str = cell(5, 1);
our_boot_auc = boot_auc(:, 5);
for m = 1:5
    if m == 5
        delta_auc_str{m} = 'Reference';
    else
        delta_vals = our_boot_auc - boot_auc(:, m);
        mean_delta = mean(delta_vals);
        delta_auc_str{m} = sprintf('+%.3f', mean_delta);
    end
end

% --- 4. Print Final Performance Summary Table ---
disp('================================================================================');
disp('                               RESULTS ');
disp('================================================================================');
fprintf('%-10s | %-8s | %-9s | %-6s | %-8s | %-5s | %-15s | %-12s\n', ...
    'Method', 'Accuracy', 'Precision', 'Recall', 'F1-Score', 'AUC', '95% CI (AUC)', 'Delta AUC');
disp('--------------------------------------------------------------------------------');
for m = 1:5
    method_name = m_names{m};
    acc = perf_results(m, 1);
    prec = perf_results(m, 2);
    rec = perf_results(m, 3);
    f1 = perf_results(m, 4);
    auc_val = perf_results(m, 5);
    fprintf('%-10s | %-8.3f | %-9.3f | %-6.3f | %-8.3f | %-5.3f | %-15s | %-12s\n', ...
        method_name, acc, prec, rec, f1, auc_val, ci_auc{m}, delta_auc_str{m});
end
disp('================================================================================');


% --- 5. Generate and Save Professional ROC Plot ---
figure('Color', 'w', 'Position', [100, 100, 800, 600]);
hold on;

colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880]};
lineStyles = {'--', '--', ':', '-.', '-'};
m_labels = {'Traditional SNR (T-SNR)', 'Robust SNR (R-SNR)', 'Pont SNR (P-SNR)', 'BLS Power (B-SNR)', 'Log-SNRAS (Proposed)'};

for m = 1:5
    current_col = csv_cols{m};
    scores_v = data_csv.(current_col)(valid);
    
    [X, Y, ~, auc_val] = perfcurve(y_true_v, scores_v, 1);
    
    if m == 5
        plot(X, Y, 'Color', colors{m}, 'LineStyle', lineStyles{m}, 'LineWidth', 3, ...
            'DisplayName', sprintf('\\textbf{%s (AUC = %.3f)}', m_labels{m}, auc_val));
    else
        plot(X, Y, 'Color', colors{m}, 'LineStyle', lineStyles{m}, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('%s (AUC = %.3f)', m_labels{m}, auc_val));
    end
end

plot([0 1], [0 1], 'k--', 'LineWidth', 1, 'DisplayName', 'Random Classifier (AUC = 0.500)');

xlabel('False Positive Rate (1 - Specificity)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('True Positive Rate (Sensitivity)', 'FontSize', 12, 'FontWeight', 'bold');
title('ROC Curve Analysis: Statistical Benchmarking (N=151)', 'FontSize', 14);
legend('Location', 'southeast', 'FontSize', 10, 'Interpreter', 'latex');
grid on;
axis square;
hold off;

savefig('ROC_Analysis_Final.fig');
saveas(gcf, 'ROC_Analysis_Final.png');
fprintf('Successfully saved ROC Plot to ROC_Analysis_Final.png\n');
