%% REPRODUCE_TABLE3_ROC Reproduces Table 1 & Figure 3 (ROC) from the manuscript
% =========================================================================
% Log-SNRAS: A Computationally Efficient Variance-Stabilized Metric 
% for Vetting Heteroscedastic Light Curves
%
% This script loads the curated multi-host benchmark evaluation results,
% computes the AUC for all metrics (including shape-corrected Log-SNRAS),
% computes stratified bootstrap 95% confidence intervals (B = 2000),
% and generates the publication ROC curve.
% =========================================================================

function reproduce_table3_roc()
    fprintf('====================================================================\n');
    fprintf(' REPRODUCING BENCHMARK CLASSIFICATION PERFORMANCE (TABLE 1 & ROC)\n');
    fprintf('====================================================================\n\n');

    % --- 1. Locate and Load Evaluation Dataset ---
    csv_file = fullfile('data', 'multi_host_evaluation_results.csv');
    if ~isfile(csv_file)
        csv_file = 'multi_host_evaluation_results.csv';
    end
    if ~isfile(csv_file)
        error('Evaluation results not found. Run "python scripts/evaluate_multi_host_benchmark.py" first.');
    end

    data = readtable(csv_file);
    fprintf('Loaded %d targets from: %s\n\n', height(data), csv_file);

    % Ground truth labels: 1 = Confirmed Planet, 0 = Non-Planetary Artifact
    y_true = data.label;
    
    % --- 2. Define Metrics to Evaluate ---
    % Note: Lower penalty = more likely planetary transit -> score = -penalty
    methods = {
        'Shape-Corrected Penalty', ...
        'Uncorrected Penalty', ...
        'Traditional SNR (T-SNR)', ...
        'Robust SNR (MAD)', ...
        'Pont SNR (2006)', ...
        'BLS SNR Proxy', ...
        'Log-SNRAS Composite', ...
        'Inverse Depth (1/delta)'
    };

    scores_matrix = [
        -data.penalty_corr, ...
        -data.penalty_raw, ...
        data.t_snr, ...
        data.r_snr, ...
        data.p_snr, ...
        data.b_snr, ...
        data.l_snras_corr, ...
        -data.depth_measured_ppm
    ];

    n_methods = length(methods);
    point_aucs = zeros(n_methods, 1);
    ci_lower = zeros(n_methods, 1);
    ci_upper = zeros(n_methods, 1);

    % --- 3. Compute AUC and Stratified Bootstrap 95% CIs (B = 2000) ---
    rng(42); % Fixed random seed for exact reproducibility
    B = 2000;
    n_samples = length(y_true);

    for m = 1:n_methods
        s = scores_matrix(:, m);
        
        % Point AUC using Wilcoxon rank-sum / trapezoidal integration
        point_aucs(m) = compute_simple_auc(y_true, s);

        % Stratified Bootstrap
        boot_aucs = zeros(B, 1);
        pos_idx = find(y_true == 1);
        neg_idx = find(y_true == 0);

        for b = 1:B
            b_pos = pos_idx(randi(length(pos_idx), length(pos_idx), 1));
            b_neg = neg_idx(randi(length(neg_idx), length(neg_idx), 1));
            b_idx = [b_pos; b_neg];
            boot_aucs(b) = compute_simple_auc(y_true(b_idx), s(b_idx));
        end

        ci_lower(m) = prctile(boot_aucs, 2.5);
        ci_upper(m) = prctile(boot_aucs, 97.5);
    end

    % --- 4. Display Results Table ---
    fprintf('%-30s | %-8s | %-18s\n', 'Method', 'AUC', '95% Bootstrap CI');
    fprintf('--------------------------------------------------------------------\n');
    for m = 1:n_methods
        fprintf('%-30s | %6.3f   | [%5.3f, %5.3f]\n', ...
            methods{m}, point_aucs(m), ci_lower(m), ci_upper(m));
    end
    fprintf('====================================================================\n\n');

    % --- 5. Generate and Save ROC Plot ---
    fig = figure('Visible', 'off', 'Position', [100, 100, 700, 550]);
    hold on; box on; grid on;

    colors = lines(n_methods);
    for m = [1, 2, 3, 7, 8]
        [fpr, tpr] = compute_roc_curve(y_true, scores_matrix(:, m));
        if m == 1
            plot(fpr, tpr, 'Color', [0, 0.2, 0.6], 'LineWidth', 2.5, ...
                'DisplayName', sprintf('%s (AUC = %.3f)', methods{m}, point_aucs(m)));
        else
            plot(fpr, tpr, 'Color', colors(m, :), 'LineWidth', 1.5, ...
                'DisplayName', sprintf('%s (AUC = %.3f)', methods{m}, point_aucs(m)));
        end
    end
    plot([0, 1], [0, 1], 'k--', 'LineWidth', 1.2, 'DisplayName', 'Random Chance (AUC = 0.500)');

    xlabel('False Positive Rate (1 - Specificity)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('True Positive Rate (Sensitivity)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Multi-Host Benchmark ROC Analysis (N = 17)', 'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'southeast', 'FontSize', 9);
    xlim([0, 1]); ylim([0, 1.02]);

    out_fig = fullfile('figures', 'Figure3_ROC_reproduced.png');
    saveas(fig, out_fig);
    close(fig);
    fprintf('Saved reproduced ROC curve to: %s\n', out_fig);
end

% Helper function: AUC calculation
function auc = compute_simple_auc(labels, scores)
    pos_scores = scores(labels == 1);
    neg_scores = scores(labels == 0);
    n_pos = length(pos_scores);
    n_neg = length(neg_scores);
    if n_pos == 0 || n_neg == 0
        auc = 0.5;
        return;
    end
    count = 0;
    for i = 1:n_pos
        count = count + sum(pos_scores(i) > neg_scores) + 0.5 * sum(pos_scores(i) == neg_scores);
    end
    auc = count / (n_pos * n_neg);
end

% Helper function: ROC curve calculation
function [fpr, tpr] = compute_roc_curve(labels, scores)
    thresholds = sort(unique(scores), 'descend');
    thresholds = [-inf; thresholds; inf];
    n_pos = sum(labels == 1);
    n_neg = sum(labels == 0);
    tpr = zeros(length(thresholds), 1);
    fpr = zeros(length(thresholds), 1);
    for t = 1:length(thresholds)
        pred = scores >= thresholds(t);
        tpr(t) = sum(pred & (labels == 1)) / n_pos;
        fpr(t) = sum(pred & (labels == 0)) / n_neg;
    end
end
