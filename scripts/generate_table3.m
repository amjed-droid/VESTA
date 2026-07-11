%% GENERATE_TABLE3 Computes AUC and Bootstrap performance results.
%
%   GENERATE_TABLE3() reads evaluation_dataset_v2.csv, applies quality filtering,
%   computes point AUCs and stratified bootstrap 95% CIs, and prints 
%   the results as plain text to the Command Window.
%
%   SEE ALSO
%   --------
%   analyze_performance, perfcurve, datasample
function generate_table3()
    % --- 1. Load and Filter Data ---
    csv_file = 'evaluation_dataset_v2.csv';
    if ~isfile(csv_file)
        error('File %s not found. Ensure it exists.', csv_file);
    end
    results = readtable(csv_file);

    % Apply quality filtering
    valid_mask = ~isnan(results.Label) & ~isnan(results.T_SNR) & ...
                 ~isnan(results.L_SNRAS) & (results.N_in >= 30);
    results_valid = results(valid_mask, :);
    y_true_v = results_valid.Label;
    n_valid = sum(valid_mask);

    % --- 2. Setup Metrics ---
    raw_methods = {'T-SNR', 'R-SNR', 'P-SNR', 'B-SNR', 'Log-SNRAS'};
    raw_cols    = {'T_SNR', 'R_SNR', 'P_SNR', 'B_SNR', 'L_SNRAS'};
    n_raw = length(raw_methods);

    % Compute point AUCs for Raw metrics
    raw_auc = zeros(n_raw,1);
    raw_scores_cell = cell(n_raw, 1);
    for m = 1:n_raw
        scores = results_valid.(raw_cols{m});
        scores(isnan(scores)) = 0;
        raw_scores_cell{m} = scores;
        [~,~,~,raw_auc(m)] = perfcurve(y_true_v, scores, 1);
    end

    % Compute point AUCs for Scale-invariant metrics
    penalty_scores = results_valid.Penalty_pct;
    psi_scores     = results_valid.Psi;
    [~,~,~,auc_penalty] = perfcurve(y_true_v, -penalty_scores, 1);
    [~,~,~,auc_psi]     = perfcurve(y_true_v, -psi_scores, 1);

    % --- 3. Stratified Bootstrap (B = 2000) ---
    B_boot = 2000;
    fprintf('Running stratified bootstrap (B = %d)...\n', B_boot);
    pos_idx = find(y_true_v == 1);
    neg_idx = find(y_true_v == 0);

    boot_raw = zeros(B_boot, n_raw);
    boot_pen = zeros(B_boot, 1);
    boot_psi = zeros(B_boot, 1);

    for b = 1:B_boot
        sample_pos = datasample(pos_idx, length(pos_idx));
        sample_neg = datasample(neg_idx, length(neg_idx));
        boot_idx = [sample_pos(:); sample_neg(:)];
        yb = y_true_v(boot_idx);

        for m = 1:n_raw
            [~,~,~,boot_raw(b,m)] = perfcurve(yb, raw_scores_cell{m}(boot_idx), 1);
        end
        [~,~,~,boot_pen(b)] = perfcurve(yb, -penalty_scores(boot_idx), 1);
        [~,~,~,boot_psi(b)] = perfcurve(yb, -psi_scores(boot_idx), 1);
    end

    % --- 4. Compute Confidence Intervals ---
    ci_raw = zeros(n_raw, 2);
    for m = 1:n_raw
        ci_raw(m,:) = prctile(boot_raw(:,m), [2.5, 97.5]);
    end
    ci_pen = prctile(boot_pen, [2.5, 97.5]);
    ci_psi = prctile(boot_psi, [2.5, 97.5]);

    % --- 5. Print Results as Plain Text ---
    fprintf('\nTable 3: Performance Benchmarking (N = %d)\n', n_valid);
    fprintf('Raw SNR-based (reference only):\n');
    for m = 1:n_raw
        fprintf('%s: AUC = %.3f, CI = [%.3f, %.3f]\n', ...
            raw_methods{m}, raw_auc(m), ci_raw(m,1), ci_raw(m,2));
    end
    fprintf('\nScale-invariant (primary result):\n');
    fprintf('Penalty (-P): AUC = %.3f, CI = [%.3f, %.3f]\n', auc_penalty, ci_pen(1), ci_pen(2));
    fprintf('Dispersion (-Psi): AUC = %.3f, CI = [%.3f, %.3f]\n', auc_psi, ci_psi(1), ci_psi(2));
end
