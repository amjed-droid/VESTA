function [y_true_v, scores_inv, raw_auc, boot_pen, auc_penalty, n_valid] = analyze_performance(results)
%ANALYZE_PERFORMANCE Statistical validation pipeline for Log-SNRAS benchmarking.
%
%   [Y_TRUE_V, SCORES_INV, RAW_AUC, BOOT_PEN, AUC_PENALTY, N_VALID] = ...
%       ANALYZE_PERFORMANCE(RESULTS) performs the complete statistical
%   validation procedure underlying Table 3 and Section 3.4 of the
%   manuscript "Log-SNRAS: A Computationally Efficient Variance-Stabilized
%   Metric for Vetting Heteroscedastic Light Curves."
%
%   The function evaluates the discriminative performance of five
%   candidate vetting metrics (Traditional SNR, Robust SNR, Pont SNR,
%   BLS-based SNR, and the composite Log-SNRAS score) against an
%   independent ephemeris/literature ground truth, and additionally
%   evaluates the isolated logarithmic penalty and its underlying
%   dispersion contrast as scale-invariant discriminators.
%
%   METHODOLOGY
%   -----------
%   1. Sample selection: Only segments with an unambiguous ground-truth
%      label, valid T-SNR/L-SNRAS values, and at least 30 in-transit
%      data points (N_in >= 30) are retained, per the sampling-
%      consistency requirement derived in Appendix B.
%   2. Point estimates: The Area Under the ROC Curve (AUC) is computed
%      for each of the five raw metrics, and for the negated penalty
%      and dispersion contrast (so that higher values indicate a
%      greater likelihood of a genuine transit).
%   3. Stratified bootstrap (B = 2000): Confirmed and Artifact classes
%      are resampled independently with replacement, preserving class
%      proportions, to construct 95% confidence intervals (2.5th and
%      97.5th percentiles) for every metric.
%   4. Jackknife stability: The AUC of the isolated penalty is
%      recomputed after sequentially removing each segment, to assess
%      whether discriminative performance is evenly distributed across
%      the sample or driven by a small subset of points.
%   5. Exact DeLong test: The ROC curve of the isolated penalty is
%      formally compared against each of the five candidate metrics via
%      the analytic DeLong statistic, yielding a z-value and two-sided
%      p-value per comparison.
%   6. Permutation test (10,000 shuffles): Class labels are randomly
%      permuted while penalty scores are held fixed, generating a null
%      distribution of AUC values from which an empirical p-value is
%      derived.
%
%   INPUT
%   -----
%   results : table
%       Unified evaluation table (as produced by the master pipeline),
%       expected to contain at least the following variables: Label,
%       T_SNR, R_SNR, P_SNR, B_SNR, L_SNRAS, Psi, Penalty_pct, N_in.
%
%   OUTPUTS
%   -------
%   y_true_v    : Ground-truth labels (1 = Confirmed, 0 = Artifact) for
%                 the quality-filtered valid subset.
%   scores_inv  : Negated penalty scores (-Penalty_pct) used as the
%                 primary scale-invariant discriminator.
%   raw_auc     : AUC values for the five raw SNR-based metrics, in the
%                 order {T-SNR, R-SNR, P-SNR, B-SNR, Log-SNRAS}.
%   boot_pen    : Bootstrap distribution (B x 1) of the AUC of the
%                 isolated penalty, used for confidence interval
%                 estimation and downstream diagnostic plots.
%   auc_penalty : Point estimate of the AUC of the isolated penalty.
%   n_valid     : Number of segments retained in the quality-filtered
%                 evaluation subset.
%
%   See also: DELONGEXACT, PERFCURVE

    % --- 1. Sample selection: quality-filtered valid subset ---
    valid = ~isnan(results.Label) & ~isnan(results.T_SNR) & ...
            ~isnan(results.L_SNRAS) & results.N_in >= 30;

    fprintf('Confirmed (valid): %d | Artifact (valid): %d\n', ...
        sum(results.Label(valid)==1), sum(results.Label(valid)==0));

    y_true_v = results.Label(valid);
    n_valid  = sum(valid);
    fprintf('\n=== TABLE 3: Performance Benchmarking (N=%d valid) ===\n', n_valid);

    % --- 2. Point-estimate AUC for the five raw SNR-based metrics ---
    raw_methods = {'T-SNR','R-SNR','P-SNR','B-SNR','Log-SNRAS'};
    raw_cols    = {'T_SNR','R_SNR','P_SNR','B_SNR','L_SNRAS'};
    n_raw = length(raw_methods);
    raw_auc = zeros(n_raw,1);
    for m = 1:n_raw
        scores = results.(raw_cols{m})(valid);
        scores(isnan(scores)) = 0;
        [~,~,~,raw_auc(m)] = perfcurve(y_true_v, scores, 1);
    end

    % --- 3. Point-estimate AUC for the scale-invariant discriminators ---
    penalty_scores = results.Penalty_pct(valid);
    psi_scores     = results.Psi(valid);
    [~,~,~,auc_penalty] = perfcurve(y_true_v, -penalty_scores, 1);
    [~,~,~,auc_psi]     = perfcurve(y_true_v, -psi_scores, 1);

    fprintf('\n-- Raw SNR-based AUC --\n');
    for m = 1:n_raw, fprintf('%-12s : AUC = %.3f\n', raw_methods{m}, raw_auc(m)); end
    fprintf('\n-- Scale-invariant --\n');
    fprintf('%-12s : AUC = %.3f\n', 'Penalty (-P)', auc_penalty);
    fprintf('%-12s : AUC = %.3f\n', 'Dispersion(-Psi)', auc_psi);

    % --- 4. Stratified bootstrap (B = 2000) for confidence intervals ---
    cfg_B = 2000;
    fprintf('\nRunning stratified bootstrap (B=%d)...\n', cfg_B);
    pos_idx = find(y_true_v == 1);
    neg_idx = find(y_true_v == 0);
    boot_raw = zeros(cfg_B, n_raw);
    boot_pen = zeros(cfg_B, 1);
    boot_psi = zeros(cfg_B, 1);

    for b = 1:cfg_B
        sample_pos = datasample(pos_idx, length(pos_idx));
        sample_neg = datasample(neg_idx, length(neg_idx));
        boot_idx = [sample_pos(:); sample_neg(:)];
        yb = y_true_v(boot_idx);

        for m = 1:n_raw
            s_all = results.(raw_cols{m})(valid);
            s_all(isnan(s_all)) = 0;
            [~,~,~,boot_raw(b,m)] = perfcurve(yb, s_all(boot_idx), 1);
        end
        [~,~,~,boot_pen(b)] = perfcurve(yb, -penalty_scores(boot_idx), 1);
        [~,~,~,boot_psi(b)] = perfcurve(yb, -psi_scores(boot_idx), 1);
    end

    % --- 5. Bootstrap 95% confidence intervals ---
    ci_raw = cell(n_raw,1);
    fprintf('\n-- Bootstrap 95%% CIs: raw SNR --\n');
    for m = 1:n_raw
        ci = prctile(boot_raw(:,m), [2.5, 97.5]);
        ci_raw{m} = sprintf('[%.3f, %.3f]', ci(1), ci(2));
        fprintf('%-12s : point=%.3f | CI=%s\n', raw_methods{m}, raw_auc(m), ci_raw{m});
    end

    ci_pen = prctile(boot_pen, [2.5, 97.5]);
    ci_psi = prctile(boot_psi, [2.5, 97.5]);
    fprintf('\n-- Bootstrap 95%% CIs: scale-invariant --\n');
    fprintf('Penalty (-P)     : point=%.3f | CI=[%.3f, %.3f]\n', auc_penalty, ci_pen(1), ci_pen(2));
    fprintf('Dispersion (-Psi): point=%.3f | CI=[%.3f, %.3f]\n', auc_psi, ci_psi(1), ci_psi(2));

    % --- 6. Jackknife stability analysis ---
    fprintf('\n=== Jackknife Stability Analysis ===\n');
    n_segments = length(y_true_v);
    auc_jack = zeros(n_segments, 1);
    scores_inv = -penalty_scores;
    for i = 1:n_segments
        idx = true(n_segments, 1);
        idx(i) = false;
        [~,~,~,auc_jack(i)] = perfcurve(y_true_v(idx), scores_inv(idx), 1);
    end
    fprintf('Jackknife AUC: mean = %.4f, std = %.6f\n', mean(auc_jack), std(auc_jack));

    % --- 7. Exact DeLong test: isolated penalty vs. each raw metric ---
    fprintf('\n=== DeLong Exact Test: Penalty vs. Raw SNR ===\n');
    methods = {'T_SNR','R_SNR','P_SNR','B_SNR','L_SNRAS'};
    for m = 1:length(methods)
        scores_raw = results.(methods{m})(valid);
        scores_raw(isnan(scores_raw)) = 0;
        [p_val_delong, z_val] = deLongExact(y_true_v, scores_inv, scores_raw);
        [~,~,~,auc_raw] = perfcurve(y_true_v, scores_raw, 1);
        fprintf('%s vs Penalty: diff AUC=%.3f, z=%.3f, p=%.4f\n', ...
            methods{m}, auc_penalty - auc_raw, z_val, p_val_delong);
    end

    % --- 8. Stratified permutation test (10,000 shuffles) ---
    fprintf('\n=== Permutation Test ===\n');
    n_perm = 10000;
    auc_perm = zeros(n_perm, 1);
    for i = 1:n_perm
        y_perm = y_true_v(randperm(length(y_true_v)));
        [~,~,~,auc_perm(i)] = perfcurve(y_perm, scores_inv, 1);
    end
    p_val_perm = mean(auc_perm >= auc_penalty);
    fprintf('Permutation p-value = %.4f\n', p_val_perm);

end
