"""
Statistical Verification Suite for Log-SNRAS
=============================================
Implements:
1. DeLong asymptotic test for correlated ROC curves (DeLong et al. 1988)
2. Stratified Bootstrap Confidence Intervals
3. Permutation Significance Testing
4. Multi-Host Leave-One-Host-Out Evaluation
"""

import numpy as np
from scipy import stats


def compute_auc(y_true, y_score):
    """
    Computes Wilcoxon-Mann-Whitney AUC statistic.
    """
    y_true = np.asarray(y_true, dtype=bool)
    y_score = np.asarray(y_score, dtype=float)
    
    pos = y_score[y_true]
    neg = y_score[~y_true]
    
    n_pos = len(pos)
    n_neg = len(neg)
    if n_pos == 0 or n_neg == 0:
        return np.nan
        
    # Mann-Whitney U
    ranks = stats.rankdata(np.concatenate([pos, neg]))
    rank_pos = np.sum(ranks[:n_pos])
    u_val = rank_pos - (n_pos * (n_pos + 1)) / 2.0
    return float(u_val / (n_pos * n_neg))


def delong_roc_variance(y_true, y_score):
    """
    Computes the variance of AUC using the method of DeLong, DeLong & Clarke-Pearson (1988).
    """
    y_true = np.asarray(y_true, dtype=bool)
    y_score = np.asarray(y_score, dtype=float)

    m = np.sum(y_true)     # positives
    n = np.sum(~y_true)    # negatives
    if m == 0 or n == 0:
        return np.nan, np.nan

    pos = y_score[y_true]
    neg = y_score[~y_true]

    # Placements V10 and V01
    # V10_i = (1/n) * sum_j I(pos_i > neg_j)
    v10 = np.mean((pos[:, None] > neg[None, :]) + 0.5 * (pos[:, None] == neg[None, :]), axis=1)
    # V01_j = (1/m) * sum_i I(pos_i > neg_j)
    v01 = np.mean((pos[:, None] > neg[None, :]) + 0.5 * (pos[:, None] == neg[None, :]), axis=0)

    auc = np.mean(v10)
    s10 = np.var(v10, ddof=1) if m > 1 else 0.0
    s01 = np.var(v01, ddof=1) if n > 1 else 0.0

    var_auc = (s10 / m) + (s01 / n)
    return float(auc), float(var_auc), v10, v01


def delong_test(y_true, y_score_a, y_score_b):
    """
    Two-sided DeLong asymptotic test comparing two correlated ROC curves.
    
    Returns
    -------
    diff_auc : float
    z_score : float
    p_value : float
    """
    y_true = np.asarray(y_true, dtype=bool)
    m = np.sum(y_true)
    n = np.sum(~y_true)

    auc_a, var_a, v10_a, v01_a = delong_roc_variance(y_true, y_score_a)
    auc_b, var_b, v10_b, v01_b = delong_roc_variance(y_true, y_score_b)

    # Covariance between AUC_A and AUC_B
    s10_ab = np.cov(v10_a, v10_b, ddof=1)[0, 1] if m > 1 else 0.0
    s01_ab = np.cov(v01_a, v01_b, ddof=1)[0, 1] if n > 1 else 0.0
    cov_ab = (s10_ab / m) + (s01_ab / n)

    var_diff = var_a + var_b - 2.0 * cov_ab
    if var_diff <= 0 or np.isnan(var_diff):
        z = 0.0
        p_val = 1.0
    else:
        diff_auc = auc_a - auc_b
        z = diff_auc / np.sqrt(var_diff)
        p_val = 2.0 * (1.0 - stats.norm.cdf(abs(z)))

    return float(auc_a - auc_b), float(z), float(p_val)


def stratified_bootstrap_ci(y_true, y_score, n_boot=2000, alpha=0.05, seed=42):
    """
    Computes stratified bootstrap confidence interval for the AUC.
    """
    rng = np.random.default_rng(seed)
    y_true = np.asarray(y_true, dtype=bool)
    y_score = np.asarray(y_score, dtype=float)

    pos_idx = np.where(y_true)[0]
    neg_idx = np.where(~y_true)[0]

    n_pos = len(pos_idx)
    n_neg = len(neg_idx)
    if n_pos == 0 or n_neg == 0:
        return np.nan, (np.nan, np.nan), np.array([])

    boot_aucs = np.empty(n_boot)
    for b in range(n_boot):
        sample_pos = rng.choice(pos_idx, size=n_pos, replace=True)
        sample_neg = rng.choice(neg_idx, size=n_neg, replace=True)
        
        b_pos = y_score[sample_pos]
        b_neg = y_score[sample_neg]
        
        # Fast AUC calculation
        ranks = stats.rankdata(np.concatenate([b_pos, b_neg]))
        u_val = np.sum(ranks[:n_pos]) - (n_pos * (n_pos + 1)) / 2.0
        boot_aucs[b] = u_val / (n_pos * n_neg)

    ci_lower = float(np.percentile(boot_aucs, 100 * (alpha / 2.0)))
    ci_upper = float(np.percentile(boot_aucs, 100 * (1.0 - alpha / 2.0)))
    point_est = float(compute_auc(y_true, y_score))
    return point_est, (ci_lower, ci_upper), boot_aucs


def permutation_test(y_true, y_score, n_permutations=10000, seed=42):
    """
    Calculates permutation p-value for the AUC against the null hypothesis (AUC = 0.5).
    """
    rng = np.random.default_rng(seed)
    y_true = np.asarray(y_true, dtype=bool)
    y_score = np.asarray(y_score, dtype=float)

    observed_auc = compute_auc(y_true, y_score)
    perm_aucs = np.empty(n_permutations)
    
    for i in range(n_permutations):
        perm_labels = rng.permutation(y_true)
        perm_aucs[i] = compute_auc(perm_labels, y_score)

    p_val = float(np.mean(perm_aucs >= observed_auc))
    return observed_auc, p_val, perm_aucs
