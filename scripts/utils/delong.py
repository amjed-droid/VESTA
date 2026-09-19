"""
delong.py — Exact DeLong test for comparing two correlated ROC curves.

Implements the U-statistics-based method of:
    DeLong, E. R., DeLong, D. M., & Clarke-Pearson, D. L. (1988).
    Comparing the areas under two or more correlated receiver operating
    characteristic curves: a nonparametric approach. Biometrics, 44(3), 837-845.

Usage
-----
    from scripts.utils.delong import delong_test
    z, p = delong_test(y_true, scores1, scores2)
"""

import numpy as np
from scipy import stats


def _placement_matrix(scores: np.ndarray,
                      idx_pos: np.ndarray,
                      idx_neg: np.ndarray) -> np.ndarray:
    """
    Compute the (n_pos × n_neg) U-statistic placement matrix V where:
        V[i, j] = 1    if scores[pos_i] > scores[neg_j]
                  0.5  if scores[pos_i] == scores[neg_j]
                  0    otherwise
    """
    s_pos = scores[idx_pos][:, np.newaxis]   # (n_pos, 1)
    s_neg = scores[idx_neg][np.newaxis, :]   # (1, n_neg)
    V = (s_pos > s_neg).astype(float) + 0.5 * (s_pos == s_neg).astype(float)
    return V                                 # (n_pos, n_neg)


def delong_test(labels: np.ndarray,
                scores1: np.ndarray,
                scores2: np.ndarray) -> tuple[float, float]:
    """
    Exact DeLong test comparing AUC(scores1) vs AUC(scores2) on the same
    binary classification problem.

    Parameters
    ----------
    labels  : array-like of int, shape (n,)
        Ground-truth binary labels (1 = positive, 0 = negative).
    scores1 : array-like of float, shape (n,)
        Prediction scores for classifier 1.
    scores2 : array-like of float, shape (n,)
        Prediction scores for classifier 2.

    Returns
    -------
    z_score : float
        DeLong z-statistic (positive means scores1 has higher AUC).
    p_value : float
        Two-sided p-value from the standard normal distribution.
    """
    labels  = np.asarray(labels,  dtype=float).ravel()
    scores1 = np.asarray(scores1, dtype=float).ravel()
    scores2 = np.asarray(scores2, dtype=float).ravel()

    if len(np.unique(labels)) < 2:
        raise ValueError("labels must contain at least two distinct classes.")
    if not (len(labels) == len(scores1) == len(scores2)):
        raise ValueError("labels, scores1, and scores2 must have the same length.")

    idx_pos = np.where(labels == 1)[0]
    idx_neg = np.where(labels == 0)[0]
    n_pos, n_neg = len(idx_pos), len(idx_neg)

    if n_pos == 0 or n_neg == 0:
        raise ValueError("Both positive and negative samples must be present.")

    # Compute placement matrices
    V1 = _placement_matrix(scores1, idx_pos, idx_neg)  # (n_pos, n_neg)
    V2 = _placement_matrix(scores2, idx_pos, idx_neg)

    # Point AUCs
    auc1 = V1.mean()
    auc2 = V2.mean()

    # Per-sample averages for covariance estimation
    V1_pos = V1.mean(axis=1)  # (n_pos,) — averaged over neg for each pos
    V1_neg = V1.mean(axis=0)  # (n_neg,) — averaged over pos for each neg
    V2_pos = V2.mean(axis=1)
    V2_neg = V2.mean(axis=0)

    # 2×2 covariance matrices (between the two classifiers)
    # cov_pos[i,j] = Cov(V1_pos[i], V2_pos[j])
    S_pos = np.cov(np.stack([V1_pos, V2_pos])) / n_pos   # (2,2)
    S_neg = np.cov(np.stack([V1_neg, V2_neg])) / n_neg   # (2,2)
    S = S_pos + S_neg                                      # combined

    delta = auc1 - auc2
    var_delta = S[0, 0] + S[1, 1] - 2.0 * S[0, 1]

    if var_delta <= 0:
        if delta == 0.0:
            return 0.0, 1.0
        return (np.inf if delta > 0 else -np.inf), 0.0

    se = np.sqrt(var_delta)
    z_score = delta / se
    p_value = 2.0 * stats.norm.sf(abs(z_score))

    return z_score, p_value


def compute_auc(labels: np.ndarray, scores: np.ndarray) -> float:
    """
    Compute AUC via the Wilcoxon–Mann–Whitney U-statistic.

    Parameters
    ----------
    labels : array-like of int  (1 = positive, 0 = negative)
    scores : array-like of float

    Returns
    -------
    auc : float in [0, 1]
    """
    labels = np.asarray(labels, dtype=float).ravel()
    scores = np.asarray(scores, dtype=float).ravel()
    idx_pos = np.where(labels == 1)[0]
    idx_neg = np.where(labels == 0)[0]
    if len(idx_pos) == 0 or len(idx_neg) == 0:
        return 0.5
    V = _placement_matrix(scores, idx_pos, idx_neg)
    return float(V.mean())
