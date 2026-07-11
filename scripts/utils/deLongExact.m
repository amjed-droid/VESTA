function [p_value, z_score] = deLongExact(labels, scores1, scores2)
% DELONGEXACT Computes the exact DeLong test for comparing two correlated ROC curves.
%
%   [P, Z] = DELONGEXACT(LABELS, SCORES1, SCORES2) performs the DeLong 
%   exact test to compare the Area Under the ROC Curve (AUC) of two 
%   classifiers evaluated on the same sample.
%
%   METHODOLOGY (U-statistics)
%   --------------------------
%   The DeLong test is based on the covariance matrix of the AUC estimates
%   using the theory of U-statistics. It does not rely on assumptions of
%   normality of the scores, making it highly robust for evaluating
%   correlated classification metrics. The two-sided p-value is derived
%   from the standard normal distribution of the z-statistic:
%
%       z = (AUC1 - AUC2) / SE_delta
%
%   INPUTS
%   ------
%   labels   : (N x 1) Vector of ground-truth binary labels (1 = positive, 0 = negative).
%   scores1  : (N x 1) Vector of scores for the first classifier.
%   scores2  : (N x 1) Vector of scores for the second classifier.
%
%   OUTPUTS
%   -------
%   p_value  : Two-sided p-value of the DeLong test.
%   z_score  : Z-statistic of the DeLong test.
%
%   SEE ALSO
%   --------
%   perfcurve, generate_table4


if nargin < 3
    error('Inputs labels, scores1, and scores2 are all required.');
end
if length(unique(labels)) < 2
    error('Labels must contain at least two distinct classes (positive and negative).');
end

% Ensure all inputs are column vectors
labels = labels(:);
scores1 = scores1(:);
scores2 = scores2(:);

% Identify positive and negative indices
pos = 1; neg = 0;
idx_pos = find(labels == pos);
idx_neg = find(labels == neg);
n_pos = length(idx_pos);
n_neg = length(idx_neg);

% Helper to compute the placement matrix (U-statistic kernel)
function V = computeV(scores)
    s_pos = scores(idx_pos);
    s_neg = scores(idx_neg);
    V = (s_pos > s_neg') + 0.5 * (s_pos == s_neg');
end

% Compute placement matrices
V1 = computeV(scores1);
V2 = computeV(scores2);

% Compute AUCs
auc1 = mean(V1(:));
auc2 = mean(V2(:));

% Compute covariance components
% Average across negative samples for each positive sample
V1_pos_avg = mean(V1, 2); 
V2_pos_avg = mean(V2, 2);

% Average across positive samples for each negative sample
V1_neg_avg = mean(V1, 1)'; 
V2_neg_avg = mean(V2, 1)';

% Covariances (2x2 matrices)
cov_pos = cov(V1_pos_avg, V2_pos_avg);
cov_neg = cov(V1_neg_avg, V2_neg_avg);

% Combined covariance matrix
S = (cov_pos ./ n_pos) + (cov_neg ./ n_neg);

% Test statistic
delta = auc1 - auc2;
se = sqrt(S(1,1) + S(2,2) - 2*S(1,2));

if se == 0
    if delta == 0
        z_score = 0;
        p_value = 1.0;
    else
        z_score = Inf;
        p_value = 0.0;
    end
else
    z_score = delta / se;
    p_value = 2 * normcdf(-abs(z_score));
end
end
