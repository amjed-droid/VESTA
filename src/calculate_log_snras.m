function [log_snras, penalty_term, psi] = calculate_log_snras(flux, snr_trad, transit_mask)
% CALCULATE_LOG_SNRAS computes the Variance-Stabilized Vetting Metric.
%
%   [log_snras, penalty_term, psi] = calculate_log_snras(flux, snr_trad, transit_mask)
%
%   This function implements the Log-SNRAS metric as described in:
%   "Log-SNRAS: A Computationally Efficient Variance-Stabilized Metric for 
%    Vetting Heteroscedastic Light Curves" (Submitted to Astronomy and Computing).
%
%   INPUTS:
%       flux         : (Vector) Normalized light curve flux array.
%       snr_trad     : (Scalar) The traditional Signal-to-Noise Ratio (Depth/Sigma_out * sqrt(N)).
%       transit_mask : (Boolean Vector) Logical mask where True indicates in-transit points.
%
%   OUTPUTS:
%       log_snras    : (Scalar) The final penalized detection metric.
%       penalty_term : (Scalar) The denominator scaling factor (1 + ln(1+psi)).
%       psi          : (Scalar) The normalized dispersion contrast.
%
%   COMPLEXITY:
%       Time Complexity: O(N) - Linear relative to light curve length.
%
%   AUTHOR:
%       Ahmed Sattar Jabbar
%       Department of Statistics, Mustansiriyah University.

    %% 1. Input Validation
    if nargin < 3
        error('Error: Transit mask is required to separate In-Transit from Out-of-Transit regions.');
    end
    
    %% 2. Partitioning (Data Splitting)
    f_in = flux(transit_mask);   
    f_out = flux(~transit_mask); 
    
    % Edge case handling: Ensure we have data in both segments
    if isempty(f_in) || isempty(f_out)
        warning('Empty segment detected. Returning 0 score.');
        log_snras = 0; penalty_term = 1; psi = 0;
        return;
    end
    
    %% 3. Robust Statistics Estimation
    % Using 'omitnan' ensures robustness against missing data points
    sigma_in = std(f_in, 'omitnan');
    sigma_out = std(f_out, 'omitnan');
    
    %% 4. Variance Contrast (Psi Calculation)
    % Calculates the normalized dispersion contrast
    if sigma_out == 0 || isnan(sigma_out)
        psi = 0;
    else
        psi = abs(sigma_in - sigma_out) / sigma_out;
    end
    
    %% 5. Logarithmic Penalty (Entropy-Based)
    % Based on the Lomax Prior assumption for heteroscedastic noise
    penalty_log_term = log(1 + psi);
    
    % The denominator scaling factor
    penalty_term = 1 + penalty_log_term;
    
    %% 6. Final Log-SNRAS Calculation
    % Apply the penalty to the traditional SNR
    log_snras = snr_trad / penalty_term;
    
end