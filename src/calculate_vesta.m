function results = calculate_vesta(flux, transit_mask, varargin)
% CALCULATE_VESTA - Variance-Evaluated Stability for Transit Analysis
%
% Syntax:
%   results = calculate_vesta(flux, transit_mask)
%   results = calculate_vesta(flux, transit_mask, 'Time', time, 'Period', P, 'T0', T0)
%
% Inputs:
%   flux          - 1D array of photometric flux values
%   transit_mask  - Logical array (true = in-transit, false = out-of-transit)
%
% Optional Parameters:
%   'Time'        - 1D array of time timestamps (same length as flux)
%   'Period'      - Orbital period (days)
%   'T0'          - Transit epoch midpoint (days)
%   'NBins'       - Number of phase bins for shape detrending (default: 25)
%
% Outputs:
%   results       - Struct containing:
%       .vesta        - Variance-stabilized detection score
%       .snr_trad     - Classical matched-filter SNR
%       .penalty      - Logarithmic stability penalty P = ln(1 + psi)
%       .psi          - Relative dispersion contrast
%       .sigma_in     - Residual in-transit standard deviation
%       .sigma_out    - Out-of-transit baseline standard deviation
%       .depth        - Measured transit depth
%       .tier         - Statistical vetting tier ('Tier 1', 'Tier 2', 'Tier 3')
%
% Reference:
%   Jabbar, A. S. (2026), "Population-Scale Heteroscedastic Transit Vetting:
%   A Depth-Invariant Variance-Stabilized Framework (VESTA) Benchmarked Across
%   Four Decades in Kepler Light Curves", Astronomy and Computing.

    p = inputParser;
    addRequired(p, 'flux', @isnumeric);
    addRequired(p, 'transit_mask', @islogical);
    addParameter(p, 'Time', [], @isnumeric);
    addParameter(p, 'Period', [], @isnumeric);
    addParameter(p, 'T0', [], @isnumeric);
    addParameter(p, 'NBins', 25, @isnumeric);
    parse(p, flux, transit_mask, varargin{:});

    time = p.Results.Time;
    period = p.Results.Period;
    t0 = p.Results.T0;
    n_bins = p.Results.NBins;

    % Validate masks
    f_in = flux(transit_mask);
    f_out = flux(~transit_mask);

    n_in = length(f_in);
    n_out = length(f_out);

    if n_in < 4 || n_out < 10
        error('VESTA:InsufficientData', 'Insufficient transit or baseline cadences.');
    end

    % Out-of-transit dispersion and baseline median
    sigma_out = std(f_out, 'omitnan');
    baseline_level = median(f_out, 'omitnan');
    in_level = median(f_in, 'omitnan');
    depth = max(0.0, baseline_level - in_level);

    % Non-parametric phase shape detrending if time & ephemeris are provided
    if ~isempty(time) && ~isempty(period) && ~isempty(t0) && period > 0
        t_in = time(transit_mask);
        phases = mod((t_in - t0) / period, 1.0);
        phases(phases > 0.5) = phases(phases > 0.5) - 1.0;

        bin_edges = linspace(min(phases), max(phases), n_bins + 1);
        bin_centers = 0.5 * (bin_edges(1:end-1) + bin_edges(2:end));
        bin_medians = zeros(1, n_bins);

        for b = 1:n_bins
            in_b = phases >= bin_edges(b) & phases < bin_edges(b+1);
            if any(in_b)
                bin_medians(b) = median(f_in(in_b), 'omitnan');
            else
                bin_medians(b) = in_level;
            end
        end

        % Linear interpolation of transit profile
        f_model = interp1(bin_centers, bin_medians, phases, 'linear', in_level);
        residuals = f_in - f_model;
        sigma_in = std(residuals, 'omitnan');
    else
        sigma_in = std(f_in, 'omitnan');
    end

    % Normalized dispersion contrast and logarithmic penalty
    psi = abs(sigma_in - sigma_out) / max(sigma_out, eps);
    penalty = log(1.0 + psi);

    % Traditional matched-filter SNR
    snr_trad = (depth / max(sigma_out, eps)) * sqrt(n_in);

    % Composite VESTA score
    vesta_score = snr_trad / (1.0 + penalty);

    % Calibrated F-distribution operational tier
    if penalty <= 0.15
        tier = 'Tier 1 (Clean Candidate)';
    elseif penalty <= 0.60
        tier = 'Tier 2 (Review Required)';
    else
        tier = 'Tier 3 (Automated Veto)';
    end

    results.vesta = vesta_score;
    results.snr_trad = snr_trad;
    results.penalty = penalty;
    results.psi = psi;
    results.sigma_in = sigma_in;
    results.sigma_out = sigma_out;
    results.depth = depth;
    results.n_in = n_in;
    results.n_out = n_out;
    results.tier = tier;
end
