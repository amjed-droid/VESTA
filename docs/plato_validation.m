%% run_plato_advanced.m
% Advanced batch analysis of simulated PLATO light curves using Log-SNRAS
% Includes injection tests: V-shaped eclipses and intra-transit noise
% Now also computes P-SNR (Pont+ 2006) and T-SNR alongside Traditional SNR
% Measures execution time of Log-SNRAS core function
% Author: Ahmed Sattar Jabbar
% Date: May 2026

%% Setup and global parameters
clear; clc;

% --- Fixed transit parameters (from PLATO Zenodo dataset description https://doi.org/10.5281/zenodo.13939257) ---
T0_DAYS = 0.5;                % First transit centred at 0.5 days
PERIOD_DAYS = 1.0;            % Transit every 24 hours exactly
DURATION_HRS = 4.0;           % Typical transit duration in hours
DURATION_DAYS = DURATION_HRS / 24;

% --- Injection parameters ---
NOISE_LEVEL = 0.0002;         % Intra-transit noise level (Gaussian std)
ECLIPSE_DEPTH = 0.001;        % Depth for V-shaped eclipse injection

%% Locate all CSV files in the current folder
file_list = dir('*.csv');
num_files = length(file_list);

if num_files == 0
    error('No CSV files found in the current directory.');
end

%% Prepare tables to store results for all scenarios
% Original data results (one row per file)
res_orig = table();
res_orig.Filename = strings(num_files, 1);
res_orig.Trad_SNR = zeros(num_files, 1);
res_orig.T_SNR   = zeros(num_files, 1);
res_orig.P_SNR   = zeros(num_files, 1);
res_orig.Log_SNRAS = zeros(num_files, 1);
res_orig.Penalty   = zeros(num_files, 1);
res_orig.Psi       = zeros(num_files, 1);
res_orig.Time_ms   = zeros(num_files, 1);   

% Noisy transit results
res_noisy = res_orig;

% V-shaped eclipse results
res_vshape = res_orig;

fprintf('Starting advanced analysis of %d files...\n\n', num_files);

%% Main analysis loop
for k = 1:num_files
    filename = file_list(k).name;
    fprintf('=== Processing: %s ===\n', filename);
    
    try
        % --- 1. Load and clean the data ---
        data = readtable(filename, 'ReadVariableNames', false);
        time_seconds = data{:, 1};
        flux_raw = data{:, 2};
        
        % Remove outliers via 5-sigma clipping
        z_scores = (flux_raw - mean(flux_raw, 'omitnan')) / std(flux_raw, 'omitnan');
        clean_mask = abs(z_scores) < 5;
        time_seconds = time_seconds(clean_mask);
        flux_raw = flux_raw(clean_mask);
        
        % Normalise flux around 1
        flux = flux_raw / median(flux_raw, 'omitnan');
        
        % Convert time from seconds to days (needed for P-SNR binning)
        time_days = time_seconds / (24 * 3600);
        
        % --- 2. Build the transit mask (same for all scenarios) ---
        max_time = max(time_days);
        num_transits = floor((max_time - T0_DAYS) / PERIOD_DAYS) + 1;
        transit_mask = false(size(time_days));
        
        for i = 0:(num_transits - 1)
            current_t0 = T0_DAYS + i * PERIOD_DAYS;
            transit_mask = transit_mask | ...
                           (time_days >= current_t0 - DURATION_DAYS/2) & ...
                           (time_days <= current_t0 + DURATION_DAYS/2);
        end
        
        % =============================================
        % SCENARIO 1: Original clean data
        % =============================================
        flux_orig = flux;
        [trad_snr_orig, t_snr_orig, p_snr_orig, ...
         log_snras_orig, penalty_orig, psi_orig, time_orig] = ...
            compute_all_metrics_timed(flux_orig, transit_mask, time_days);
        
        fprintf('  [Original] Psi = %.4f | Log-SNRAS = %.2f | Time = %.2f ms\n', ...
            psi_orig, log_snras_orig, time_orig);
        
        % =============================================
        % SCENARIO 2: Add Gaussian noise inside transits
        % =============================================
        flux_noisy = flux;
        for i = 0:(num_transits - 1)
            current_t0 = T0_DAYS + i * PERIOD_DAYS;
            in_transit_now = (time_days >= current_t0 - DURATION_DAYS/2) & ...
                             (time_days <= current_t0 + DURATION_DAYS/2);
            flux_noisy(in_transit_now) = flux_noisy(in_transit_now) + ...
                                         NOISE_LEVEL * randn(sum(in_transit_now), 1);
        end
        [trad_snr_noisy, t_snr_noisy, p_snr_noisy, ...
         log_snras_noisy, penalty_noisy, psi_noisy, time_noisy] = ...
            compute_all_metrics_timed(flux_noisy, transit_mask, time_days);
        
        fprintf('  [Noisy  ] Psi = %.4f | Log-SNRAS = %.2f | Time = %.2f ms\n', ...
            psi_noisy, log_snras_noisy, time_noisy);
        
        % =============================================
        % SCENARIO 3: Inject V-shaped eclipsing binary
        % =============================================
        flux_vshape = flux;
        for i = 0:(num_transits - 1)
            current_t0 = T0_DAYS + i * PERIOD_DAYS;
            in_transit_now = (time_days >= current_t0 - DURATION_DAYS/2) & ...
                             (time_days <= current_t0 + DURATION_DAYS/2);
            t_in = time_days(in_transit_now);
            t_center = current_t0;
            depth_v = ECLIPSE_DEPTH * (1 - abs(t_in - t_center) / (DURATION_DAYS/2));
            flux_vshape(in_transit_now) = flux_vshape(in_transit_now) - depth_v;
        end
        [trad_snr_vshape, t_snr_vshape, p_snr_vshape, ...
         log_snras_vshape, penalty_vshape, psi_vshape, time_vshape] = ...
            compute_all_metrics_timed(flux_vshape, transit_mask, time_days);
        
        fprintf('  [V-shape] Psi = %.4f | Log-SNRAS = %.2f | Time = %.2f ms\n\n', ...
            psi_vshape, log_snras_vshape, time_vshape);
        
        % --- Store results for Original ---
        res_orig.Filename(k)   = string(filename);
        res_orig.Trad_SNR(k)   = trad_snr_orig;
        res_orig.T_SNR(k)      = t_snr_orig;
        res_orig.P_SNR(k)      = p_snr_orig;
        res_orig.Log_SNRAS(k)  = log_snras_orig;
        res_orig.Penalty(k)    = penalty_orig;
        res_orig.Psi(k)        = psi_orig;
        res_orig.Time_ms(k)    = time_orig;
        
        % --- Store results for Noisy ---
        res_noisy.Filename(k)  = string(filename);
        res_noisy.Trad_SNR(k)  = trad_snr_noisy;
        res_noisy.T_SNR(k)     = t_snr_noisy;
        res_noisy.P_SNR(k)     = p_snr_noisy;
        res_noisy.Log_SNRAS(k) = log_snras_noisy;
        res_noisy.Penalty(k)   = penalty_noisy;
        res_noisy.Psi(k)       = psi_noisy;
        res_noisy.Time_ms(k)   = time_noisy;
        
        % --- Store results for V-shape ---
        res_vshape.Filename(k)  = string(filename);
        res_vshape.Trad_SNR(k)  = trad_snr_vshape;
        res_vshape.T_SNR(k)     = t_snr_vshape;
        res_vshape.P_SNR(k)     = p_snr_vshape;
        res_vshape.Log_SNRAS(k) = log_snras_vshape;
        res_vshape.Penalty(k)   = penalty_vshape;
        res_vshape.Psi(k)       = psi_vshape;
        res_vshape.Time_ms(k)   = time_vshape;
        
    catch ME
        fprintf('Failed! %s\n', ME.message);
    end
end

%% Display and export results
fprintf('\n========== SPEED RESULTS (Log-SNRAS Core Only) ==========\n');
speed_table = table();
speed_table.Filename = res_orig.Filename;
speed_table.Time_Original_ms = res_orig.Time_ms;
speed_table.Time_Noisy_ms   = res_noisy.Time_ms;
speed_table.Time_VShape_ms  = res_vshape.Time_ms;
disp(speed_table);

fprintf('\n========== COMPARISON OF SCENARIOS (PSI VALUES) ==========\n');
comparison_psi = table();
comparison_psi.Filename = res_orig.Filename;
comparison_psi.Psi_Orig = res_orig.Psi;
comparison_psi.Psi_Noisy = res_noisy.Psi;
comparison_psi.Psi_VShape = res_vshape.Psi;
disp(comparison_psi);

% Export all results (including timing)
writetable(res_orig,   'PLATO_Original_Results.csv');
writetable(res_noisy,  'PLATO_Noisy_Results.csv');
writetable(res_vshape, 'PLATO_VShape_Results.csv');
writetable(speed_table, 'PLATO_Speed_Results.csv');

fprintf('\nResults exported to CSV files.\n');

%% ===== Local helper functions =====

function [trad_snr, t_snr, p_snr, log_snras, penalty, psi, elapsed_ms] = ...
    compute_all_metrics_timed(flux, transit_mask, time_days)
    % Compute all four metrics and measure Log-SNRAS core execution time.

    f_in  = flux(transit_mask);
    f_out = flux(~transit_mask);
    t_out = time_days(~transit_mask);
    
    sigma_out = std(f_out, 'omitnan');
    depth = 1 - median(f_in, 'omitnan');
    N_in = sum(transit_mask);
    
    trad_snr = (depth / sigma_out) * sqrt(N_in);
    t_snr    = trad_snr;
    p_snr    = compute_p_snr(f_out, t_out, depth, sigma_out, N_in);
    
    t_start = tic;
    [log_snras, penalty, psi] = calculate_log_snras(flux, trad_snr, transit_mask);
    elapsed_ms = toc(t_start) * 1000;   
end

function p_snr = compute_p_snr(f_out, t_out, depth, sigma_out, N_in)
    % P-SNR according to Pont et al. (2006), considering red noise.
    min_bin_days = 0.01;
    max_bin_days = 0.5;
    n_steps = 20;
    bin_sizes = logspace(log10(min_bin_days), log10(max_bin_days), n_steps);
    
    beta_values = zeros(size(bin_sizes));
    valid = false(size(bin_sizes));
    
    for i = 1:n_steps
        bin_days = bin_sizes(i);
        t_min = min(t_out);
        t_max = max(t_out);
        bin_edges = t_min:bin_days:t_max;
        if length(bin_edges) < 3, continue; end
        
        [~, ~, bin_idx] = histcounts(t_out, bin_edges);
        bin_means = accumarray(bin_idx(bin_idx>0), f_out(bin_idx>0), [], @mean);
        if length(bin_means) < 3, continue; end
        
        sigma_actual = std(bin_means, 'omitnan');
        avg_points_per_bin = length(f_out) / (length(bin_edges)-1);
        if avg_points_per_bin < 1, continue; end
        sigma_white = sigma_out / sqrt(avg_points_per_bin);
        
        beta_values(i) = sigma_actual / sigma_white;
        valid(i) = true;
    end
    
    if any(valid)
        beta = max(beta_values(valid));
        beta = max(beta, 1.0);
    else
        beta = 1.0;
    end
    
    sigma_red = beta * sigma_out;
    p_snr = (depth / sigma_red) * sqrt(N_in);
end