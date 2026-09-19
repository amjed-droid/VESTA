function results = compute_and_store(results, row, fName, sourceType, flux, ...
    label, label_source, transit_expected, t_start, t_end, tier1_cut, tier2_cut, ...
    time_clean, T0, Period, T_dur)
%% COMPUTE_AND_STORE Computes candidate vetting metrics and populates catalog table.

    median_flux = median(flux);
    f_norm = flux / median_flux;
    name_lower = lower(fName);
    
    % Uniform Ephemeris Phasing
    if ~isempty(time_clean) && length(time_clean) == length(flux) && ~isnan(T0) && ~isnan(Period)
        phase = mod(time_clean - T0, Period) / Period;
        phase(phase > 0.5) = phase(phase > 0.5) - 1;
        half_dur_phase = (T_dur/2) / Period;
        transit_mask = abs(phase) <= half_dur_phase;
    else
        % Default fallback for unknown ephemeris
        transit_mask = false(size(f_norm));
    end
    
    f_out = f_norm(~transit_mask);
    f_in  = f_norm(transit_mask);
    N_in = length(f_in); N_out = length(f_out);
    
    if N_in < 4 || N_out < 10
        results.Filename(row) = string(fName);
        results.SourceType(row) = string(sourceType);
        results.Label(row) = label;
        results.Label_Source(row) = string(label_source);
        results.Transit_Expected(row) = double(transit_expected);
        results.T_start_BTJD(row) = t_start;
        results.T_end_BTJD(row) = t_end;
        results.N_total(row) = length(flux);
        results.T_SNR(row) = NaN; results.L_SNRAS(row) = NaN;
        results.Tier(row) = "N/A";
        return;
    end
    
    depth = abs(mean(f_out) - mean(f_in));
    var_out = var(f_out) + 1e-12;
    sig_out = std(f_out) + 1e-12;
    snr_trad = (depth / sig_out) * sqrt(N_in);
    
    % Robust SNR
    sigma_mad = 1.4826 * mad(f_out, 1) + 1e-12;
    snr_robust = (depth / sigma_mad) * sqrt(N_in);
    
    % Pont SNR
    f_sm = movmean(f_out, 12);
    snr_pont = depth / sqrt((sig_out^2/N_in) + std(f_sm)^2 + 1e-12);
    
    % BLS Proxy
    snr_bls = (depth^2 / var_out) * sqrt(N_in);
    
    % Log-SNRAS with transit shape correction
    if N_in >= 6
        x_in = linspace(-1, 1, N_in)';
        p_fit = polyfit(x_in, f_in, 2);
        f_model = polyval(p_fit, x_in);
        sig_in = std(f_in - f_model) + 1e-12;
    else
        sig_in = std(f_in) + 1e-12;
    end
    
    psi = abs(sig_in - sig_out) / sig_out;
    penalty = log(1 + psi);
    log_snras = snr_trad / (1 + penalty);
    suppression_pct = (1 - log_snras / snr_trad) * 100;
    
    % Tier assignment
    if penalty <= 0.15
        tier_label = 'Tier 1';
    elseif penalty <= 0.60
        tier_label = 'Tier 2';
    else
        tier_label = 'Tier 3';
    end
    
    results.Filename(row) = string(fName);
    results.SourceType(row) = string(sourceType);
    results.Label(row) = label;
    results.Label_Source(row) = string(label_source);
    results.N_total(row) = length(flux);
    results.N_in(row) = N_in; results.N_out(row) = N_out;
    results.Depth_ppm(row) = depth * 1e6;
    results.T_SNR(row) = snr_trad; results.R_SNR(row) = snr_robust;
    results.P_SNR(row) = snr_pont; results.B_SNR(row) = snr_bls;
    results.L_SNRAS(row) = log_snras;
    results.Psi(row) = psi; results.Penalty_pct(row) = penalty;
    results.Suppression_pct(row) = suppression_pct;
    results.Tier(row) = string(tier_label);
    results.Transit_Expected(row) = double(transit_expected);
    results.T_start_BTJD(row) = t_start;
    results.T_end_BTJD(row) = t_end;
end
