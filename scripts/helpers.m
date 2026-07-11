% =========================================================
%  Helper functions 
%  =========================================================

function [label, label_source, transit_expected] = get_label(name_lower, ...
    time_btjd, T0, Period, T_margin, artifact_kw, complex_kw, other_confirmed_kw)
    transit_expected = 0;
    if contains(name_lower, '261136679')
        if ~isempty(time_btjd) && any(~isnan(time_btjd))
            t_clean = time_btjd(~isnan(time_btjd));
            t_start = min(t_clean); t_end = max(t_clean);
            n_start = floor((t_start - T0) / Period);
            transit_times = T0 + (n_start-1:n_start+ceil((t_end-t_start)/Period)+1)' * Period;
            in_window = transit_times >= (t_start - T_margin) & transit_times <= (t_end + T_margin);
            transit_expected = double(any(in_window));
            if transit_expected == 1
                label = 1; label_source = 'ephemeris_confirmed';
            else
                label = 0; label_source = 'ephemeris_no_transit';
            end
        else
            label = NaN; label_source = 'name_no_time_excluded'; transit_expected = NaN;
        end
        return;
    end
    if any(cellfun(@(k) contains(name_lower, k), artifact_kw))
        label = 0; label_source = 'literature_artifact'; return;
    end
    if any(cellfun(@(k) contains(name_lower, k), complex_kw))
        label = NaN; label_source = 'complex_excluded'; return;
    end
    if any(cellfun(@(k) contains(name_lower, k), other_confirmed_kw))
        label = 1; label_source = 'literature_confirmed'; return;
    end
    label = NaN; label_source = 'unknown_excluded';
end

function results = compute_and_store(results, row, fName, sourceType, flux, ...
    label, label_source, transit_expected, t_start, t_end, tier1_cut, tier2_cut, ...
    time_clean, T0, Period, T_dur)
    median_flux = median(flux);
    f_norm = flux / median_flux;
    name_lower = lower(fName);
    is_known_planet = contains(name_lower, '261136679');
    
    if is_known_planet && ~isempty(time_clean) && length(time_clean) == length(flux)
        phase = mod(time_clean - T0, Period) / Period;
        phase(phase > 0.5) = phase(phase > 0.5) - 1;
        half_dur_phase = (T_dur/2) / Period;
        transit_mask = abs(phase) <= half_dur_phase;
        if sum(transit_mask) < 5 || sum(~transit_mask) < 5
            sigma_global = std(f_norm);
            transit_mask = f_norm < (1 - 3 * sigma_global);
            if sum(transit_mask) < 5
                [~, sortIdx] = sort(f_norm);
                num_points = ceil(0.01 * length(f_norm));
                transit_mask = false(size(f_norm));
                transit_mask(sortIdx(1:num_points)) = true;
            end
        end
    else
        sigma_global = std(f_norm);
        transit_mask = f_norm < (1 - 3 * sigma_global);
        if sum(transit_mask) < 5
            [~, sortIdx] = sort(f_norm);
            num_points = ceil(0.01 * length(f_norm));
            transit_mask = false(size(f_norm));
            transit_mask(sortIdx(1:num_points)) = true;
        end
    end
    
    f_out = f_norm(~transit_mask);
    f_in  = f_norm(transit_mask);
    N_in = length(f_in); N_out = length(f_out);
    
    if N_in == 0 || N_out == 0
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
    sigma_mad = 1.4826 * mad(f_out,1) + 1e-12;
    snr_robust = (depth / sigma_mad) * sqrt(N_in);
    f_sm = movmean(f_out, 12);
    snr_pont = depth / sqrt((sig_out^2/N_in) + std(f_sm)^2 + 1e-12);
    snr_bls = (depth^2 / var_out) * sqrt(N_in);
    
    % External function call
    [log_snras, penalty_term, psi] = calculate_log_snras(f_norm, snr_trad, transit_mask);
    penalty_pct = (penalty_term - 1) * 100;
    suppression_pct = (1 - log_snras / snr_trad) * 100;
    
    if penalty_pct <= tier1_cut, tier_label = 'Tier 1';
    elseif penalty_pct <= tier2_cut, tier_label = 'Tier 2';
    else, tier_label = 'Tier 3'; end
    
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
    results.Psi(row) = psi; results.Penalty_pct(row) = penalty_pct;
    results.Suppression_pct(row) = suppression_pct;
    results.Tier(row) = string(tier_label);
    results.Transit_Expected(row) = double(transit_expected);
    results.T_start_BTJD(row) = t_start;
    results.T_end_BTJD(row) = t_end;
end

function [p_value, z_score] = deLongExact(labels, scores1, scores2)
    labels = labels(:);
    scores1 = scores1(:);
    scores2 = scores2(:);
    pos = 1; neg = 0;
    idx_pos = find(labels == pos);
    idx_neg = find(labels == neg);
    n_pos = length(idx_pos);
    n_neg = length(idx_neg);
    function V = computeV(s)
        s_pos = s(idx_pos);
        s_neg = s(idx_neg);
        V = (s_pos > s_neg') + 0.5 * (s_pos == s_neg');
    end
    V1 = computeV(scores1);
    V2 = computeV(scores2);
    auc1 = mean(V1(:));
    auc2 = mean(V2(:));
    V1_pos_avg = mean(V1, 2);
    V2_pos_avg = mean(V2, 2);
    V1_neg_avg = mean(V1, 1)';
    V2_neg_avg = mean(V2, 1)';
    cov_pos = cov(V1_pos_avg, V2_pos_avg);
    cov_neg = cov(V1_neg_avg, V2_neg_avg);
    S = cov_pos / n_pos + cov_neg / n_neg;
    delta = auc1 - auc2;
    se = sqrt(S(1,1) + S(2,2) - 2*S(1,2));
    if se == 0
        if delta == 0, z_score = 0; p_value = 1.0;
        else, z_score = Inf; p_value = 0.0; end
    else
        z_score = delta / se;
        p_value = 2 * normcdf(-abs(z_score));
    end
end
