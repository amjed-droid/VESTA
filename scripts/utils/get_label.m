function [label, label_source, transit_expected] = get_label(name_lower, ...
    time_btjd, T0, Period, T_margin, artifact_kw, complex_kw, other_confirmed_kw)
%% GET_LABEL Assigns independent ground-truth label based on literature/ephemeris.

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
