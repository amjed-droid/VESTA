%% GENERATE_TABLE2 Generates Table 2: Abridged benchmark catalog (plain text output).
%
%   GENERATE_TABLE2() reads the evaluation_dataset_v2.csv and extracts
%   representative targets across the three classification tiers (Tier 1, 
%   Tier 2, and Tier 3) as presented in the manuscript's Table 2.
%   It handles missing values (displayed as '---') and marks excluded
%   ambiguous targets with a dagger symbol (†).
%
%   METHODOLOGY
%   -----------
%   1. Loads the unified evaluation catalog.
%   2. Filters the required targets.
%   3. Computes short readable names.
%   4. Replaces NaN values with '---' for display.
%   5. Prints a plain text table to the Command Window (no LaTeX markup).
%   6. Exports the results to Table2_AbridgedCatalog.csv.
%
%   DEPENDENCIES
%   ------------
%   evaluation_dataset_v2.csv : Must be in the current directory.
%
%   SEE ALSO
%   --------
%   main_pipeline, generate_table5
function generate_table2(target_list, outfile)
    % Set default arguments
    if nargin < 1
        target_list = {'sector_9', 'sector_1', 'kepler_1625', ...
                       '12557548', 'wd_1145', 'toi_201', ...
                       'v723_mon', '5812701', 'quarter_10', '12644769'};
    end
    if nargin < 2
        outfile = 'Table2_AbridgedCatalog.csv';
    end

    % --- 1. Load the evaluation dataset ---
    csv_file = 'evaluation_dataset_v2.csv';
    if ~isfile(csv_file)
        error('File %s not found. Please ensure it exists in the current directory.', csv_file);
    end
    data = readtable(csv_file);
    data.Filename = string(data.Filename);
    data.Tier = string(data.Tier);

    % --- 2. Search for each target and format names ---
    table2_rows = table();
    fprintf('\n=== Extracting Table 2 Abridged Catalog ===\n');
    
    for k = 1:length(target_list)
        kw = target_list{k};
        mask = contains(lower(data.Filename), lower(kw));
        sub = data(mask, :);
        if height(sub) > 0
            row = sub(1, :);
            
            % Generate short names (plain text, no LaTeX)
            if contains(row.Filename, 'sector_9')
                row.ShortName = string('TIC 261136679 (S9)');
            elseif contains(row.Filename, 'sector_1')
                row.ShortName = string('TIC 261136679 (S1)');
            elseif contains(row.Filename, 'Kepler_1625')
                row.ShortName = string('Kepler-1625');
            elseif contains(row.Filename, '12557548')
                row.ShortName = string('KIC 12557548†');
            elseif contains(row.Filename, 'wd_1145')
                row.ShortName = string('WD 1145+017†');
            elseif contains(row.Filename, 'TOI')
                row.ShortName = string('TOI-201');
            elseif contains(row.Filename, 'V723')
                row.ShortName = string('V723 Mon');
            elseif contains(row.Filename, 'quarter_10')
                row.ShortName = string('KIC 8462852 (Q10)');
            else
                % Extract prefix before dot for other KIC targets
                fname = char(row.Filename);
                [~, id] = strtok(fname, '.');
                if isempty(id) || isequal(id, '.fits')
                    id = fname;
                else
                    id = extractBefore(fname, '.');
                end
                row.ShortName = string(id);
            end
            
            table2_rows = [table2_rows; row];
            fprintf('  [%s] -> Matched: %s (Tier %s)\n', ...
                char(kw), char(row.ShortName), char(row.Tier));
        else
            fprintf('  [%s] -> No match found\n', char(kw));
        end
    end

    % --- 3. Display formatted table in Command Window (plain text) ---
    if height(table2_rows) == 0
        warning('No rows matched the target list. Output file will be empty.');
    else
        % Print header (no LaTeX)
        fprintf('\n');
        fprintf(repmat('=', 1, 115));
        fprintf('\n');
        fprintf('                         TABLE 2: ABRIDGED BENCHMARK CATALOG\n');
        fprintf(repmat('=', 1, 115));
        fprintf('\n');
        fprintf('%-28s %-8s %-8s %-8s %-8s %-10s %-6s\n', ...
            'Target', 'T-SNR', 'R-SNR', 'P-SNR', 'B-SNR', 'L-SNRAS', 'Tier');
        fprintf(repmat('-', 1, 115));
        fprintf('\n');

        % Group by Tier order: Tier 1 -> Tier 2 -> Tier 3
        tier_titles = {'Tier 1 (P <= 0.15)', ...
                       'Tier 2 (0.15 < P <= 0.60)', ...
                       'Tier 3 (P > 0.60) --- Veto'};

        for t = 1:3
            % Determine filter condition based on P (Penalty_pct/100)
            if t == 1
                filt = (table2_rows.Penalty_pct / 100) <= 0.15;
            elseif t == 2
                filt = ((table2_rows.Penalty_pct / 100) > 0.15) & ((table2_rows.Penalty_pct / 100) <= 0.60);
            else
                filt = (table2_rows.Penalty_pct / 100) > 0.60;
            end
            
            sub_tier = table2_rows(filt, :);
            if height(sub_tier) > 0
                fprintf('\n--- %s ---\n', tier_titles{t});
                fprintf(repmat('-', 1, 115));
                fprintf('\n');
                
                for i = 1:height(sub_tier)
                    row = sub_tier(i, :);
                    
                    % Handle missing values (NaN -> '---')
                    if isnan(row.P_SNR), p_snr_str = '---   '; else, p_snr_str = sprintf('%-8.2f', row.P_SNR); end
                    if isnan(row.B_SNR), b_snr_str = '---   '; else, b_snr_str = sprintf('%-8.2f', row.B_SNR); end
                    
                    % Remove extra "Tier " from display for the numerical column
                    tier_display = char(row.Tier);
                    tier_display = strrep(tier_display, 'Tier ', '');
                    
                    fprintf('%-28s %-8.2f %-8.2f %s %s %-10.2f %-6s\n', ...
                        char(row.ShortName), ...
                        row.T_SNR, row.R_SNR, ...
                        p_snr_str, b_snr_str, ...
                        row.L_SNRAS, tier_display);
                end
            end
        end
        fprintf('\n');
        fprintf(repmat('=', 1, 115));
        fprintf('\n');

        % --- 4. Export to CSV ---
        disp_cols = {'ShortName','T_SNR','R_SNR','P_SNR','B_SNR','L_SNRAS','Tier'};
        writetable(table2_rows(:, disp_cols), outfile);
        fprintf('\n✅ Table 2 saved to: %s\n', outfile);
    end
end
