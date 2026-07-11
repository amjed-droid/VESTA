%% GENERATE_SCATTER_PLOT_SCRIPT
%  Standalone script to generate the metric stability scatter plot.
%  This script reads 'evaluation_dataset_v2.csv' directly and produces
%  a log-log scatter plot with three theoretical tier boundaries.
%
%  METHODOLOGY
%  -----------
%  1. Loads the evaluation dataset from the CSV file.
%  2. Classifies segments into Confirmed Planets (Label=1),
%     Heteroscedastic Artifacts (Label=0), and Excluded (NaN).
%  3. Overlays the Unity Line (L = T), Tier 1/2 boundary (P=0.15),
%     and Tier 2/3 boundary (P=0.60) to illustrate penalty compression.
%  4. Renders the scatter plot with log-log scaling and exports
%     the figure as a high-resolution PNG (600 DPI).
%
%  SEE ALSO
%  --------
%  readtable, scatter, exportgraphics
clear; clc; close all;

%% --- 1. LOAD DATA ---
% The script automatically looks for the file in the current directory.
csv_file = 'evaluation_dataset_v2.csv';
if ~isfile(csv_file)
    error('File %s not found. Please ensure it exists in the current directory.', csv_file);
end
data = readtable(csv_file);

% Extract required variables
T_SNR   = data.T_SNR;
L_SNRAS = data.L_SNRAS;
Label   = data.Label;

%% --- 2. CLASSIFY SEGMENTS ---
% Ground truth classification based on the master pipeline labels
is_confirmed = (Label == 1);
is_artifact  = (Label == 0);
is_excluded  = isnan(Label);

%% --- 3. CREATE THE FIGURE ---
figure('Color', 'w', 'Position', [100, 100, 900, 750]);
hold on; grid on; box on;

% (a) Unity Line (L = T) - Black dashed
plot([0.1 2000], [0.1 2000], 'k--', 'LineWidth', 1.5, ...
    'DisplayName', 'Unity Line ($L = T$)');

% (b) Tier 1/2 Boundary (P = 0.15) -> L = T / 1.15 - Blue dotted
T_range = logspace(-1, 4, 100);
plot(T_range, T_range/1.15, 'b:', 'LineWidth', 1.5, ...
    'DisplayName', 'Tier 1/2 ($\mathcal{P}=0.15$)');

% (c) Tier 2/3 Boundary (P = 0.60) -> L = T / 1.60 - Red dash-dot
plot(T_range, T_range/1.60, 'r-.', 'LineWidth', 1.5, ...
    'DisplayName', 'Tier 2/3 ($\mathcal{P}=0.60$)');

% (d) Plot Data Points
% Excluded (Complex/Ambiguous) - Gray triangles
scatter(T_SNR(is_excluded), L_SNRAS(is_excluded), ...
    60, '^', 'MarkerFaceColor', [0.7 0.7 0.7], ...
    'MarkerEdgeColor', 'k', 'LineWidth', 0.5, ...
    'DisplayName', 'Excluded (Complex/Ambiguous)');

% Confirmed Planets - White circles
scatter(T_SNR(is_confirmed), L_SNRAS(is_confirmed), ...
    80, 'o', 'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'Confirmed Planets');

% Heteroscedastic Artifacts - Gray squares
scatter(T_SNR(is_artifact), L_SNRAS(is_artifact), ...
    80, 's', 'MarkerFaceColor', [0.5 0.5 0.5], ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'Heteroscedastic Artifacts');

%% --- 4. FORMAT AXES AND TITLES ---
set(gca, 'XScale', 'log', 'YScale', 'log');

% Add labels and title (English)
xlabel('Traditional SNR (T-SNR)', 'FontSize', 13);
ylabel('Log-SNRAS Stability Score (L-SNRAS)', 'FontSize', 13);
title('Metric Stability Mapping with Tiered Vetting Boundaries', ...
    'FontSize', 14, 'FontWeight', 'bold');

% Adjust axis limits dynamically
xlim([0.9, max(T_SNR(~isnan(T_SNR))) * 1.3]);
ylim([0.9, max(L_SNRAS(~isnan(L_SNRAS))) * 1.3]);

% Add legend with LaTeX interpreter
legend('Location', 'northwest', 'FontSize', 10, 'Interpreter', 'latex');

% Style adjustments
ax = gca;
ax.TickDir = 'in';
ax.LineWidth = 1.2;
ax.FontSize = 11;

%% --- 5. EXPORT FIGURE ---
output_filename = 'fig5_scatter_FINAL.png';
exportgraphics(gcf, output_filename, 'Resolution', 600);
fprintf('Figure saved successfully as: %s\n', output_filename);
