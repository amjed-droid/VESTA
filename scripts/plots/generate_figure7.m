%% =========================================================================
%  Script: generate_figure7_toi201_diagnostic.m
%  Purpose: Diagnostic decomposition of in-transit dispersion for TOI-201,
%           isolating the deterministic transit-shape component from
%           genuine photometric noise (Section 3.5 of the manuscript).
%
%  This script demonstrates that TOI-201's misclassification as Tier 3
%  (Veto) by the isolated logarithmic penalty is driven by the transit's
%  own geometric shape (ingress-bottom-egress curvature) inflating
%  sigma_in, rather than genuine heteroscedastic noise instability.
%  A LOESS-smoothed shape fit is subtracted from the in-transit points,
%  and sigma_in is recomputed from the residuals only.
%
%  INPUT
%  -----
%  TOI_201.fits : raw TESS light curve (binary table), expected fields:
%                 Field 1 = TIME (BTJD), Field 2 = FLUX (SAP or PDCSAP)
%
%  OUTPUT
%  ------
%  Figure7_TOI201_ShapeDiagnostic.png : 600 DPI, two-panel diagnostic figure
%
%  DEPENDENCIES
%  ------------
%  Requires the Curve Fitting Toolbox (for the smooth() function with the
%  'loess' method).
%
%  AUTHOR
%  ------
%  Ahmed Sattar Jabbar
% =========================================================================

clear; clc; close all;

%% --- 1. Load raw TESS light curve ---
fits_path = fullfile('data_folder', 'TOI_201.fits');  % adjust path as needed
flux_data = fitsread(fits_path, 'binarytable');

time = double(flux_data{1});
flux = double(flux_data{2});

valid = ~isnan(flux) & ~isnan(time);
time  = time(valid);
flux  = flux(valid);
f_norm = flux / median(flux);

%% --- 2. Identify in-transit window via 3-sigma auto-detection ---
% (Consistent with the transit-mask logic used throughout the main
% pipeline; see process_files.m)
sigma_global  = std(f_norm);
transit_mask  = f_norm < (1 - 3*sigma_global);

time_in = time(transit_mask);
f_in    = f_norm(transit_mask);

[time_in_sorted, sortIdx] = sort(time_in);
f_in_sorted = f_in(sortIdx);

%% --- 3. Fit and remove the deterministic transit shape (LOESS) ---
win_size = max(3, round(0.15 * length(f_in_sorted)));
if mod(win_size, 2) == 0
    win_size = win_size + 1;  % smooth() requires an odd window size
end
shape_fit = smooth(f_in_sorted, win_size, 'loess');
residuals = f_in_sorted - shape_fit;

%% --- 4. Report raw vs. shape-corrected sigma_in (for Table in Section 3.5) ---
sigma_in_raw       = std(f_in_sorted);
sigma_in_corrected = std(residuals);
fprintf('Raw sigma_in:              %.6f\n', sigma_in_raw);
fprintf('Shape-corrected sigma_in:  %.6f\n', sigma_in_corrected);
fprintf('Reduction: %.1f%%\n', 100*(1 - sigma_in_corrected/sigma_in_raw));

%% --- 5. Generate two-panel diagnostic figure ---
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 9, 4]);

subplot(1,2,1);
plot(time_in_sorted, f_in_sorted, '.', 'Color', [0.8 0.1 0.1], 'MarkerSize', 10);
hold on;
plot(time_in_sorted, shape_fit, 'b-', 'LineWidth', 2);
title('(a) In-Transit Points and Fitted Transit Shape', 'FontSize', 11);
xlabel('Time (BTJD)');
ylabel('Normalized Flux');
legend('In-transit data', 'LOESS shape fit', 'Location', 'best');
grid on; box on;

subplot(1,2,2);
plot(time_in_sorted, residuals, 'k.', 'MarkerSize', 10);
hold on;
yline(0, 'g--', 'LineWidth', 1.5);
title(sprintf('(b) Residuals After Shape Removal (\\sigma=%.5f)', sigma_in_corrected), ...
    'FontSize', 11);
xlabel('Time (BTJD)');
ylabel('Residual Flux');
grid on; box on;

sgtitle('TOI-201: Isolating Geometric Shape from Genuine Photometric Noise', ...
    'FontSize', 13, 'FontWeight', 'bold');

%% --- 6. Export high-resolution figure for manuscript/repository ---
exportgraphics(fig, 'Figure7_TOI201_ShapeDiagnostic.png', 'Resolution', 600);
fprintf('\nFigure saved successfully as Figure7_TOI201_ShapeDiagnostic.png\n');
