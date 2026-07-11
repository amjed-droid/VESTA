% =========================================================================
% Script: plot_3d_surface.m
% Purpose: Regenerate the 3D Theoretical Response Surface for Log-SNRAS
% =========================================================================
clear; clc; close all;

% 1. Define Parameters for the Theoretical Surface
% Assuming a normalized transit depth (delta) of 1 and N_in = 100 for visualization
delta = 1.0;
N_in = 100;

% 2. Create the Grid for In-Transit and Out-of-Transit Sigma
[S_O, S_I] = meshgrid(linspace(0.5, 3, 50), linspace(0.5, 3, 50));

% 3. Apply the Log-SNRAS Mathematical Formulation
% P_surf: Logarithmic Penalty based on dispersion contrast
P_surf = log(1 + abs(S_I - S_O) ./ S_O);

% Z_surf: The final variance-stabilized Log-SNRAS score
Z_surf = (delta * sqrt(N_in)) ./ (S_O .* (1 + P_surf));

% 4. Visualization Setup
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 7, 6]);

% Plot the surface
surf(S_O, S_I, Z_surf, 'EdgeColor', 'none');
colormap(parula); 

% Configure Colorbar
c = colorbar;
c.Label.String = 'Log-SNRAS Score';
c.Label.FontSize = 11;

% Configure Axes and Labels
xlabel('\sigma_{out} (Out-of-Transit Noise)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('\sigma_{in} (In-Transit Noise)', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('Detection Score', 'FontSize', 12, 'FontWeight', 'bold');
title('Log-SNRAS Response Surface', 'FontSize', 14, 'FontWeight', 'bold');

% Set viewing angle for optimal 3D perspective
view(135, 30);
grid on;

% 5. Export high-resolution image for publication/repository
exportgraphics(fig, 'Log_SNRAS_Response_Surface.png', 'Resolution', 600);
fprintf('Figure saved successfully as Log_SNRAS_Response_Surface.png\n');