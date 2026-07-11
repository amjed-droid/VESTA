%% GENERATE_BASELINE_FIGURE Creates the baseline validation figure for Pi Mensae.
%
%   GENERATE_BASELINE_FIGURE(SNR_TRAD, SNR_LOG, PENALTY, OUTPUT_FILENAME)
%   generates a two-panel figure illustrating the baseline validation and
%   metric transparency for Pi Mensae (Sector 1).
%
%   The top panel displays a stable homoscedastic light curve with a
%   shallow planetary transit, overlaid with a blue transit model and
%   labeled to indicate the noise characteristics.
%   The bottom panel presents a bar chart comparing the Traditional SNR
%   versus the Log-SNRAS stability score, with the metric transparency
%   (penalty percentage) noted in the subtitle.
%
%   METHODOLOGY
%   -----------
%   1. Synthetic Light Curve: A realistic out-of-transit baseline is
%      generated with low homoscedastic Gaussian noise. An in-transit
%      dip is added using a smooth exponential kernel to mimic a
%      planetary transit.
%   2. Visual Consistency: The phase range (0.3 to 0.7), noise scaling,
%      and transit depth are tuned to match the visual aesthetics of
%      the manuscript's Figure 5 while accepting user-defined statistics.
%   3. Bar Chart: Two bars are rendered using user-supplied numerical
%      values. The bars are colored (gray and blue) to be clearly
%      distinguishable even in grayscale printing.
%   4. Export: The figure is exported as a high-resolution PNG file
%      (600 DPI) using exportgraphics.
%
%   INPUTS
%   ------
%   snr_trad        : (Scalar) Traditional SNR value for the bar chart.
%                     Default: 40.63.
%   snr_log         : (Scalar) Log-SNRAS value for the bar chart.
%                     Default: 36.58.
%   penalty         : (Scalar) Penalty percentage for the subtitle.
%                     Default: 11.1.
%   output_filename : (String) Name of the output PNG file.
%                     Default: 'Baseline_Validation_Pi_Mensae.png'.
%
%   OUTPUTS
%   -------
%   PNG Figure : Saved to the specified OUTPUT_FILENAME at 600 DPI.
%
%   SEE ALSO
%   --------
%   exportgraphics, scatter, bar
function generate_baseline_figure(snr_trad, snr_log, penalty, output_filename)
    % Set default values if not provided
    if nargin < 1, snr_trad = 40.63; end
    if nargin < 2, snr_log  = 36.58; end
    if nargin < 3, penalty  = 11.1;  end
    if nargin < 4, output_filename = 'Baseline_Validation_Pi_Mensae.png'; end

    % --- 1. Generate Synthetic Light Curve Data ---
    rng(123); % Fix random seed for reproducibility

    % Generate out-of-transit points (Phase: 0.3 to 0.45 and 0.55 to 0.7)
    n_out = 600;
    x_out1 = 0.3 + 0.15 * rand(n_out/2, 1);
    x_out2 = 0.55 + 0.15 * rand(n_out/2, 1);
    x_out = [x_out1; x_out2];

    % Generate in-transit points (Phase: 0.45 to 0.55)
    n_in = 150;
    x_in = 0.45 + 0.10 * rand(n_in, 1);

    % Merge phases
    x_phase = [x_out; x_in];

    % Generate flux with homoscedastic noise and transit dip
    sigma_noise = 0.0006; 
    flux_out = 1.0 + sigma_noise * randn(n_out, 1);
    flux_in  = 0.9965 + sigma_noise * randn(n_in, 1);
    flux_vals = [flux_out; flux_in];

    % Generate the blue transit model (smooth curve)
    x_model = linspace(0.3, 0.7, 1000);
    y_model = 1 - 0.0035 * exp(-((x_model - 0.5) / 0.02).^6);

    % --- 2. Create the Figure ---
    figure('Color', 'w', 'Units', 'inches', 'Position', [2, 2, 8, 10]);

    % ==================== Top Panel: Light Curve ====================
    subplot(2, 1, 1);
    hold on; box on; grid on;

    % 1. Plot scattered data points
    scatter(x_phase, flux_vals, 12, 'k', 'filled', 'MarkerFaceAlpha', 0.6);

    % 2. Plot the blue transit model
    plot(x_model, y_model, 'b-', 'LineWidth', 2.5);

    % 3. Add explanatory text
    text(0.31, 0.9966, 'Homoscedastic Noise ($\sigma_{in} \approx \sigma_{out}$)', ...
        'FontSize', 11, 'FontWeight', 'bold', 'Color', 'k', 'Interpreter', 'latex');

    % Format axes
    xlim([0.3, 0.7]);
    ylim([0.995, 1.002]);
    xlabel('Phase', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Normalized Flux', 'FontSize', 12, 'FontWeight', 'bold');
    title('(a) Stable Baseline: Pi Mensae (Sector 1)', 'FontSize', 14, 'FontWeight', 'bold');

    % ==================== Bottom Panel: Bar Chart ====================
    subplot(2, 1, 2);
    hold on; box on; grid on;

    % Plot bars with user-defined values
    b = bar([snr_trad, snr_log], 0.6);

    % Apply colors (Gray and Blue)
    b.FaceColor = 'flat';
    b.CData = [0.5, 0.5, 0.5; 0, 0.45, 0.74];

    % Add text labels above bars
    text(1, snr_trad + (max([snr_trad, snr_log])*0.05), sprintf('%.2f', snr_trad), ...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
    text(2, snr_log + (max([snr_trad, snr_log])*0.05), sprintf('%.2f', snr_log), ...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');

    % Format axes (Fixed ylim to [0, 50] as per your request)
    ylim([0, 50]); 
    set(gca, 'XTickLabel', {'Traditional SNR', 'Log-SNRAS'}, 'XTick', 1:2, 'FontSize', 12);
    ylabel('Detection Score ($\sigma$)', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
    title(sprintf('(b) Metric Transparency (Penalty: %.1f%%)', penalty), ...
        'FontSize', 14, 'FontWeight', 'bold');

    % --- 3. Overall Title and Export ---
    sgtitle('Baseline Validation using Pi Mensae (Sector 1)', 'FontSize', 17, 'FontWeight', 'bold');

    % Export with 600 DPI
    exportgraphics(gcf, output_filename, 'Resolution', 600);
    fprintf('Baseline validation figure saved as: %s\n', output_filename);
end
