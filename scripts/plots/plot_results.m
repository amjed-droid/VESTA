function plot_results(y_true_v, scores_inv, raw_auc, boot_pen, auc_penalty, n_valid, raw_scores_cell)
%% PLOT_RESULTS Generates the ROC and Bootstrap distribution figures.
%
%   PLOT_RESULTS(Y_TRUE_V, SCORES_INV, RAW_AUC, BOOT_PEN, AUC_PENALTY, ...
%                N_VALID, RAW_SCORES_CELL) produces the two primary
%   diagnostic figures for the Log-SNRAS manuscript:
%
%       1. Multi-class ROC curve comparing the five raw SNR metrics with
%          the isolated penalty (-Penalty_pct) and random chance.
%       2. Bootstrap distribution histogram of the penalty-based AUC,
%          overlaid with the point estimate and appropriate axes labels.
%
%   METHODOLOGY
%   -----------
%   ROC Figure:
%       - The `perfcurve` function is used to generate (FPR, TPR) pairs.
%       - The raw metrics are plotted with differentiating colors and
%         line styles (dashed, dash-dot, dotted) to ensure clarity in
%         both color and grayscale publication formats.
%       - The isolated penalty is highlighted with a thicker green line.
%       - The legend reports the computed point AUC for each metric.
%
%   Bootstrap Distribution Figure:
%       - A 50-bin histogram of the bootstrapped AUC values (BOOT_PEN).
%       - A vertical red dashed line marks the point estimate AUC_PENALTY.
%       - Axes are labeled to reflect the Bootstrap AUC of the Penalty metric.
%
%   All figures are exported via `exportgraphics` at 600 DPI resolution
%   to ensure publication-quality vector/raster output.
%
%   INPUTS
%   ------
%   y_true_v        : Ground-truth labels for the valid subset.
%   scores_inv      : Negated penalty scores (-Penalty_pct).
%   raw_auc         : Point AUC of the five raw metrics.
%   boot_pen        : Bootstrap distribution of penalty AUC (B x 1).
%   auc_penalty     : Point AUC of the isolated penalty.
%   n_valid         : Number of segments in the valid subset.
%   raw_scores_cell : Cell array of raw scores for the five SNR metrics.
%
%   OUTPUTS
%   -------
%   PNG Figures: ROC_EphemerisGroundTruth.png,
%                Bootstrap_PenaltyDistribution.png
%
%   SEE ALSO
%   --------
%   main_pipeline, analyze_performance, perfcurve, exportgraphics

    raw_methods_plot = {'T-SNR','R-SNR','P-SNR','B-SNR','Log-SNRAS'};
    raw_colors = {[0.8,0.2,0.2],[0.2,0.2,0.7],[0.5,0.5,0.5],[0.3,0.3,0.3],[0.9,0.5,0.1]};
    raw_styles = {'--','-.',':', '--','-.'};

    % --- FIGURE: ROC ---
    figure('Color','w','Units','inches','Position',[1,1,8,7]);
    hold on; box on; grid on;
    
    [X_pen, Y_pen] = perfcurve(y_true_v, scores_inv, 1);
    h_pen = plot(X_pen, Y_pen, 'Color',[0,0.6,0.2], 'LineWidth',3.5);
    
    h_raw_plot = zeros(length(raw_methods_plot),1);
    for m = 1:length(raw_methods_plot)
        [X, Y] = perfcurve(y_true_v, raw_scores_cell{m}, 1);
        h_raw_plot(m) = plot(X, Y, 'Color', raw_colors{m}, 'LineStyle', raw_styles{m}, 'LineWidth', 1.5);
    end
    
    h_chance = plot([0 1],[0 1],'k:','LineWidth',1.2);
    
    all_handles = [h_pen; h_raw_plot; h_chance];
    all_names = {
        sprintf('Penalty $-\\mathcal{P}$ (AUC = %.3f)', auc_penalty);
        sprintf('T-SNR (AUC = %.3f)', raw_auc(1));
        sprintf('R-SNR (AUC = %.3f)', raw_auc(2));
        sprintf('P-SNR (AUC = %.3f)', raw_auc(3));
        sprintf('B-SNR (AUC = %.3f)', raw_auc(4));
        sprintf('Log-SNRAS (AUC = %.3f)', raw_auc(5));
        'Random Chance (0.500)'
    };
    
    legend(all_handles, all_names, 'Location','southeast', 'Interpreter','latex', 'FontSize',11);
    set(gca,'FontSize',13,'LineWidth',1.2);
    xlabel('False Positive Rate','FontSize',15,'FontWeight','bold');
    ylabel('True Positive Rate','FontSize',15,'FontWeight','bold');
    title(sprintf('ROC: Ephemeris/Literature Ground Truth (N=%d)', n_valid), 'FontSize',16,'FontWeight','bold');
    axis square; hold off;
    exportgraphics(gcf,'ROC_EphemerisGroundTruth.png','Resolution',600);
    fprintf('\nSaved: ROC_EphemerisGroundTruth.png\n');

    % --- Bootstrap distribution plot for Penalty AUC ---
    figure('Color','w','Units','inches','Position',[1,1,6,4.5]);
    histogram(boot_pen, 50, 'FaceColor',[0,0.6,0.2],'EdgeColor','none');
    hold on;
    xline(auc_penalty, 'r--', 'LineWidth', 2, 'Label', sprintf('Point estimate = %.3f', auc_penalty));
    xlabel('Bootstrap AUC (Penalty)', 'FontSize', 13);
    ylabel('Frequency', 'FontSize', 13);
    title('Bootstrap Distribution: Penalty-based AUC', 'FontSize', 14);
    grid on;
    exportgraphics(gcf, 'Bootstrap_PenaltyDistribution.png', 'Resolution', 600);
    fprintf('Saved: Bootstrap_PenaltyDistribution.png\n');
end
