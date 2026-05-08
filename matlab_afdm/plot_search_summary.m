function plot_search_summary(results_table, outputDir, timestamp)
%PLOT_SEARCH_SUMMARY 绘制 GPS 独有近退秩搜索摘要图。

    if nargin < 2 || isempty(outputDir)
        outputDir = pwd;
    end
    if nargin < 3 || isempty(timestamp)
        timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    end

    valid = ~isnan(results_table.sigma_ratio_GPS);
    T = results_table(valid, :);
    if isempty(T)
        warning('No Phi-verified rows available for plotting.');
        return;
    end

    [~, order] = sort(T.score1, 'descend');
    top = T(order(1:min(20, numel(order))), :);

    figure('Name', 'GPS unique rank-loss score', 'Color', 'w');
    bar([top.sigma_ratio_base, top.sigma_ratio_GPS, top.sigma_ratio_prop]);
    grid on;
    ylabel('\sigma_{min}(\Phi)/\sigma_{max}(\Phi)');
    title('Top candidates by baseline/proposed vs GPS sigma-ratio score');
    legend({'Baseline', 'GPS', 'Proposed'}, 'Location', 'best');
    xticklabels(compact_labels(top));
    xtickangle(45);
    saveas(gcf, fullfile(outputDir, ['gps_unique_sigma_ratio_top_' timestamp '.png']));

    figure('Name', 'GPS column correlation top candidates', 'Color', 'w');
    bar([top.col_corr_base, top.col_corr_GPS, top.col_corr_prop]);
    grid on;
    ylim([0, 1]);
    ylabel('Column correlation of \Phi(\delta)');
    title('Column correlation on top GPS near-rank-loss candidates');
    legend({'Baseline', 'GPS', 'Proposed'}, 'Location', 'best');
    xticklabels(compact_labels(top));
    xtickangle(45);
    saveas(gcf, fullfile(outputDir, ['gps_unique_col_corr_top_' timestamp '.png']));
end

function labels = compact_labels(T)
    labels = strings(height(T), 1);
    for idx = 1:height(T)
        labels(idx) = sprintf('N%d V%d L%d %s', T.N(idx), T.V(idx), T.L(idx), T.delta_type(idx));
    end
end
