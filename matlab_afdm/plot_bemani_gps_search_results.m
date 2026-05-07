function plot_bemani_gps_search_results(results_table, outputDir, timestamp, phi_plot)
%PLOT_BEMANI_GPS_SEARCH_RESULTS Generate summary plots for key-equation search.

    if nargin < 2 || isempty(outputDir)
        outputDir = pwd;
    end
    if nargin < 3 || isempty(timestamp)
        timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    end
    if nargin < 4
        phi_plot = struct();
    end

    figure('Name', 'Bemani key equation mean error', 'Color', 'w');
    loglog(results_table.mean_E_base, results_table.mean_E_GPS, '.', 'MarkerSize', 6);
    hold on;
    if ismember('mean_E_proposed', results_table.Properties.VariableNames)
        loglog(results_table.mean_E_base, results_table.mean_E_proposed, '.', 'MarkerSize', 6);
    end
    finiteVals = [results_table.mean_E_base; results_table.mean_E_GPS];
    finiteVals = finiteVals(isfinite(finiteVals) & finiteVals > 0);
    if isempty(finiteVals)
        finiteVals = [1e-16; 1];
    end
    low = min(finiteVals);
    high = max(finiteVals);
    loglog([low high], [low high], 'k--', 'LineWidth', 1.5);
    grid on;
    xlabel('mean E baseline');
    ylabel('mean E compared scheme');
    title('Bemani key-equation mean error');
    if ismember('mean_E_proposed', results_table.Properties.VariableNames)
        legend({'GPS', 'Proposed', 'y=x'}, 'Location', 'best');
    end
    saveas(gcf, fullfile(outputDir, ['bemani_gps_mean_error_scatter_' timestamp '.png']));

    top = top_by(results_table, 'improvement_mean', 'descend', 10);
    figure('Name', 'Top mean-error comparison', 'Color', 'w');
    if ismember('mean_E_proposed', top.Properties.VariableNames)
        bar([top.mean_E_base, top.mean_E_GPS, top.mean_E_proposed, top.mean_E_zero, top.mean_E_ocdm]);
        legend({'baseline', 'GPS', 'proposed', 'zero', 'OCDM-like'}, 'Location', 'best');
    else
        bar([top.mean_E_base, top.mean_E_GPS, top.mean_E_zero, top.mean_E_ocdm]);
        legend({'baseline', 'GPS', 'zero', 'OCDM-like'}, 'Location', 'best');
    end
    grid on;
    ylabel('Mean equation error');
    title('Top suspicious cases by improvement mean');
    xticklabels(top.case_id);
    xtickangle(45);
    saveas(gcf, fullfile(outputDir, ['bemani_gps_top_mean_error_bar_' timestamp '.png']));

    phi_rows = results_table(~isnan(results_table.sigma_ratio_base) & ~isnan(results_table.sigma_ratio_GPS), :);
    if ~isempty(phi_rows)
        top_phi = top_by_ratio(phi_rows, 10);
        figure('Name', 'Phi sigma ratio comparison', 'Color', 'w');
        if ismember('sigma_ratio_proposed', top_phi.Properties.VariableNames)
            bar([top_phi.sigma_ratio_base, top_phi.sigma_ratio_GPS, top_phi.sigma_ratio_proposed]);
            legend({'baseline', 'GPS', 'proposed'}, 'Location', 'best');
        else
            bar([top_phi.sigma_ratio_base, top_phi.sigma_ratio_GPS]);
            legend({'baseline', 'GPS'}, 'Location', 'best');
        end
        grid on;
        ylabel('sigma min / sigma max');
        title('Phi(delta) sigma-ratio comparison');
        xticklabels(compact_case_labels(top_phi));
        xtickangle(45);
        saveas(gcf, fullfile(outputDir, ['bemani_gps_phi_sigma_ratio_bar_' timestamp '.png']));
    end

    if isfield(phi_plot, 'phi1') && isfield(phi_plot, 'phi2') && ~isempty(phi_plot.phi1)
        figure('Name', 'GPS Phi column proportionality', 'Color', 'w');
        subplot(3, 1, 1);
        plot(abs(phi_plot.phi1), 'LineWidth', 1.5);
        grid on;
        ylabel('|phi1|');
        title(phi_plot.title, 'Interpreter', 'none');

        subplot(3, 1, 2);
        plot(abs(phi_plot.phi2), 'LineWidth', 1.5);
        grid on;
        ylabel('|phi2|');

        subplot(3, 1, 3);
        plot(angle(phi_plot.phi1 ./ phi_plot.phi2), 'LineWidth', 1.5);
        grid on;
        ylabel('angle(phi1/phi2)');
        xlabel('Index');
        saveas(gcf, fullfile(outputDir, ['bemani_gps_phi_columns_' timestamp '.png']));
    end
end

function labels = compact_case_labels(T)
    labels = strings(height(T), 1);
    for idx = 1:height(T)
        deltaText = char(T.delta_type(idx));
        deltaText = regexprep(deltaText, '^delta\d+_', '');
        deltaText = regexprep(deltaText, '_random_', '_r');
        if strlength(deltaText) > 18
            deltaText = extractBefore(deltaText, 19);
        end
        labels(idx) = sprintf('%s L%d %s %s', ...
            T.case_id(idx), T.L(idx), T.gps_pattern_name(idx), string(deltaText));
    end
end

function top = top_by(T, fieldName, direction, n)
    [~, order] = sort(T.(fieldName), direction);
    order = order(1:min(n, numel(order)));
    top = T(order, :);
end

function top = top_by_ratio(T, n)
    ratio = T.sigma_ratio_GPS ./ max(T.sigma_ratio_base, eps);
    [~, order] = sort(ratio, 'ascend');
    order = order(1:min(n, numel(order)));
    top = T(order, :);
end
