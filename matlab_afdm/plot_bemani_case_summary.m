function plot_bemani_case_summary(resultsSource, outputDir)
%PLOT_BEMANI_CASE_SUMMARY Plot unambiguous Case A-D summary figures.
%   Each x-axis group is one named two-path case. For each case, the row
%   with the largest GPS improvement_mean is selected from results_table.

    if nargin < 1 || isempty(resultsSource)
        resultsSource = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
            'results', 'bemani_gps_key_equation_search.mat');
    end
    if nargin < 2 || isempty(outputDir)
        outputDir = fileparts(resultsSource);
    end

    loaded = load(resultsSource, 'results_table');
    results_table = loaded.results_table;

    caseNames = ["Case_A", "Case_B", "Case_C", "Case_D"];
    caseLabels = ["A: L=16", "B: L=32", "C: L=15", "D: L=17"];
    summary = select_case_rows(results_table, caseNames);

    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    plot_mean_error(summary, caseLabels, outputDir, timestamp);
    plot_phi_sigma_ratio(summary, caseLabels, outputDir, timestamp);
    plot_phi_col_corr(summary, caseLabels, outputDir, timestamp);
    write_summary_table(summary, outputDir, timestamp);
end

function summary = select_case_rows(results_table, caseNames)
    summary = table();
    for idx = 1:numel(caseNames)
        rows = results_table(results_table.case_id == caseNames(idx), :);
        if isempty(rows)
            error('No rows found for %s.', caseNames(idx));
        end
        [~, bestIdx] = max(rows.improvement_mean);
        summary = [summary; rows(bestIdx, :)]; %#ok<AGROW>
    end
end

function plot_mean_error(summary, caseLabels, outputDir, timestamp)
    values = [summary.mean_E_base, summary.mean_E_GPS, summary.mean_E_proposed, ...
        summary.mean_E_zero, summary.mean_E_ocdm];

    figure('Name', 'Bemani Case A-D mean equation error', 'Color', 'w');
    bar(values);
    grid on;
    set(gca, 'YScale', 'log');
    xticklabels(caseLabels);
    ylabel('Mean key-equation error');
    xlabel('Named two-path case');
    title('Bemani key-equation mismatch by case');
    legend({'Baseline', 'GPS', 'Proposed', 'c2=0', 'c2=1/(2N)'}, 'Location', 'northwest');
    saveas(gcf, fullfile(outputDir, ['bemani_case_summary_mean_error_' timestamp '.png']));
end

function plot_phi_sigma_ratio(summary, caseLabels, outputDir, timestamp)
    values = [summary.sigma_ratio_base, summary.sigma_ratio_GPS, summary.sigma_ratio_proposed];

    figure('Name', 'Bemani Case A-D Phi sigma ratio', 'Color', 'w');
    bar(values);
    grid on;
    xticklabels(caseLabels);
    ylabel('\sigma_{min}(\Phi) / \sigma_{max}(\Phi)');
    xlabel('Named two-path case');
    title('Phi(\delta) conditioning by case');
    legend({'Baseline', 'GPS', 'Proposed'}, 'Location', 'southwest');
    saveas(gcf, fullfile(outputDir, ['bemani_case_summary_phi_sigma_ratio_' timestamp '.png']));
end

function plot_phi_col_corr(summary, caseLabels, outputDir, timestamp)
    values = [summary.col_corr_base, summary.col_corr_GPS, summary.col_corr_proposed];

    figure('Name', 'Bemani Case A-D Phi column correlation', 'Color', 'w');
    bar(values);
    grid on;
    ylim([0, 1]);
    xticklabels(caseLabels);
    ylabel('|(H_1\delta)^H(H_2\delta)| / (||H_1\delta|| ||H_2\delta||)');
    xlabel('Named two-path case');
    title('Phi(\delta) column correlation by case');
    legend({'Baseline', 'GPS', 'Proposed'}, 'Location', 'northwest');
    saveas(gcf, fullfile(outputDir, ['bemani_case_summary_phi_col_corr_' timestamp '.png']));
end

function write_summary_table(summary, outputDir, timestamp)
    cols = {'case_id', 'l2', 'alpha2', 'L', 'gps_pattern_name', ...
        'proposed_pattern_name', 'delta_type', 'mean_E_base', 'mean_E_GPS', ...
        'mean_E_proposed', 'sigma_ratio_base', 'sigma_ratio_GPS', ...
        'sigma_ratio_proposed', 'col_corr_base', 'col_corr_GPS', ...
        'col_corr_proposed'};
    writetable(summary(:, cols), ...
        fullfile(outputDir, ['bemani_case_summary_' timestamp '.csv']));
end
