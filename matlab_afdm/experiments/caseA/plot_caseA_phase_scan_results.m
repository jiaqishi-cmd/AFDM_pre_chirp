function plot_caseA_phase_scan_results(scan, outputDir)
%PLOT_CASEA_PHASE_SCAN_RESULTS 缁樺埗 Case A 鐩镐綅鎵弿缁撴灉銆?
    if nargin < 2 || isempty(outputDir)
        outputDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results');
    end
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    T = scan.phase_scan_table;
    best = scan.best_case;
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');

    figure('Name', 'Case A dmin2 vs theta', 'Color', 'w');
    plot(T.theta, 10 * log10(T.dmin2_base + eps), 'LineWidth', 1.8); hold on;
    plot(T.theta, 10 * log10(T.dmin2_GPS + eps), 'LineWidth', 1.8);
    plot(T.theta, 10 * log10(T.dmin2_prop + eps), 'LineWidth', 1.8);
    xline(best.theta, '--k', 'best \theta', 'LineWidth', 1.2);
    grid on;
    xlabel('\theta (rad)');
    ylabel('d_{min}^2 (dB)');
    title('Case A finite-delta minimum distance vs path-gain phase');
    legend({'Baseline', 'GPS', 'Proposed'}, 'Location', 'best');
    saveas(gcf, fullfile(outputDir, ['caseA_phase_scan_dmin2_' timestamp '.png']));

    figure('Name', 'Case A GPS/base dmin ratio', 'Color', 'w');
    plot(T.theta, 10 * log10(T.ratio_GPS_base + eps), 'LineWidth', 1.8); hold on;
    yline(0, ':k', 'LineWidth', 1.1);
    xline(best.theta, '--k', 'best \theta', 'LineWidth', 1.2);
    grid on;
    xlabel('\theta (rad)');
    ylabel('10log10(dmin2 GPS / dmin2 baseline)');
    title('Case A GPS minimum-distance disadvantage vs baseline');
    saveas(gcf, fullfile(outputDir, ['caseA_phase_scan_ratio_gps_base_' timestamp '.png']));

    figure('Name', 'Case A dmin gaps', 'Color', 'w');
    plot(T.theta, T.gap_base_GPS_dB, 'LineWidth', 1.8); hold on;
    plot(T.theta, T.gap_prop_GPS_dB, 'LineWidth', 1.8);
    xline(best.theta, '--k', 'best \theta', 'LineWidth', 1.2);
    grid on;
    xlabel('\theta (rad)');
    ylabel('Gap over GPS (dB)');
    title('Case A baseline/proposed distance gap over GPS');
    legend({'Baseline-GPS', 'Proposed-GPS'}, 'Location', 'best');
    saveas(gcf, fullfile(outputDir, ['caseA_phase_scan_gaps_' timestamp '.png']));

    h1 = 1 / sqrt(2);
    h2 = exp(1i * best.theta) / sqrt(2);
    delta = scan.best_delta_GPS;
    Hbase = h1 * scan.paths.base.H1 + h2 * scan.paths.base.H2;
    Hgps = h1 * scan.paths.gps.H1 + h2 * scan.paths.gps.H2;
    Hprop = h1 * scan.paths.prop.H1 + h2 * scan.paths.prop.H2;

    figure('Name', 'Case A best-theta delta response', 'Color', 'w');
    plot(abs(Hbase * delta), 'LineWidth', 1.5); hold on;
    plot(abs(Hgps * delta), 'LineWidth', 1.5);
    plot(abs(Hprop * delta), 'LineWidth', 1.5);
    grid on;
    xlabel('Subcarrier index');
    ylabel('|H_{tot}\delta_{GPS-best}|');
    title(sprintf('Case A response on GPS-best delta, theta=%.4f rad', best.theta));
    legend({'Baseline', 'GPS', 'Proposed'}, 'Location', 'best');
    saveas(gcf, fullfile(outputDir, ['caseA_besttheta_delta_response_' timestamp '.png']));
end
