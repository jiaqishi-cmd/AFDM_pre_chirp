function plot_caseA_besttheta_ber(ber_results, theta, outputDir)
%PLOT_CASEA_BESTTHETA_BER 姹囨€荤粯鍒?best theta 鐨?BPSK/QPSK BER銆?
    if nargin < 3 || isempty(outputDir)
        outputDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results');
    end
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    fields = fieldnames(ber_results);
    for idx = 1:numel(fields)
        result = ber_results.(fields{idx});
        figure('Name', ['Case A BER ' fields{idx}], 'Color', 'w');
        semilogy(result.snr_values, result.ber, 'o-', 'LineWidth', 2);
        grid on;
        xlabel('SNR (dB)');
        ylabel('BER');
        title(sprintf('Case A %s fixed theta=%.6f pi', upper(fields{idx}), theta / pi));
        legend(result.schemes, 'Location', 'southwest');
        saveas(gcf, fullfile(outputDir, sprintf('caseA_besttheta_%s_ber_%s.png', fields{idx}, timestamp)));
    end
end
