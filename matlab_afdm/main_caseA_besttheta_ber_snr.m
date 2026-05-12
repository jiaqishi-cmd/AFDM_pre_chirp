% MAIN_CASEA_BESTTHETA_BER_SNR
% 读取 Case A 相位扫描结果，用 best theta 和 theta=pi 分别跑 BPSK/QPSK BER。
rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
setup_paths(rootDir);

outputDir = fullfile(fileparts(rootDir), 'results');
scanPath = fullfile(outputDir, 'caseA_phase_scan_dmin_results.mat');
if ~exist(scanPath, 'file')
    fprintf('Phase scan result not found. Running main_caseA_phase_scan_dmin first...\n');
    run(fullfile(rootDir, 'main_caseA_phase_scan_dmin.m'));
end

loaded = load(scanPath, 'best_case', 'best_pattern', 'phase_scan_table');
bestTheta = loaded.best_case.theta;
gpsPattern = loaded.best_pattern;

opts = struct();
opts.snr_values = 0:2:34;
opts.num_frames = 2000;
opts.modulation_list = [2 4];
opts.seed = 20260508;

fprintf('========== Case A best-theta BER ==========\n');
fprintf('Best theta %.8f rad = %.6f*pi\n', bestTheta, bestTheta / pi);
opts.label_prefix = 'caseA_besttheta';
besttheta_results = run_ber_snr_fixed_channel_caseA(bestTheta, gpsPattern, opts); %#ok<NASGU>

fprintf('\n========== Case A theta=pi BER control ==========\n');
opts.label_prefix = 'caseA_theta_pi';
theta_pi_results = run_ber_snr_fixed_channel_caseA(pi, gpsPattern, opts); %#ok<NASGU>

comparison_table = build_theta_comparison_table(loaded.phase_scan_table, bestTheta);
save(fullfile(outputDir, 'caseA_besttheta_ber_snr_results.mat'), ...
    'besttheta_results', 'theta_pi_results', 'comparison_table', 'bestTheta', 'gpsPattern');
writetable(comparison_table, fullfile(outputDir, 'caseA_besttheta_vs_pi_comparison.csv'));

plot_caseA_besttheta_ber(besttheta_results, bestTheta, outputDir);
plot_caseA_besttheta_ber(theta_pi_results, pi, outputDir);
disp(comparison_table);

function comparison = build_theta_comparison_table(T, bestTheta)
    [~, bestIdx] = min(abs(T.theta - bestTheta));
    [~, piIdx] = min(abs(T.theta - pi));
    rows = T([bestIdx, piIdx], :);
    theta_name = ["best_theta"; "theta_pi"];
    comparison = table(theta_name, rows.theta, rows.theta_over_pi, ...
        rows.dmin2_base, rows.dmin2_GPS, rows.dmin2_prop, ...
        rows.gap_base_GPS_dB, rows.gap_prop_GPS_dB, rows.score, ...
        rows.delta_type_GPS, ...
        'VariableNames', {'theta_name','theta','theta_over_pi','dmin2_base','dmin2_GPS', ...
        'dmin2_prop','gap_base_GPS_dB','gap_prop_GPS_dB','score','delta_type_GPS'});
end
