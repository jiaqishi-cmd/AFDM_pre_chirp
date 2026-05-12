% MAIN_FOCUSED_BESTCASE_BER_SNR
% 读取 focused Bemani/GPS 搜索得到的代表性 case，并做轻量 BER-SNR 验证。
% 默认选择 Phi(delta) 指标最差的 case，因为它更接近“近退秩”证据。
rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
setup_paths(rootDir);

outputDir = fullfile(fileparts(rootDir), 'results');
bestCasePath = fullfile(outputDir, 'focused_bemani_gps_best_cases.mat');

if ~exist(bestCasePath, 'file')
    fprintf('Focused best-case file not found. Running focused search first...\n');
    run(fullfile(rootDir, 'main_focused_bemani_gps_search.m'));
end

loaded = load(bestCasePath, 'focused_best_cases');

% Case A 的 key-equation mismatch 差距最大，更适合先看 BER 曲线是否拉开。
% 如需切回 Phi 指标最差的 Case D，把 case_kind 改成 'phi'。
case_kind = 'mismatch';
switch case_kind
    case 'mismatch'
        best_case = loaded.focused_best_cases.mismatch_case;
    case 'phi'
        best_case = loaded.focused_best_cases.phi_case;
    otherwise
        error('Unknown case_kind: %s', case_kind);
end

fprintf('========== Focused best-case BER-SNR ==========\n');
fprintf('Using %s case: Case %s, L=%d, l2=%d, alpha2=%d, pattern=%s, delta=%s\n', ...
    case_kind, best_case.case_name, best_case.L, best_case.l2, best_case.alpha2, ...
    best_case.gps_pattern, best_case.delta_type);
fprintf('Focused metrics: mean_E baseline %.4g, GPS %.4g | sigma_ratio baseline %.4g, GPS %.4g | col_corr baseline %.4g, GPS %.4g\n', ...
    best_case.mean_E_base, best_case.mean_E_GPS, ...
    best_case.sigma_ratio_base, best_case.sigma_ratio_GPS, ...
    best_case.col_corr_base, best_case.col_corr_GPS);

opts = struct();
opts.snr_values = 0:2:34;
opts.num_frames = 1000;
opts.seed = 20260508;
opts.gain_mode = 'fixed';
opts.fixed_theta = pi;
opts.label = sprintf('focused_%s_case_%s', case_kind, best_case.case_name);
results = run_bestcase_ber_snr(best_case, opts); %#ok<NASGU>

save(fullfile(outputDir, 'focused_bestcase_ber_snr.mat'), 'results', 'best_case');
