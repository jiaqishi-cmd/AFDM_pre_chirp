% MAIN_GPS_BESTCASE_BER_SNR
% 读取第一阶段 best_case 并运行 BER-SNR 验证。

rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
setup_paths(rootDir);

bestCasePath = fullfile(fileparts(rootDir), 'results', 'gps_unique_best_case.mat');
if ~exist(bestCasePath, 'file')
    error('Best case file not found. Run main_gps_unique_rank_loss_search first.');
end

loaded = load(bestCasePath, 'best_case');
opts.snr_values = 0:2:30;
opts.num_frames = 500;
results = run_bestcase_ber_snr(loaded.best_case, opts); %#ok<NASGU>
