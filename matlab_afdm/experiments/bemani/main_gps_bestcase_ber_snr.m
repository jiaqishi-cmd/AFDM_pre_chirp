% MAIN_GPS_BESTCASE_BER_SNR
% 璇诲彇绗竴闃舵 best_case 骞惰繍琛?BER-SNR 楠岃瘉銆?
rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
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
