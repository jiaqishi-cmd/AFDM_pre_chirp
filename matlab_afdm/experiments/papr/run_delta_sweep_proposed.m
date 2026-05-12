% RUN_DELTA_SWEEP_PROPOSED
% 鍒嗘瀽 proposed 涓夌偣鎵板姩闆嗗悎 {c2-delta, c2, c2+delta} 涓?delta 澶у皬鏃讹紝
% 瀵?PAPR銆佷俊閬撴棤鍏崇粨鏋勯闄╁拰 Case A fixed-channel BER 鐨勫奖鍝嶃€?
rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
addpath(rootDir);
setup_paths(rootDir);

% ========================
% 鍙皟鍙傛暟
% ========================
rng(1, 'twister');
delta_ratio_list = [0.01 0.02 0.05 0.1 0.2 0.5 1.0];
numPaprFrames = 1000;
paprTargetCcdf = 1e-3;
paprThresholds = 0:0.05:12;
SNR_dB = 28;
berMaxBits = 64000;
berMinErrTarget = 100;
theta_caseA = pi;
seedBase = 20260509;

% ========================
% 鍩烘湰绯荤粺鍙傛暟
% ========================
baseConfig = afdm_config();
N = baseConfig.waveform.NumSubcarriers;
V = baseConfig.pre_chirp.num_groups;
c2_base = baseConfig.pre_chirp.base_c2;
c1_caseA = 7 / (2 * 64);
gpsPattern = [2 2 1 1]; %#ok<NASGU> 鐢ㄤ簬缁撴瀯鍙傝€冿紝涓嶅奖鍝?greedy GPS PAPR銆?
metricsCfg = struct();
metricsCfg.N = N;
metricsCfg.M = N;
metricsCfg.modType = 'BPSK';
metricsCfg.phase_tol = pi / 24;
metricsCfg.num_bins = 64;
metricsCfg.struct_weights = [0.3 0.2 0.5];

outputDir = fullfile(fileparts(rootDir), 'results');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

fprintf('========== Proposed delta sweep ==========\n');
fprintf('N=%d, V=%d, c2=%.6g, PAPR frames=%d, BER SNR=%g dB\n', ...
    N, V, c2_base, numPaprFrames, SNR_dB);

% ========================
% 鍙傝€冿細baseline 鍜?GPS 鐨?PAPR 涓?BER
% ========================
fprintf('\nSampling baseline/GPS PAPR references...\n');
[paprBaseSamples, ~] = sample_scheme_papr(baseConfig, 'baseline', numPaprFrames, seedBase);
[paprGpsSamples, gpsSelection] = sample_scheme_papr(baseConfig, 'paper_grouping', numPaprFrames, seedBase);
papr_base_at_target = papr_at_ccdf(paprBaseSamples, paprTargetCcdf);
papr_gps_at_target = papr_at_ccdf(paprGpsSamples, paprTargetCcdf);

fprintf('Running baseline/GPS Case A BER references...\n');
ber_base_ref = run_caseA_scheme_ber(c2_base, SNR_dB, theta_caseA, c1_caseA, berMaxBits, berMinErrTarget, seedBase + 11);
gpsProfile = build_pre_chirp_profile('paper_grouping', N, baseConfig.pre_chirp);
c2_gps_ref = build_pattern_c2_from_candidate_set(gpsProfile.candidate_set, gpsProfile.group_index, mode_selection(gpsSelection, V));
ber_gps_ref = run_caseA_scheme_ber(c2_gps_ref, SNR_dB, theta_caseA, c1_caseA, berMaxBits, berMinErrTarget, seedBase + 12);

% ========================
% delta sweep 涓诲惊鐜?% ========================
numDelta = numel(delta_ratio_list);
papr_prop_at_target = zeros(numDelta, 1);
ber_prop = zeros(numDelta, 1);
ber_prop_err = zeros(numDelta, 1);
ber_prop_bits = zeros(numDelta, 1);
R_struct = zeros(numDelta, 1);
R_phase = zeros(numDelta, 1);
alignment_ratio = zeros(numDelta, 1);
diagonal_perturb_dist = zeros(numDelta, 1);
R_dev = zeros(numDelta, 1);
select_minus_ratio = zeros(numDelta, 1);
select_zero_ratio = zeros(numDelta, 1);
select_plus_ratio = zeros(numDelta, 1);
mean_papr_prop = zeros(numDelta, 1);
p99_papr_prop = zeros(numDelta, 1);

selectedPatterns = cell(numDelta, 1);
paprPropSamplesAll = zeros(numPaprFrames, numDelta);

for deltaIdx = 1:numDelta
    deltaRatio = delta_ratio_list(deltaIdx);
    delta = deltaRatio * c2_base;
    cfg = baseConfig;
    cfg.pre_chirp.delta = delta;

    fprintf('\nDelta ratio %.4g (delta=%.6g)\n', deltaRatio, delta);
    [paprSamples, selectionMat, c2SelectedList] = sample_scheme_papr(cfg, 'proposed_grouping', numPaprFrames, seedBase + 100 * deltaIdx);
    paprPropSamplesAll(:, deltaIdx) = paprSamples;
    selectedPatterns{deltaIdx} = selectionMat;
    papr_prop_at_target(deltaIdx) = papr_at_ccdf(paprSamples, paprTargetCcdf);
    mean_papr_prop(deltaIdx) = mean(paprSamples);
    p99_papr_prop(deltaIdx) = percentile_by_sort(paprSamples, 0.99);

    [select_minus_ratio(deltaIdx), select_zero_ratio(deltaIdx), select_plus_ratio(deltaIdx)] = ...
        candidate_selection_ratios(selectionMat);

    metricRows = compute_selected_structural_metrics(c2SelectedList, c2_base * ones(N, 1), metricsCfg);
    R_struct(deltaIdx) = mean([metricRows.R_struct]);
    R_phase(deltaIdx) = mean([metricRows.R_phase]);
    alignment_ratio(deltaIdx) = mean([metricRows.alignment_ratio]);
    diagonal_perturb_dist(deltaIdx) = mean([metricRows.diagonal_perturb_dist]);
    R_dev(deltaIdx) = mean([metricRows.R_dev]);

    representativeSelection = mode_selection(selectionMat, V);
    propProfile = build_pre_chirp_profile('proposed_grouping', N, cfg.pre_chirp);
    c2PropBer = build_pattern_c2_from_candidate_set(propProfile.candidate_set, propProfile.group_index, representativeSelection);
    [ber_prop(deltaIdx), ber_prop_err(deltaIdx), ber_prop_bits(deltaIdx)] = ...
        run_caseA_scheme_ber(c2PropBer, SNR_dB, theta_caseA, c1_caseA, berMaxBits, berMinErrTarget, seedBase + 200 * deltaIdx);

    fprintf('PAPR@1e-3 %.3f dB | R_struct %.3f | BER %.3e (%d/%d) | selection [- 0 +]=[%.2f %.2f %.2f]\n', ...
        papr_prop_at_target(deltaIdx), R_struct(deltaIdx), ber_prop(deltaIdx), ...
        ber_prop_err(deltaIdx), ber_prop_bits(deltaIdx), ...
        select_minus_ratio(deltaIdx), select_zero_ratio(deltaIdx), select_plus_ratio(deltaIdx));
end

% ========================
% 淇濆瓨鍜岀粯鍥?% ========================
results_table = table(delta_ratio_list(:), papr_prop_at_target, mean_papr_prop, p99_papr_prop, ...
    R_struct, R_phase, alignment_ratio, diagonal_perturb_dist, R_dev, ...
    ber_prop, ber_prop_err, ber_prop_bits, ...
    select_minus_ratio, select_zero_ratio, select_plus_ratio, ...
    'VariableNames', {'delta_ratio','papr_at_ccdf_1e3','mean_papr','p99_papr', ...
    'R_struct','R_phase','alignment_ratio','diagonal_perturb_dist','R_dev', ...
    'ber_prop','ber_prop_err','ber_prop_bits', ...
    'select_minus_ratio','select_zero_ratio','select_plus_ratio'});

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
save(fullfile(outputDir, 'results_delta_sweep.mat'), ...
    'results_table', 'delta_ratio_list', 'papr_base_at_target', 'papr_gps_at_target', ...
    'ber_base_ref', 'ber_gps_ref', 'paprBaseSamples', 'paprGpsSamples', ...
    'paprPropSamplesAll', 'selectedPatterns', 'SNR_dB', 'theta_caseA', ...
    'paprTargetCcdf', 'numPaprFrames', 'berMaxBits');
writetable(results_table, fullfile(outputDir, ['delta_sweep_proposed_summary_' timestamp '.csv']));

plot_delta_sweep_results(results_table, delta_ratio_list, papr_base_at_target, papr_gps_at_target, ...
    ber_base_ref, ber_gps_ref, SNR_dB, outputDir, timestamp);

fprintf('\nSaved delta sweep results to %s\n', fullfile(outputDir, 'results_delta_sweep.mat'));
disp(results_table);

function [paprSamples, selectionMat, c2SelectedList] = sample_scheme_papr(baseConfig, scheme, numFrames, seedBase)
    N = baseConfig.waveform.NumSubcarriers;
    bitsPerFrame = N * log2(baseConfig.modulation.M_mod);
    cfg0 = apply_pre_chirp_scheme(baseConfig, scheme);
    numGroups = max(cfg0.pre_chirp.profile.group_index);
    paprSamples = zeros(numFrames, 1);
    selectionMat = ones(numFrames, numGroups);
    c2SelectedList = zeros(N, numFrames);

    for frameIdx = 1:numFrames
        rng(seedBase + frameIdx, 'twister');
        txBits = randi([0, 1], bitsPerFrame, 1);
        cfg = apply_pre_chirp_scheme(baseConfig, scheme);
        cfg.tx.bits = txBits;
        [~, papr, ~, txState] = afdm_tx_engine(cfg);
        paprSamples(frameIdx) = papr;
        c2SelectedList(:, frameIdx) = expand_c2(txState.c2, N);
        if isfield(txState.pre_chirp_profile, 'selection') && ...
                isfield(txState.pre_chirp_profile.selection, 'selected_candidate_index')
            selectionMat(frameIdx, :) = txState.pre_chirp_profile.selection.selected_candidate_index(:).';
        end
    end
end

function c2Vec = expand_c2(c2, N)
    if isscalar(c2)
        c2Vec = c2 * ones(N, 1);
    else
        c2Vec = c2(:);
    end
end

function value = papr_at_ccdf(samples, targetCcdf)
    % CCDF=P(PAPR>x)銆倀arget=1e-3 瀵瑰簲 99.9 percentile銆?    value = percentile_by_sort(samples, 1 - targetCcdf);
end

function value = percentile_by_sort(samples, q)
    sorted = sort(samples(:));
    idx = max(1, min(numel(sorted), ceil(q * numel(sorted))));
    value = sorted(idx);
end

function [ber, err, bits] = run_caseA_scheme_ber(c2Vec, snrDb, theta, c1, maxBits, minErrTarget, seedBase)
    N = 64;
    cfg = afdm_config();
    cfg.waveform.NumSubcarriers = N;
    cfg.waveform.CPPLength = 2;
    cfg.waveform.c1 = c1;
    cfg.waveform.c2 = c2Vec;
    cfg.modulation.M_mod = 2;
    cfg.modulation.modType = 'psk';
    cfg.channel.multipath = true;
    cfg.channel.add_noise = true;
    cfg.channel.snr_db = snrDb;
    cfg.channel.delay_taps = [0, 2];
    cfg.channel.doppler_taps = [0, 2];
    cfg.channel.doppler_freq = cfg.channel.doppler_taps / N;
    cfg.channel.chan_coef = [1, exp(1i * theta)] / sqrt(2);
    cfg.pre_chirp.scheme = 'baseline';
    cfg.pre_chirp.profile.scheme = 'baseline';
    cfg.pre_chirp.profile.c2 = c2Vec;

    err = 0;
    bits = 0;
    frameIdx = 0;
    while err < minErrTarget && bits < maxBits
        frameIdx = frameIdx + 1;
        frame = simulate_frame(cfg, seedBase + frameIdx);
        err = err + frame.err_bits;
        bits = bits + frame.total_bits;
    end
    ber = err / max(bits, 1);
end

function metrics = compute_selected_structural_metrics(c2SelectedList, c2BaseVec, cfg)
    numFrames = size(c2SelectedList, 2);
    maxMetricFrames = min(numFrames, 200);
    sampleIdx = round(linspace(1, numFrames, maxMetricFrames));
    m = (0:size(c2SelectedList, 1)-1).';
    phiBase = exp(1i * 2 * pi .* c2BaseVec(:) .* (m.^2));
    metrics = repmat(struct('R_struct', 0, 'R_phase', 0, 'alignment_ratio', 0, ...
        'diagonal_perturb_dist', 0, 'R_dev', 0), 1, numel(sampleIdx));
    for idx = 1:numel(sampleIdx)
        c2Vec = c2SelectedList(:, sampleIdx(idx));
        item = afdm.metrics.c2_structural(c2Vec, cfg);
        phi = exp(1i * 2 * pi .* c2Vec(:) .* (m.^2));
        diagonalDist = mean(abs(phi - phiBase));
        metrics(idx).R_struct = item.R_struct;
        metrics(idx).R_phase = item.phase_degeneracy_risk;
        metrics(idx).alignment_ratio = item.constellation_alignment_ratio;
        metrics(idx).diagonal_perturb_dist = diagonalDist;
        metrics(idx).R_dev = diagonalDist / 2; % 褰掍竴鍒?[0,1] 闄勮繎鐨勭浉浣?mask 鍋忕椋庨櫓銆?    end
end

function [minusRatio, zeroRatio, plusRatio] = candidate_selection_ratios(selectionMat)
    total = numel(selectionMat);
    % proposed_profile 褰撳墠鍊欓€夐『搴忎负 [c2, c2-delta, c2+delta]銆?    zeroRatio = sum(selectionMat(:) == 1) / total;
    minusRatio = sum(selectionMat(:) == 2) / total;
    plusRatio = sum(selectionMat(:) == 3) / total;
end

function selection = mode_selection(selectionMat, numGroups)
    selection = ones(1, numGroups);
    for groupIdx = 1:numGroups
        selection(groupIdx) = mode(selectionMat(:, groupIdx));
    end
end

function c2Vec = build_pattern_c2_from_candidate_set(candidateSet, groupIndex, selection)
    c2Vec = candidateSet(:, 1);
    for groupIdx = 1:numel(selection)
        indices = groupIndex == groupIdx;
        c2Vec(indices) = candidateSet(indices, selection(groupIdx));
    end
end

function plot_delta_sweep_results(T, deltaRatios, paprBase, paprGps, berBase, berGps, snrDb, outputDir, timestamp)
    figure('Name', 'Delta sweep PAPR@CCDF', 'Color', 'w');
    semilogx(deltaRatios, T.papr_at_ccdf_1e3, 'o-', 'LineWidth', 2); hold on;
    yline(paprBase, '--', 'Baseline', 'LineWidth', 1.4);
    yline(paprGps, '--', 'GPS', 'LineWidth', 1.4);
    grid on;
    xlabel('\delta / c_2');
    ylabel('PAPR at CCDF=10^{-3} (dB)');
    title('Proposed delta sweep: PAPR tail');
    legend({'Proposed', 'Baseline ref', 'GPS ref'}, 'Location', 'best');
    saveas(gcf, fullfile(outputDir, ['delta_sweep_papr_' timestamp '.png']));

    figure('Name', 'Delta sweep structural risk', 'Color', 'w');
    semilogx(deltaRatios, T.R_struct, 'o-', 'LineWidth', 2); hold on;
    semilogx(deltaRatios, T.alignment_ratio, 's-', 'LineWidth', 2);
    semilogx(deltaRatios, T.R_dev, '^-', 'LineWidth', 2);
    grid on;
    xlabel('\delta / c_2');
    ylabel('Metric value');
    title('Proposed delta sweep: structural metrics');
    legend({'R\_struct', 'alignment ratio', 'R\_dev'}, 'Location', 'best');
    saveas(gcf, fullfile(outputDir, ['delta_sweep_structural_' timestamp '.png']));

    figure('Name', 'Delta sweep fixed BER', 'Color', 'w');
    semilogy(deltaRatios, max(T.ber_prop, 0.5 ./ max(T.ber_prop_bits, 1)), 'o-', 'LineWidth', 2); hold on;
    yline(max(berBase, min(T.ber_prop(T.ber_prop > 0), [], 'omitnan') / 10), '--', 'Baseline', 'LineWidth', 1.4);
    yline(max(berGps, min(T.ber_prop(T.ber_prop > 0), [], 'omitnan') / 10), '--', 'GPS', 'LineWidth', 1.4);
    grid on;
    xlabel('\delta / c_2');
    ylabel(sprintf('BER at SNR=%g dB', snrDb));
    title('Proposed delta sweep: Case A fixed-channel BER');
    legend({'Proposed', 'Baseline ref', 'GPS ref'}, 'Location', 'best');
    saveas(gcf, fullfile(outputDir, ['delta_sweep_ber_' timestamp '.png']));

    figure('Name', 'Delta sweep candidate selection', 'Color', 'w');
    bar(deltaRatios, [T.select_minus_ratio, T.select_zero_ratio, T.select_plus_ratio], 'stacked');
    set(gca, 'XScale', 'log');
    grid on;
    xlabel('\delta / c_2');
    ylabel('Selection ratio');
    title('Proposed greedy candidate selection ratios');
    legend({'c2-\delta', 'c2', 'c2+\delta'}, 'Location', 'eastoutside');
    saveas(gcf, fullfile(outputDir, ['delta_sweep_selection_' timestamp '.png']));
end
