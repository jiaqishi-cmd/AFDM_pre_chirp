% RUN_CHOI_STYLE_RANDOM_CHANNEL_BER
% 鍙傝€?Choi 璁烘枃涓父鐢ㄧ殑 M=64銆丩cpp=8 AFDM 璁剧疆锛岄噰鐢ㄩ殢鏈?P-path
% doubly selective Rayleigh channel锛屾瘮杈?baseline/GPS/proposed 鐨勫钩鍧?BER銆?%
% 鏈疄楠屽亣璁炬帴鏀剁 perfect CSI锛岀敤浜庨殧绂讳笉鍚?pre-chirp / c2 pattern
% 璁捐鏈韩鐨勫奖鍝嶃€傚彂灏勭涓嶄娇鐢?CSI锛屽彧鏍规嵁褰撳墠鏁版嵁鍋?PAPR selection銆?% GPS/proposed 鐨?selected pattern 鏆傛椂鍋囪閫氳繃鐞嗘兂 SI 鍛婄煡鎺ユ敹绔紝SI 璁捐涓嶆槸鏈疄楠岄噸鐐广€?%
% 鏈疄楠屼笌 Case A fixed-channel stress test 浜掕ˉ锛?% - random channel BER 璇存槑鏅€氶殢鏈轰俊閬撲笅鐨勫钩鍧囨€ц兘锛?% - fixed Case A BER 璇存槑 GPS 鐨?worst-case vulnerability銆?% 濡傛灉 GPS 鍦ㄩ殢鏈轰俊閬撲笅涓嶅嚭鐜?BER floor锛屽苟涓嶅惁瀹?fixed stress test 鐨勬剰涔夈€?
rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
addpath(rootDir);
setup_paths(rootDir);

% ========================
% 绯荤粺涓庝豢鐪熷弬鏁?% ========================
rng(1, 'twister');
if ~exist('M', 'var'), M = 64; end
if ~exist('N', 'var'), N = M; end
if ~exist('Lcpp', 'var'), Lcpp = 8; end
if ~exist('c1', 'var'), c1 = 1 / M; end
if ~exist('modulation', 'var'), modulation = 'BPSK'; end
if ~exist('V', 'var'), V = 8; end
if ~exist('W', 'var'), W = 2; end %#ok<NASGU>
if ~exist('grouping', 'var'), grouping = 'contiguous'; end
if ~exist('delta_ratio_list', 'var'), delta_ratio_list = [0.02 0.05 0.1 0.2 0.5]; end

if ~exist('SNR_dB_list', 'var'), SNR_dB_list = 0:2:30; end
if ~exist('minErrTarget', 'var'), minErrTarget = 300; end
if ~exist('maxBits', 'var'), maxBits = 2e5; end
if ~exist('maxFrames', 'var'), maxFrames = 5000; end
if ~exist('targetSNR', 'var'), targetSNR = 24; end

% Choi-style medium random channel profile.
if ~exist('channel_profile', 'var'), channel_profile = 'medium'; end
[numPaths, lmax, alpha_max] = channel_profile_params(channel_profile, Lcpp);

% ========================
% 鏋勯€犲熀纭€閰嶇疆
% ========================
baseCfg = afdm_config();
baseCfg.waveform.NumSubcarriers = N;
baseCfg.waveform.CPPLength = Lcpp;
baseCfg.waveform.c1 = c1;
baseCfg.modulation.M_mod = 2;
baseCfg.modulation.modType = 'psk';
baseCfg.pre_chirp.base_c2 = sqrt(2) / (10 * N);
baseCfg.pre_chirp.num_groups = V;
baseCfg.pre_chirp.num_candidates = 2;
baseCfg.pre_chirp.grouping = grouping;
baseCfg.channel.multipath = true;
baseCfg.channel.add_noise = true;
baseCfg.simulation.refresh_channel_per_frame = false;

numSnr = numel(SNR_dB_list);
numDelta = numel(delta_ratio_list);

ber_base = zeros(1, numSnr);
ber_gps = zeros(1, numSnr);
ber_prop_all = zeros(numDelta, numSnr);
ber_base_plot = zeros(1, numSnr);
ber_gps_plot = zeros(1, numSnr);
ber_prop_plot_all = zeros(numDelta, numSnr);

err_base = zeros(1, numSnr);
err_gps = zeros(1, numSnr);
err_prop_all = zeros(numDelta, numSnr);
bit_base = zeros(1, numSnr);
bit_gps = zeros(1, numSnr);
bit_prop_all = zeros(numDelta, numSnr);
frames_used = zeros(1, numSnr);

fprintf('========== Choi-style random channel BER ==========\n');
fprintf('Channel profile=%s, P=%d, lmax=%d, alpha_max=%d\n', ...
    channel_profile, numPaths, lmax, alpha_max);
fprintf('N=%d, Lcpp=%d, c1=%.6g, V=%d, modulation=%s\n', N, Lcpp, c1, V, modulation);

for snrIdx = 1:numSnr
    snrDb = SNR_dB_list(snrIdx);
    fprintf('\n--- SNR %.1f dB ---\n', snrDb);

    errBase = 0;
    errGps = 0;
    errProp = zeros(numDelta, 1);
    bitBase = 0;
    bitGps = 0;
    bitProp = zeros(numDelta, 1);
    frameIdx = 0;

    while frameIdx < maxFrames && ~all_done(errBase, bitBase, errGps, bitGps, errProp, bitProp, minErrTarget, maxBits)
        frameIdx = frameIdx + 1;
        frameSeed = 20260509 + 1000000 * snrIdx + frameIdx;
        rng(frameSeed, 'twister');
        txBits = randi([0, 1], N, 1);
        ch = generate_choi_style_channel(N, numPaths, lmax, alpha_max);

        if ~(errBase >= minErrTarget || bitBase >= maxBits)
            cfg = prepare_frame_cfg(baseCfg, 'baseline', snrDb, txBits, ch);
            frame = simulate_frame(cfg, frameSeed + 10);
            errBase = errBase + frame.err_bits;
            bitBase = bitBase + frame.total_bits;
        end

        if ~(errGps >= minErrTarget || bitGps >= maxBits)
            cfg = prepare_frame_cfg(baseCfg, 'paper_grouping', snrDb, txBits, ch);
            frame = simulate_frame(cfg, frameSeed + 20);
            errGps = errGps + frame.err_bits;
            bitGps = bitGps + frame.total_bits;
        end

        for deltaIdx = 1:numDelta
            if errProp(deltaIdx) >= minErrTarget || bitProp(deltaIdx) >= maxBits
                continue;
            end
            cfg = baseCfg;
            cfg.pre_chirp.delta = delta_ratio_list(deltaIdx) * baseCfg.pre_chirp.base_c2;
            cfg = prepare_frame_cfg(cfg, 'proposed_grouping', snrDb, txBits, ch);
            frame = simulate_frame(cfg, frameSeed + 100 + deltaIdx);
            errProp(deltaIdx) = errProp(deltaIdx) + frame.err_bits;
            bitProp(deltaIdx) = bitProp(deltaIdx) + frame.total_bits;
        end

        if mod(frameIdx, 250) == 0
            fprintf('frame %4d | base %.2e GPS %.2e prop0.1 %.2e\n', ...
                frameIdx, errBase / max(bitBase, 1), errGps / max(bitGps, 1), ...
                errProp(min(3, numDelta)) / max(bitProp(min(3, numDelta)), 1));
        end
    end

    frames_used(snrIdx) = frameIdx;
    err_base(snrIdx) = errBase;
    err_gps(snrIdx) = errGps;
    err_prop_all(:, snrIdx) = errProp;
    bit_base(snrIdx) = bitBase;
    bit_gps(snrIdx) = bitGps;
    bit_prop_all(:, snrIdx) = bitProp;

    ber_base(snrIdx) = errBase / max(bitBase, 1);
    ber_gps(snrIdx) = errGps / max(bitGps, 1);
    ber_prop_all(:, snrIdx) = errProp ./ max(bitProp, 1);
    ber_base_plot(snrIdx) = plot_floor(ber_base(snrIdx), errBase, bitBase);
    ber_gps_plot(snrIdx) = plot_floor(ber_gps(snrIdx), errGps, bitGps);
    for deltaIdx = 1:numDelta
        ber_prop_plot_all(deltaIdx, snrIdx) = plot_floor(ber_prop_all(deltaIdx, snrIdx), errProp(deltaIdx), bitProp(deltaIdx));
    end

    fprintf('Baseline BER %.3e (%d/%d)\n', ber_base(snrIdx), errBase, bitBase);
    fprintf('GPS      BER %.3e (%d/%d)\n', ber_gps(snrIdx), errGps, bitGps);
    for deltaIdx = 1:numDelta
        fprintf('Proposed delta/c2=%.3f BER %.3e (%d/%d)\n', ...
            delta_ratio_list(deltaIdx), ber_prop_all(deltaIdx, snrIdx), errProp(deltaIdx), bitProp(deltaIdx));
    end
end

% ========================
% 淇濆瓨涓庣粯鍥?% ========================
outputDir = fullfile(fileparts(rootDir), 'results');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

cfg_summary = struct();
cfg_summary.M = M;
cfg_summary.N = N;
cfg_summary.Lcpp = Lcpp;
cfg_summary.c1 = c1;
cfg_summary.modulation = modulation;
cfg_summary.V = V;
cfg_summary.grouping = grouping;
cfg_summary.channel_profile = channel_profile;
cfg_summary.numPaths = numPaths;
cfg_summary.lmax = lmax;
cfg_summary.alpha_max = alpha_max;
cfg_summary.minErrTarget = minErrTarget;
cfg_summary.maxBits = maxBits;
cfg_summary.maxFrames = maxFrames;

save(fullfile(outputDir, 'results_choi_style_random_channel_ber.mat'), ...
    'SNR_dB_list', 'delta_ratio_list', ...
    'ber_base', 'ber_gps', 'ber_prop_all', ...
    'ber_base_plot', 'ber_gps_plot', 'ber_prop_plot_all', ...
    'err_base', 'err_gps', 'err_prop_all', ...
    'bit_base', 'bit_gps', 'bit_prop_all', ...
    'frames_used', 'cfg_summary');

plot_random_channel_ber(SNR_dB_list, ber_base_plot, ber_gps_plot, ber_prop_plot_all, ...
    delta_ratio_list, outputDir, timestamp);
plot_random_channel_ber_vs_delta(SNR_dB_list, ber_base_plot, ber_gps_plot, ber_prop_plot_all, ...
    delta_ratio_list, targetSNR, outputDir, timestamp);

fprintf('\nChoi-style random channel BER summary\n');
fprintf('Channel: P=%d, lmax=%d, alpha_max=%d\n', numPaths, lmax, alpha_max);
[~, idxTarget] = min(abs(SNR_dB_list - targetSNR));
fprintf('At SNR = %.1f dB:\n', SNR_dB_list(idxTarget));
fprintf('Baseline BER = %.3e\n', ber_base(idxTarget));
fprintf('GPS BER      = %.3e\n', ber_gps(idxTarget));
for deltaIdx = 1:numDelta
    fprintf('Proposed delta/c2=%.3f BER=%.3e\n', ...
        delta_ratio_list(deltaIdx), ber_prop_all(deltaIdx, idxTarget));
end

function cfg = prepare_frame_cfg(baseCfg, scheme, snrDb, txBits, ch)
    cfg = apply_pre_chirp_scheme(baseCfg, scheme);
    cfg.channel.snr_db = snrDb;
    cfg.channel.delay_taps = ch.delays(:).';
    cfg.channel.doppler_taps = ch.dopplers(:).';
    cfg.channel.doppler_freq = ch.dopplers(:).' / cfg.waveform.NumSubcarriers;
    cfg.channel.chan_coef = ch.gains(:).';
    cfg.tx.bits = txBits;
end

function ch = generate_choi_style_channel(N, numPaths, lmax, alphaMax) %#ok<INUSD>
    if numPaths <= lmax + 1
        delays = randperm(lmax + 1, numPaths) - 1;
    else
        delays = randi([0, lmax], 1, numPaths);
    end
    dopplers = randi([-alphaMax, alphaMax], 1, numPaths);
    gains = (randn(1, numPaths) + 1i * randn(1, numPaths)) / sqrt(2 * numPaths);
    ch.delays = delays;
    ch.dopplers = dopplers;
    ch.gains = gains;
    ch.numPaths = numPaths;
end

function [numPaths, lmax, alphaMax] = channel_profile_params(profile, Lcpp)
    switch lower(profile)
        case 'mild'
            numPaths = 2;
            lmax = 3;
            alphaMax = 1;
        case 'medium'
            numPaths = 4;
            lmax = Lcpp - 1;
            alphaMax = 2;
        case 'severe'
            numPaths = 6;
            lmax = Lcpp - 1;
            alphaMax = 3;
        otherwise
            error('Unknown channel_profile: %s', profile);
    end
end

function tf = all_done(errBase, bitBase, errGps, bitGps, errProp, bitProp, minErrTarget, maxBits)
    baseDone = errBase >= minErrTarget || bitBase >= maxBits;
    gpsDone = errGps >= minErrTarget || bitGps >= maxBits;
    propDone = all(errProp >= minErrTarget | bitProp >= maxBits);
    tf = baseDone && gpsDone && propDone;
end

function value = plot_floor(ber, err, bits)
    if err == 0
        value = 0.5 / max(bits, 1);
    else
        value = ber;
    end
end

function plot_random_channel_ber(SNR, berBase, berGps, berProp, deltaList, outputDir, timestamp)
    preferred = [0.05 0.1 0.2];
    selected = [];
    for val = preferred
        [d, idx] = min(abs(deltaList - val));
        if d < 1e-12
            selected(end+1) = idx; %#ok<AGROW>
        end
    end
    if isempty(selected)
        selected = unique(round(linspace(1, numel(deltaList), min(3, numel(deltaList)))));
    end

    figure('Name', 'Choi-style random channel BER', 'Color', 'w');
    semilogy(SNR, berBase, 'o-', 'LineWidth', 2); hold on;
    semilogy(SNR, berGps, 's-', 'LineWidth', 2);
    legends = {'Baseline', 'GPS'};
    for idx = selected
        semilogy(SNR, berProp(idx, :), '^-', 'LineWidth', 2);
        legends{end+1} = sprintf('Proposed \\delta/c_2=%.3g', deltaList(idx)); %#ok<AGROW>
    end
    grid on;
    xlabel('SNR (dB)');
    ylabel('BER');
    title('BER-SNR under Choi-style random doubly selective channel');
    legend(legends, 'Location', 'southwest');
    saveas(gcf, fullfile(outputDir, ['fig_choi_style_random_channel_ber_' timestamp '.png']));
end

function plot_random_channel_ber_vs_delta(SNR, berBase, berGps, berProp, deltaList, targetSNR, outputDir, timestamp)
    [~, idxTarget] = min(abs(SNR - targetSNR));
    figure('Name', 'Random channel BER vs delta', 'Color', 'w');
    semilogy(deltaList, berProp(:, idxTarget), 'o-', 'LineWidth', 2); hold on;
    yline(berBase(idxTarget), '--', 'Baseline', 'LineWidth', 1.3);
    yline(berGps(idxTarget), '--', 'GPS', 'LineWidth', 1.3);
    grid on;
    xlabel('\delta / c_2');
    ylabel(sprintf('BER at %.1f dB', SNR(idxTarget)));
    title('BER versus perturbation size under random channel');
    legend({'Proposed', 'Baseline ref', 'GPS ref'}, 'Location', 'best');
    saveas(gcf, fullfile(outputDir, ['fig_choi_style_random_channel_ber_vs_delta_' timestamp '.png']));
end
