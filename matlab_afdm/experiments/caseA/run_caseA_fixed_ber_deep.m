% RUN_CASEA_FIXED_BER_DEEP
% 瀵?Case A + theta=pi 鐨?fixed-channel BER 鍋氶珮 bit 鏁伴獙璇併€?% 鐩爣锛氭妸 BER 浣滀负涓昏瘉鎹紝纭 GPS 鐨?error floor 涓?proposed/baseline 鐨勫樊寮傘€?
rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
addpath(rootDir);
setup_paths(rootDir);

% ========================
% 鍙傛暟闆嗕腑璁剧疆
% ========================
rng(1, 'twister');
SNR_dB = 0:2:34;
maxBitsPerSNR = 1e6;        % 姣忎釜 SNR銆佹瘡涓?scheme 鐨勬渶澶?bit 鏁?minErrTarget = 300;         % 璇爜杈惧埌璇ュ€煎悗鍙彁鍓嶅仠姝?theta = pi;                 % Case A worst fixed phase
delta_ratio = 0.2;          % delta sweep 涓緝绋崇殑鎶樹腑鐐?M_mod = 2;                  % BPSK

N = 64;
V = 4;
c1 = 7 / (2 * N);
c2_base = sqrt(2) / (10 * N);
delta = delta_ratio * c2_base;
gps_pattern = [2 2 1 1];    % Case A 涓?GPS 鐨勫嵄闄?pattern
prop_pattern = [2 2 1 1];   % proposed 浣跨敤鍚?group pattern 浣滃浐瀹氬鐓?
[c2_gps, ~] = build_c2m_gps_pattern(N, V, gps_pattern);
[c2_prop, ~] = build_c2m_proposed_pattern(N, V, prop_pattern, c2_base, delta);

schemes = {'baseline', 'GPS', 'proposed'};
c2_list = {c2_base, c2_gps, c2_prop};
numSchemes = numel(schemes);
numSnr = numel(SNR_dB);

ber = zeros(numSnr, numSchemes);
ber_plot = zeros(numSnr, numSchemes);
error_bits = zeros(numSnr, numSchemes);
total_bits = zeros(numSnr, numSchemes);
frames_used = zeros(numSnr, numSchemes);

fprintf('========== Deep Case A fixed-channel BER ==========\n');
fprintf('theta/pi=%.3f, delta/c2=%.3f, maxBits/SNR=%g, minErrTarget=%d\n', ...
    theta / pi, delta_ratio, maxBitsPerSNR, minErrTarget);

for snrIdx = 1:numSnr
    for schemeIdx = 1:numSchemes
        cfg = build_caseA_fixed_cfg(N, c1, c2_list{schemeIdx}, SNR_dB(snrIdx), theta, M_mod);
        err = 0;
        bits = 0;
        frameIdx = 0;
        while bits < maxBitsPerSNR && err < minErrTarget
            frameIdx = frameIdx + 1;
            seed = 20260509 + 1000000 * snrIdx + 10000 * schemeIdx + frameIdx;
            frame = simulate_frame(cfg, seed);
            err = err + frame.err_bits;
            bits = bits + frame.total_bits;
        end

        error_bits(snrIdx, schemeIdx) = err;
        total_bits(snrIdx, schemeIdx) = bits;
        frames_used(snrIdx, schemeIdx) = frameIdx;
        ber(snrIdx, schemeIdx) = err / max(bits, 1);
        if err == 0
            ber_plot(snrIdx, schemeIdx) = 0.5 / max(bits, 1);
        else
            ber_plot(snrIdx, schemeIdx) = ber(snrIdx, schemeIdx);
        end

        fprintf('%-8s | SNR %5.1f dB | errors %6d/%-8d | BER %.3e | frames %d\n', ...
            schemes{schemeIdx}, SNR_dB(snrIdx), err, bits, ber(snrIdx, schemeIdx), frameIdx);
    end
end

div_order = estimate_diversity_order_reliable(SNR_dB, ber, error_bits);

outputDir = fullfile(fileparts(rootDir), 'results');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

figure('Name', 'Deep Case A fixed BER', 'Color', 'w');
semilogy(SNR_dB, ber_plot, 'o-', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('BER');
title(sprintf('Case A fixed-channel BER, theta/pi=%.1f, BPSK', theta / pi));
legend({'Baseline', 'GPS', 'Proposed'}, 'Location', 'southwest');
saveas(gcf, fullfile(outputDir, ['caseA_fixed_ber_deep_' timestamp '.png']));

results = struct();
results.SNR_dB = SNR_dB;
results.schemes = schemes;
results.ber = ber;
results.ber_plot = ber_plot;
results.error_bits = error_bits;
results.total_bits = total_bits;
results.frames_used = frames_used;
results.div_order = div_order;
results.theta = theta;
results.delta_ratio = delta_ratio;
results.maxBitsPerSNR = maxBitsPerSNR;
results.minErrTarget = minErrTarget;

save(fullfile(outputDir, 'caseA_fixed_ber_deep.mat'), 'results');
fprintf('Diversity estimates: baseline %.3g | GPS %.3g | proposed %.3g\n', div_order);
fprintf('Saved deep BER results to %s\n', fullfile(outputDir, 'caseA_fixed_ber_deep.mat'));

function cfg = build_caseA_fixed_cfg(N, c1, c2, snrDb, theta, MMod)
    cfg = afdm_config();
    cfg.waveform.NumSubcarriers = N;
    cfg.waveform.CPPLength = 2;
    cfg.waveform.c1 = c1;
    cfg.waveform.c2 = c2;
    cfg.modulation.M_mod = MMod;
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
    cfg.pre_chirp.profile.c2 = c2;
end

function div = estimate_diversity_order_reliable(snrDb, ber, err)
    snrLinear = 10 .^ (snrDb(:) / 10);
    div = NaN(1, size(ber, 2));
    for schemeIdx = 1:size(ber, 2)
        idx = snrDb(:) >= 16 & ber(:, schemeIdx) > 0 & err(:, schemeIdx) >= 20;
        if nnz(idx) < 3
            idx = snrDb(:) >= 12 & ber(:, schemeIdx) > 0;
        end
        if nnz(idx) >= 3
            p = polyfit(log10(snrLinear(idx)), log10(ber(idx, schemeIdx)), 1);
            div(schemeIdx) = -p(1);
        end
    end
end
