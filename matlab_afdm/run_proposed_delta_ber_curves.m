% RUN_PROPOSED_DELTA_BER_CURVES
% 比较 proposed 不同 delta/c2 在 Case A fixed-channel 下的 BER-SNR 曲线。
% 该图用于观察 delta 对 BER 鲁棒性的影响，baseline/GPS 作为参考。

rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
setup_paths(rootDir);

% ========================
% 可调参数
% ========================
rng(1, 'twister');
delta_ratio_list = [0.01 0.02 0.05 0.1 0.2 0.5 1.0];
SNR_dB = 0:2:34;
theta = pi;
maxBitsPerSNR = 2e5;
minErrTarget = 150;
M_mod = 2;

N = 64;
V = 4;
c1 = 7 / (2 * N);
c2_base = sqrt(2) / (10 * N);
gps_pattern = [2 2 1 1];
prop_pattern = [2 2 1 1];

[c2_gps, ~] = build_c2m_gps_pattern(N, V, gps_pattern);

numSnr = numel(SNR_dB);
numDelta = numel(delta_ratio_list);
numCurves = 2 + numDelta;
curve_labels = strings(1, numCurves);
curve_labels(1) = "Baseline";
curve_labels(2) = "GPS";
for idx = 1:numDelta
    curve_labels(2 + idx) = sprintf("Proposed %.2g", delta_ratio_list(idx));
end

ber = zeros(numSnr, numCurves);
ber_plot = zeros(numSnr, numCurves);
error_bits = zeros(numSnr, numCurves);
total_bits = zeros(numSnr, numCurves);

fprintf('========== Proposed delta BER curves ==========\n');
fprintf('Case A theta/pi=%.2f, maxBits/SNR=%g, minErrTarget=%d\n', ...
    theta / pi, maxBitsPerSNR, minErrTarget);

c2_curves = cell(1, numCurves);
c2_curves{1} = c2_base;
c2_curves{2} = c2_gps;
for idx = 1:numDelta
    [c2_curves{2 + idx}, ~] = build_c2m_proposed_pattern( ...
        N, V, prop_pattern, c2_base, delta_ratio_list(idx) * c2_base);
end

for curveIdx = 1:numCurves
    fprintf('\nRunning %s\n', curve_labels(curveIdx));
    for snrIdx = 1:numSnr
        [ber(snrIdx, curveIdx), error_bits(snrIdx, curveIdx), total_bits(snrIdx, curveIdx)] = ...
            run_caseA_fixed_ber_once(c2_curves{curveIdx}, SNR_dB(snrIdx), theta, c1, M_mod, ...
            maxBitsPerSNR, minErrTarget, 20260509 + 100000 * curveIdx + 1000 * snrIdx);
        if error_bits(snrIdx, curveIdx) == 0
            ber_plot(snrIdx, curveIdx) = 0.5 / max(total_bits(snrIdx, curveIdx), 1);
        else
            ber_plot(snrIdx, curveIdx) = ber(snrIdx, curveIdx);
        end
        fprintf('%-14s | SNR %5.1f dB | errors %5d/%-8d | BER %.3e\n', ...
            curve_labels(curveIdx), SNR_dB(snrIdx), error_bits(snrIdx, curveIdx), ...
            total_bits(snrIdx, curveIdx), ber(snrIdx, curveIdx));
    end
end

div_order = estimate_diversity_order_reliable(SNR_dB, ber, error_bits);

outputDir = fullfile(fileparts(rootDir), 'results');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

figure('Name', 'Proposed delta BER curves', 'Color', 'w');
semilogy(SNR_dB, ber_plot, 'o-', 'LineWidth', 1.8);
grid on;
xlabel('SNR (dB)');
ylabel('BER');
title(sprintf('Case A fixed-channel BER vs proposed delta, theta/pi=%.1f', theta / pi));
legend(cellstr(curve_labels), 'Location', 'southwest');
saveas(gcf, fullfile(outputDir, ['proposed_delta_ber_curves_' timestamp '.png']));

results = struct();
results.SNR_dB = SNR_dB;
results.delta_ratio_list = delta_ratio_list;
results.curve_labels = curve_labels;
results.ber = ber;
results.ber_plot = ber_plot;
results.error_bits = error_bits;
results.total_bits = total_bits;
results.div_order = div_order;
results.maxBitsPerSNR = maxBitsPerSNR;
results.minErrTarget = minErrTarget;
results.theta = theta;
save(fullfile(outputDir, ['proposed_delta_ber_curves_' timestamp '.mat']), 'results');

summary = table(curve_labels(:), div_order(:), ...
    'VariableNames', {'curve_label', 'div_order_est'});
writetable(summary, fullfile(outputDir, ['proposed_delta_ber_summary_' timestamp '.csv']));
disp(summary);

function [ber, err, bits] = run_caseA_fixed_ber_once(c2Vec, snrDb, theta, c1, MMod, maxBits, minErrTarget, seedBase)
    N = 64;
    cfg = afdm_config();
    cfg.waveform.NumSubcarriers = N;
    cfg.waveform.CPPLength = 2;
    cfg.waveform.c1 = c1;
    cfg.waveform.c2 = c2Vec;
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
    cfg.pre_chirp.profile.c2 = c2Vec;

    err = 0;
    bits = 0;
    frameIdx = 0;
    while bits < maxBits && err < minErrTarget
        frameIdx = frameIdx + 1;
        frame = simulate_frame(cfg, seedBase + frameIdx);
        err = err + frame.err_bits;
        bits = bits + frame.total_bits;
    end
    ber = err / max(bits, 1);
end

function div = estimate_diversity_order_reliable(snrDb, ber, err)
    snrLinear = 10 .^ (snrDb(:) / 10);
    div = NaN(1, size(ber, 2));
    for curveIdx = 1:size(ber, 2)
        idx = snrDb(:) >= 16 & ber(:, curveIdx) > 0 & err(:, curveIdx) >= 20;
        if nnz(idx) < 3
            idx = snrDb(:) >= 12 & ber(:, curveIdx) > 0;
        end
        if nnz(idx) >= 3
            p = polyfit(log10(snrLinear(idx)), log10(ber(idx, curveIdx)), 1);
            div(curveIdx) = -p(1);
        end
    end
end
