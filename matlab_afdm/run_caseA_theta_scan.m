% RUN_CASEA_THETA_SCAN
% 扫描 Case A 固定两径信道的路径增益相位 theta，验证 GPS 在特定
% fixed-channel phase 下的 worst-case diversity vulnerability。
%
% 定位说明：
% 本实验不是为了证明 GPS 在 Rayleigh 平均信道下一定失败，而是验证：
% 在 Case A 这种 Bemani key-equation 危险 delay-Doppler 条件下，
% GPS 的特定 pattern 可能对某些路径相位出现明显 BER floor。

rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
setup_paths(rootDir);

% ========================
% 可调仿真参数
% ========================
rng(1, 'twister');
if ~exist('theta_list', 'var')
    theta_list = linspace(0, 2 * pi, 181);   % 0 到 2pi，每 1 度一个点
end
if ~exist('SNR_dB', 'var')
    SNR_dB = 28;                             % 高 SNR 下观察 BER floor
end
if ~exist('minErrTarget', 'var')
    minErrTarget = 100;                      % 每个 theta/scheme 至少累计的误比特目标
end
if ~exist('maxBits', 'var')
    maxBits = 2e5;                           % 每个 theta/scheme 的最大仿真 bit 数
end
if ~exist('maxFrames', 'var')
    maxFrames = ceil(maxBits / 64) + 10;     % BPSK 每帧约 N=64 bit
end

% ========================
% 固定 Case A 参数
% ========================
N = 64;
V = 4;
alpha_max = 3;
c1 = (2 * alpha_max + 1) / (2 * N);
c2_base = sqrt(2) / (10 * N);
gps_pattern = [2 2 1 1];                    % pattern_half2，前面搜索出的危险 pattern
proposed_pattern = gps_pattern;
proposed_delta = c2_base / 16;

l1 = 0;
alpha1 = 0;
l2 = 2;
alpha2 = 2;
L = 16; %#ok<NASGU>

[c2_gps, ~] = build_c2m_gps_pattern(N, V, gps_pattern);
[c2_prop, ~] = build_c2m_proposed_pattern(N, V, proposed_pattern, c2_base, proposed_delta);

schemes = {'baseline', 'GPS', 'proposed'};
c2_values = {c2_base, c2_gps, c2_prop};
numTheta = numel(theta_list);
numSchemes = numel(schemes);

ber = zeros(numTheta, numSchemes);
ber_plot = zeros(numTheta, numSchemes);
error_count = zeros(numTheta, numSchemes);
total_bits = zeros(numTheta, numSchemes);
frames_used = zeros(numTheta, numSchemes);

outputDir = fullfile(fileparts(rootDir), 'results');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
progressPath = fullfile(outputDir, 'results_caseA_theta_scan_progress.mat');

fprintf('========== Case A theta-scan BER ==========\n');
fprintf('N=%d, V=%d, Case A path2=(l=%d, alpha=%d), SNR=%g dB\n', N, V, l2, alpha2, SNR_dB);
fprintf('theta points=%d, minErrTarget=%d, maxBits=%g\n', numTheta, minErrTarget, maxBits);

for thetaIdx = 1:numTheta
    theta = theta_list(thetaIdx);
    h = [1, exp(1i * theta)] / sqrt(2);

    for schemeIdx = 1:numSchemes
        cfg = build_caseA_config(N, c1, c2_values{schemeIdx}, SNR_dB, ...
            [l1, l2], [alpha1, alpha2], h);

        err = 0;
        bits = 0;
        frameIdx = 0;
        while err < minErrTarget && bits < maxBits && frameIdx < maxFrames
            frameIdx = frameIdx + 1;
            seed = 1 + 1000000 * thetaIdx + 10000 * schemeIdx + frameIdx;
            frame = simulate_frame(cfg, seed);
            err = err + frame.err_bits;
            bits = bits + frame.total_bits;
        end

        error_count(thetaIdx, schemeIdx) = err;
        total_bits(thetaIdx, schemeIdx) = bits;
        frames_used(thetaIdx, schemeIdx) = frameIdx;
        ber(thetaIdx, schemeIdx) = err / max(bits, 1);

        % 只用于 semilogy 显示；原始 BER 仍然保存在 ber 中。
        if err == 0
            ber_plot(thetaIdx, schemeIdx) = 0.5 / max(bits, 1);
        else
            ber_plot(thetaIdx, schemeIdx) = ber(thetaIdx, schemeIdx);
        end
    end

    fprintf('theta/pi=%6.3f | base %.3e (%d/%d) | GPS %.3e (%d/%d) | prop %.3e (%d/%d)\n', ...
        theta / pi, ...
        ber(thetaIdx, 1), error_count(thetaIdx, 1), total_bits(thetaIdx, 1), ...
        ber(thetaIdx, 2), error_count(thetaIdx, 2), total_bits(thetaIdx, 2), ...
        ber(thetaIdx, 3), error_count(thetaIdx, 3), total_bits(thetaIdx, 3));

    % 逐 theta 保存进度，避免长扫描被中断后完全丢失结果。
    save(progressPath, ...
        'theta_list', 'SNR_dB', 'schemes', 'ber', 'ber_plot', ...
        'error_count', 'total_bits', 'frames_used', ...
        'minErrTarget', 'maxBits', 'gps_pattern', 'thetaIdx');
end

ber_base = ber(:, 1);
ber_gps = ber(:, 2);
ber_prop = ber(:, 3);
ber_base_plot = ber_plot(:, 1);
ber_gps_plot = ber_plot(:, 2);
ber_prop_plot = ber_plot(:, 3);

save(fullfile(outputDir, 'results_caseA_theta_scan.mat'), ...
    'theta_list', 'SNR_dB', 'schemes', 'ber', 'ber_plot', ...
    'ber_base', 'ber_gps', 'ber_prop', ...
    'ber_base_plot', 'ber_gps_plot', 'ber_prop_plot', ...
    'error_count', 'total_bits', 'frames_used', ...
    'minErrTarget', 'maxBits', 'gps_pattern');

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
figure('Name', 'Case A theta scan BER', 'Color', 'w');
semilogy(theta_list / pi, ber_base_plot, 'LineWidth', 2); hold on;
semilogy(theta_list / pi, ber_gps_plot, 'LineWidth', 2);
semilogy(theta_list / pi, ber_prop_plot, 'LineWidth', 2);
xline(1, '--k', '\theta/\pi=1', 'LineWidth', 1.2);
grid on;
xlim([0, 2]);
xlabel('\theta / \pi');
ylabel('BER');
title(sprintf('BER vs path phase theta on Case A fixed two-path channel, SNR=%g dB', SNR_dB));
legend({'baseline', 'GPS', 'proposed'}, 'Location', 'best');
saveas(gcf, fullfile(outputDir, ['caseA_theta_scan_BER_' timestamp '.png']));

% TODO: add minimum-distance or key-equation residual scan here if needed.

[worstGpsBer, worstIdx] = max(ber_gps);
[~, piIdx] = min(abs(theta_list - pi));
fprintf('\nWorst GPS theta/pi = %.3f, BER = %.3e\n', theta_list(worstIdx) / pi, worstGpsBer);
fprintf('At theta=pi: baseline=%.3e, GPS=%.3e, proposed=%.3e\n', ...
    ber_base(piIdx), ber_gps(piIdx), ber_prop(piIdx));
fprintf('Saved theta scan to %s\n', fullfile(outputDir, 'results_caseA_theta_scan.mat'));

function cfg = build_caseA_config(N, c1, c2, snrDb, delayTaps, dopplerTaps, gains)
    % 构造 Case A 的固定信道配置。接收端仍使用现有 MMSE 链路和理想 CSI。
    cfg = afdm_config();
    cfg.waveform.NumSubcarriers = N;
    cfg.waveform.CPPLength = max(delayTaps);
    cfg.waveform.c1 = c1;
    cfg.waveform.c2 = c2;
    cfg.modulation.M_mod = 2;
    cfg.modulation.modType = 'psk';
    cfg.channel.multipath = true;
    cfg.channel.add_noise = true;
    cfg.channel.snr_db = snrDb;
    cfg.channel.delay_taps = delayTaps;
    cfg.channel.doppler_taps = dopplerTaps;
    cfg.channel.doppler_freq = dopplerTaps / N;
    cfg.channel.chan_coef = gains;
    cfg.pre_chirp.scheme = 'baseline';
    cfg.pre_chirp.profile.scheme = 'baseline';
    cfg.pre_chirp.profile.c2 = c2;
end
