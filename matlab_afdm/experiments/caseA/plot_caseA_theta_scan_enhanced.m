% PLOT_CASEA_THETA_SCAN_ENHANCED
% Case A theta scan 鐨勫寮哄悗澶勭悊缁樺浘鑴氭湰銆?% 涓嶆敼鍙樺師濮嬩豢鐪熸暟鎹紝鍙仛缁樺浘涓嬮檺銆乴og-domain 骞虫粦鍜岀粺璁″睍绀恒€?
rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
addpath(rootDir);

% ========================
% 鍙皟鍙傛暟
% ========================
smooth_win = 7;          % log10(BER) 涓婄殑绉诲姩涓€肩獥鍙?default_ber_floor = 1e-5;
ber_th = 1e-3;
ratio_th1 = 10;
ratio_th2 = 100;

outputDir = fullfile(fileparts(rootDir), 'results');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% 濡傛灉宸ヤ綔鍖烘病鏈夊彉閲忥紝鍒欓粯璁よ鍙栨渶杩戜竴娆?theta scan 缁撴灉銆?if ~exist('theta_list', 'var') || ~exist('ber_base', 'var') || ...
        ~exist('ber_gps', 'var') || ~exist('ber_prop', 'var')
    dataPath = fullfile(outputDir, 'results_caseA_theta_scan.mat');
    if ~exist(dataPath, 'file')
        error('Cannot find %s. Run run_caseA_theta_scan first.', dataPath);
    end
    loaded = load(dataPath);
    theta_list = loaded.theta_list;
    ber_base = loaded.ber_base;
    ber_gps = loaded.ber_gps;
    ber_prop = loaded.ber_prop;
    if isfield(loaded, 'SNR_dB')
        SNR_dB = loaded.SNR_dB;
    else
        SNR_dB = NaN;
    end
    if isfield(loaded, 'maxBits')
        maxBits = loaded.maxBits;
    end
end

% ========================
% 鍩烘湰棰勫鐞?% ========================
theta = theta_list(:);
b0 = ber_base(:);
b1 = ber_gps(:);
b2 = ber_prop(:);

[theta, order] = sort(theta, 'ascend');
b0 = b0(order);
b1 = b1(order);
b2 = b2(order);
theta_pi = theta / pi;

if exist('maxBits', 'var') && ~isempty(maxBits)
    ber_floor = max(default_ber_floor, 0.5 / maxBits);
else
    ber_floor = default_ber_floor;
end

b0_plot = max(b0, ber_floor);
b1_plot = max(b1, ber_floor);
b2_plot = max(b2, ber_floor);

b0_smooth = smooth_log_ber(b0_plot, smooth_win);
b1_smooth = smooth_log_ber(b1_plot, smooth_win);
b2_smooth = smooth_log_ber(b2_plot, smooth_win);

eps_ratio = ber_floor;
ratio_gps = (b1_plot + eps_ratio) ./ (b0_plot + eps_ratio);
ratio_prop = (b2_plot + eps_ratio) ./ (b0_plot + eps_ratio);
ratio_gps_smooth = smooth_log_ber(ratio_gps, smooth_win);
ratio_prop_smooth = smooth_log_ber(ratio_prop, smooth_win);

[~, idx_pi] = min(abs(theta_pi - 1));
gps_high = b1 > ber_th;

timestamp = datestr(now, 'yyyymmdd_HHMMSS');

% ========================
% 鍥?1锛氬師濮?BER + 骞虫粦 BER
% ========================
figure('Name', 'Case A theta scan raw and smoothed BER', 'Color', 'w');
semilogy(theta_pi, b0_plot, '.', 'Color', [0.65 0.78 1.00], 'MarkerSize', 7, 'HandleVisibility', 'off'); hold on;
semilogy(theta_pi, b1_plot, '.', 'Color', [1.00 0.78 0.55], 'MarkerSize', 7, 'HandleVisibility', 'off');
semilogy(theta_pi, b2_plot, '.', 'Color', [0.95 0.88 0.45], 'MarkerSize', 7, 'HandleVisibility', 'off');
semilogy(theta_pi(gps_high), b1_plot(gps_high), 'o', 'Color', [0.85 0.33 0.10], ...
    'MarkerSize', 4, 'LineWidth', 1.0, 'DisplayName', 'GPS high-risk raw');
semilogy(theta_pi, b0_smooth, '-', 'Color', [0.00 0.45 0.74], 'LineWidth', 2.4, 'DisplayName', 'Baseline smooth');
semilogy(theta_pi, b1_smooth, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 2.4, 'DisplayName', 'GPS smooth');
semilogy(theta_pi, b2_smooth, '-', 'Color', [0.93 0.69 0.13], 'LineWidth', 2.4, 'DisplayName', 'Proposed smooth');
xline(1, 'k--', '\theta/\pi=1', 'LineWidth', 1.2, 'HandleVisibility', 'off');
grid on;
xlim([0, 2]);
ylim([ber_floor, max(1e-2, min(1, 2 * max([b0_plot; b1_plot; b2_plot])))]);
xlabel('\theta / \pi');
ylabel('BER');
title(sprintf('BER vs path phase theta on Case A, raw and smoothed, SNR = %.1f dB', SNR_dB));
legend('Location', 'southwest');
text(0.02, 1.25 * ber_floor, 'Zero-error points clipped to plotting floor; smoothing is on log10(BER).', ...
    'FontSize', 9, 'Color', [0.25 0.25 0.25]);
saveas(gcf, fullfile(outputDir, ['fig_caseA_theta_raw_smooth_' timestamp '.png']));

% ========================
% 鍥?2锛氱浉瀵?baseline 鐨?BER ratio
% ========================
figure('Name', 'Case A theta scan BER ratio', 'Color', 'w');
semilogy(theta_pi, ratio_gps, '.', 'Color', [1.00 0.78 0.55], 'MarkerSize', 7, 'HandleVisibility', 'off'); hold on;
semilogy(theta_pi, ratio_prop, '.', 'Color', [0.95 0.88 0.45], 'MarkerSize', 7, 'HandleVisibility', 'off');
semilogy(theta_pi, ratio_gps_smooth, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 2.4, 'DisplayName', 'GPS smooth');
semilogy(theta_pi, ratio_prop_smooth, '-', 'Color', [0.93 0.69 0.13], 'LineWidth', 2.4, 'DisplayName', 'Proposed smooth');
yline(1, 'k-', 'LineWidth', 1.0, 'HandleVisibility', 'off');
yline(10, 'k--', '10x', 'LineWidth', 1.0, 'HandleVisibility', 'off');
yline(100, 'k:', '100x', 'LineWidth', 1.0, 'HandleVisibility', 'off');
xline(1, 'k--', '\theta/\pi=1', 'LineWidth', 1.2, 'HandleVisibility', 'off');
grid on;
xlim([0, 2]);
xlabel('\theta / \pi');
ylabel('BER ratio to baseline');
title('Relative BER sensitivity versus path phase');
legend('Location', 'northwest');
saveas(gcf, fullfile(outputDir, ['fig_caseA_theta_ratio_' timestamp '.png']));

% ========================
% 鍥?3锛氬嵄闄╃浉浣嶅尯鍩熺粺璁?% ========================
max_ber_base = max(b0);
max_ber_gps = max(b1);
max_ber_prop = max(b2);
median_ber_base = median(b0);
median_ber_gps = median(b1);
median_ber_prop = median(b2);
frac_high_base = mean(b0 > ber_th);
frac_high_gps = mean(b1 > ber_th);
frac_high_prop = mean(b2 > ber_th);
frac_ratio10_gps = mean(ratio_gps > ratio_th1);
frac_ratio10_prop = mean(ratio_prop > ratio_th1);
frac_ratio100_gps = mean(ratio_gps > ratio_th2);
frac_ratio100_prop = mean(ratio_prop > ratio_th2);

figure('Name', 'Case A theta scan statistics', 'Color', 'w');
tiledlayout(1, 2, 'TileSpacing', 'compact');
nexttile;
bar(100 * [frac_high_base, frac_high_gps, frac_high_prop]);
grid on;
xticklabels({'Baseline', 'GPS', 'Proposed'});
ylabel('Fraction of theta (%)');
title(sprintf('BER > %.1e', ber_th));

nexttile;
bar(100 * [frac_ratio10_gps, frac_ratio10_prop; frac_ratio100_gps, frac_ratio100_prop]);
grid on;
xticklabels({'Ratio > 10', 'Ratio > 100'});
ylabel('Fraction of theta (%)');
title('Ratio-to-baseline high-risk fractions');
legend({'GPS', 'Proposed'}, 'Location', 'northwest');
sgtitle('High-risk theta fraction statistics');
saveas(gcf, fullfile(outputDir, ['fig_caseA_theta_statistics_' timestamp '.png']));

save(fullfile(outputDir, 'results_caseA_theta_scan_enhanced.mat'), ...
    'theta', 'theta_pi', 'b0', 'b1', 'b2', ...
    'b0_plot', 'b1_plot', 'b2_plot', ...
    'b0_smooth', 'b1_smooth', 'b2_smooth', ...
    'ratio_gps', 'ratio_prop', 'ratio_gps_smooth', 'ratio_prop_smooth', ...
    'smooth_win', 'ber_floor', 'ber_th', ...
    'frac_high_base', 'frac_high_gps', 'frac_high_prop', ...
    'frac_ratio10_gps', 'frac_ratio10_prop', ...
    'frac_ratio100_gps', 'frac_ratio100_prop');

fprintf('Theta scan summary, SNR = %.1f dB\n', SNR_dB);
fprintf('Max BER: baseline=%.3e, GPS=%.3e, proposed=%.3e\n', max_ber_base, max_ber_gps, max_ber_prop);
fprintf('Median BER: baseline=%.3e, GPS=%.3e, proposed=%.3e\n', median_ber_base, median_ber_gps, median_ber_prop);
fprintf('Fraction BER > %.1e: baseline=%.2f%%, GPS=%.2f%%, proposed=%.2f%%\n', ...
    ber_th, 100 * frac_high_base, 100 * frac_high_gps, 100 * frac_high_prop);
fprintf('Fraction ratio > 10: GPS=%.2f%%, proposed=%.2f%%\n', ...
    100 * frac_ratio10_gps, 100 * frac_ratio10_prop);
fprintf('Fraction ratio > 100: GPS=%.2f%%, proposed=%.2f%%\n', ...
    100 * frac_ratio100_gps, 100 * frac_ratio100_prop);
fprintf('At theta/pi=1: baseline=%.3e, GPS=%.3e, proposed=%.3e\n', ...
    b0(idx_pi), b1(idx_pi), b2(idx_pi));

function y = smooth_log_ber(x, win)
    logx = log10(max(x(:), realmin));
    if exist('movmedian', 'file') == 2
        y = 10 .^ movmedian(logx, win);
    else
        y = 10 .^ movmean(logx, win);
    end
end
