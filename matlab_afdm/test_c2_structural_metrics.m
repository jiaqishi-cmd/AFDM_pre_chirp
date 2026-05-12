% TEST_C2_STRUCTURAL_METRICS
% 完全信道无关地比较 baseline / Yuan GPS / proposed 的 c2 pattern 结构风险。

rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
setup_paths(rootDir);

% ========================
% 参数集中设置
% ========================
M = 64;
N = M;
V = 4;
c2_base = sqrt(2) / (10 * N);
delta = c2_base / 16;
gps_pattern = [2 2 1 1];        % pattern_half2，前面 Case A 搜索中较危险
proposed_pattern = [2 2 1 1];   % 1->c2-delta, 2->c2, 3->c2+delta

cfg = struct();
cfg.M = M;
cfg.N = N;
cfg.modType = 'BPSK';
cfg.phase_tol = pi / 24;
cfg.num_bins = 64;
cfg.struct_weights = [0.3, 0.2, 0.5];

outputDir = fullfile(fileparts(rootDir), 'results');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% ========================
% 构造三类 c2 pattern
% ========================
c2_vec_base = c2_base * ones(M, 1);
[c2_vec_gps, ~] = build_c2m_gps_pattern(N, V, gps_pattern);
[c2_vec_prop, ~] = build_c2m_proposed_pattern(N, V, proposed_pattern, c2_base, delta);

method_name = ["baseline"; "GPS"; "proposed"];
c2_list = {c2_vec_base, c2_vec_gps, c2_vec_prop};

metrics = cell(numel(c2_list), 1);
phase_entropy = zeros(numel(c2_list), 1);
R_phase = zeros(numel(c2_list), 1);
eff_bins = zeros(numel(c2_list), 1);
eff_phase_ratio = zeros(numel(c2_list), 1);
align_ratio = zeros(numel(c2_list), 1);
R_struct = zeros(numel(c2_list), 1);

for idx = 1:numel(c2_list)
    metrics{idx} = calc_c2_structural_metrics(c2_list{idx}, cfg);
    phase_entropy(idx) = metrics{idx}.phase_entropy;
    R_phase(idx) = metrics{idx}.phase_degeneracy_risk;
    eff_bins(idx) = metrics{idx}.eff_phase_bins;
    eff_phase_ratio(idx) = metrics{idx}.eff_phase_ratio;
    align_ratio(idx) = metrics{idx}.constellation_alignment_ratio;
    R_struct(idx) = metrics{idx}.R_struct;
end

results_table = table(method_name, phase_entropy, R_phase, eff_bins, ...
    eff_phase_ratio, align_ratio, R_struct);
disp(results_table);

% ========================
% 候选相关性示例：proposed 三个 group-wise 常量 candidate
% ========================
candidate_vecs = { ...
    (c2_base - delta) * ones(M, 1), ...
    c2_base * ones(M, 1), ...
    (c2_base + delta) * ones(M, 1)};
corr_metrics = calc_c2_candidate_correlation(candidate_vecs, cfg); %#ok<NASGU>
fprintf('Proposed candidate max rho offdiag = %.4f, min separation = %.4f\n', ...
    corr_metrics.max_rho_offdiag, corr_metrics.min_separation);

% ========================
% 图 1：phase histogram
% ========================
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
figure('Name', 'c2 structural phase histograms', 'Color', 'w');
for idx = 1:numel(c2_list)
    subplot(numel(c2_list), 1, idx);
    edges = metrics{idx}.phase_hist_edges / pi;
    counts = metrics{idx}.phase_hist_counts;
    centers = 0.5 * (edges(1:end-1) + edges(2:end));
    bar(centers, counts, 1.0);
    grid on;
    xlim([0, 2]);
    ylabel('Norm. count');
    title(sprintf('%s phase histogram', method_name(idx)));
    if idx == numel(c2_list)
        xlabel('phase / \pi');
    end
end
saveas(gcf, fullfile(outputDir, ['c2_structural_phase_hist_' timestamp '.png']));

% ========================
% 图 2：risk bar 对比
% ========================
figure('Name', 'c2 structural risk comparison', 'Color', 'w');
bar([R_struct, R_phase, align_ratio]);
grid on;
xticklabels(method_name);
ylabel('Risk / ratio');
title('Channel-independent c2 structural risk indicators');
legend({'R\_struct', 'R\_phase', 'alignment ratio'}, 'Location', 'northwest');
saveas(gcf, fullfile(outputDir, ['c2_structural_risk_bar_' timestamp '.png']));

save(fullfile(outputDir, 'results_c2_structural_metrics.mat'), ...
    'results_table', 'metrics', 'corr_metrics', 'cfg', ...
    'c2_vec_base', 'c2_vec_gps', 'c2_vec_prop', ...
    'gps_pattern', 'proposed_pattern');
writetable(results_table, fullfile(outputDir, 'results_c2_structural_metrics.csv'));

fprintf('Note: R_struct is channel-independent and should be validated by key-equation mismatch or BER theta scan.\n');
