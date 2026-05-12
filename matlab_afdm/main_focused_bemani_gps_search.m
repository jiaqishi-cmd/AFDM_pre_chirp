% MAIN_FOCUSED_BEMANI_GPS_SEARCH
% 聚焦 Bemani 关键等式的 GPS-AFDM 反例搜索。
% 本脚本只搜索理论上最危险的 A-F 六个两径 case，不扫描 total H_eff rank。
% 核心指标是 key-equation mismatch 和 Phi(delta)=[H1*delta,H2*delta]。

rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
setup_paths(rootDir);

N = 64;
V = 4;
W = 2; %#ok<NASGU>
M = N / V;
alpha_max = 3;
c1 = (2 * alpha_max + 1) / (2 * N);
c2_base = sqrt(2) / (10 * N);

cases = build_focused_cases();
patterns = { ...
    'pattern_alt1', [1 2 1 2]; ...
    'pattern_alt2', [2 1 2 1]; ...
    'pattern_half1', [1 1 2 2]; ...
    'pattern_half2', [2 2 1 1]; ...
    'pattern_all1', [1 1 1 1] ...
    };

fprintf('========== Focused Bemani/GPS search ==========\n');
fprintf('N=%d, V=%d, M=%d, alpha_max=%d, c1=%.10g\n', N, V, M, alpha_max, c1);
fprintf('Cases: A-F, GPS patterns: %d\n\n', size(patterns, 1));

for case_idx = 1:numel(cases)
    fprintf('Case %s: l2=%d alpha2=%d L=%d path_phase=%+.4f%+.4fi\n', ...
        cases(case_idx).name, cases(case_idx).l2, cases(case_idx).alpha2, cases(case_idx).L, ...
        real(cases(case_idx).path_phase), imag(cases(case_idx).path_phase));
end
fprintf('\n');

rows = struct([]);
row_idx = 0;

for case_idx = 1:numel(cases)
    case_def = cases(case_idx);
    H1_base = build_H_path_general_c2m(N, c1, c2_base, 0, 0);
    H2_base = build_H_path_general_c2m(N, c1, c2_base, case_def.l2, case_def.alpha2);

    for pattern_idx = 1:size(patterns, 1)
        pattern_name = patterns{pattern_idx, 1};
        pattern = patterns{pattern_idx, 2};
        [c2_gps, d_gps] = build_c2m_gps_pattern(N, V, pattern);
        H1_gps = build_H_path_general_c2m(N, c1, c2_gps, 0, 0);
        H2_gps = build_H_path_general_c2m(N, c1, c2_gps, case_def.l2, case_def.alpha2);

        delta_set = build_structured_delta_set(N, M);
        trial_rows = evaluate_delta_set(delta_set, case_def, pattern_name, d_gps, ...
            H1_base, H2_base, H1_gps, H2_gps, c2_base, N);

        % 如果结构化 delta 没有触发强信号，再追加小规模随机 delta。
        if ~has_strong_case(trial_rows)
            random_set = build_random_delta_set(N, 200, 200, 20260508 + 100 * case_idx + pattern_idx);
            trial_rows = [trial_rows, evaluate_delta_set(random_set, case_def, pattern_name, d_gps, ...
                H1_base, H2_base, H1_gps, H2_gps, c2_base, N)]; %#ok<AGROW>
        end

        for idx = 1:numel(trial_rows)
            row_idx = row_idx + 1;
            rows(row_idx).case_name = string(trial_rows(idx).case_name); %#ok<SAGROW>
            rows(row_idx).l2 = trial_rows(idx).l2;
            rows(row_idx).alpha2 = trial_rows(idx).alpha2;
            rows(row_idx).L = trial_rows(idx).L;
            rows(row_idx).path_phase = trial_rows(idx).path_phase;
            rows(row_idx).gps_pattern_name = string(trial_rows(idx).gps_pattern_name);
            rows(row_idx).delta_type = string(trial_rows(idx).delta_type);
            rows(row_idx).min_E_base = trial_rows(idx).min_E_base;
            rows(row_idx).mean_E_base = trial_rows(idx).mean_E_base;
            rows(row_idx).min_E_GPS = trial_rows(idx).min_E_GPS;
            rows(row_idx).mean_E_GPS = trial_rows(idx).mean_E_GPS;
            rows(row_idx).improvement_min = trial_rows(idx).improvement_min;
            rows(row_idx).improvement_mean = trial_rows(idx).improvement_mean;
            rows(row_idx).rank_base = trial_rows(idx).rank_base;
            rows(row_idx).rank_GPS = trial_rows(idx).rank_GPS;
            rows(row_idx).sigma_ratio_base = trial_rows(idx).sigma_ratio_base;
            rows(row_idx).sigma_ratio_GPS = trial_rows(idx).sigma_ratio_GPS;
            rows(row_idx).col_corr_base = trial_rows(idx).col_corr_base;
            rows(row_idx).col_corr_GPS = trial_rows(idx).col_corr_GPS;
        end
    end
end

results_table = struct2table(rows);
results_table.sigma_ratio_gps_over_base = results_table.sigma_ratio_GPS ./ max(results_table.sigma_ratio_base, eps);
results_table.col_corr_gain = results_table.col_corr_GPS - results_table.col_corr_base;

outputDir = fullfile(fileparts(rootDir), 'results');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
save(fullfile(outputDir, 'focused_bemani_gps_search.mat'), 'results_table');
writetable(results_table, fullfile(outputDir, 'focused_bemani_gps_search.csv'));

focused_best_cases = select_focused_best_cases(results_table, N, V, M, alpha_max, c1, c2_base); %#ok<NASGU>
save(fullfile(outputDir, 'focused_bemani_gps_best_cases.mat'), 'focused_best_cases');

print_top_tables(results_table);
plot_focused_results(results_table, outputDir);

fprintf('\nSaved focused search to %s\n', outputDir);

function cases = build_focused_cases()
    specs = { ...
        'A', 2, 2, 16; ...
        'B', 5, -3, 32; ...
        'C', 2, 1, 15; ...
        'D', 2, 3, 17; ...
        'E', 4, 3, 31; ...
        'F', 5, -2, 33 ...
        };
    N = 64;
    cases = struct('name', {}, 'l2', {}, 'alpha2', {}, 'L', {}, 'path_phase', {});
    for idx = 1:size(specs, 1)
        case_def.name = specs{idx, 1};
        case_def.l2 = specs{idx, 2};
        case_def.alpha2 = specs{idx, 3};
        case_def.L = specs{idx, 4};
        case_def.path_phase = exp(1i * 2 * pi / N * case_def.l2 * case_def.L);
        cases(end+1) = case_def; %#ok<AGROW>
    end
end

function delta_set = build_structured_delta_set(N, M)
    m = (0:N-1).';
    delta_set = struct('name', {}, 'delta', {});
    delta_set(end+1) = make_delta('delta_all1', ones(N, 1));
    delta_set(end+1) = make_delta('delta_group_alt', group_delta(N, M, [1, -1, 1, -1]));
    delta_set(end+1) = make_delta('delta_half', group_delta(N, M, [1, 1, -1, -1]));
    delta_set(end+1) = make_delta('delta_qpsk1', group_delta(N, M, [1, 1i, -1, -1i]));
    delta_set(end+1) = make_delta('delta_qpsk2', group_delta(N, M, [1, -1i, -1, 1i]));
    delta_set(end+1) = make_delta('delta_alt', (-1) .^ m);
end

function delta_set = build_random_delta_set(N, numBpsk, numQpsk, seed)
    oldRngState = rng;
    cleanup = onCleanup(@() rng(oldRngState));
    rng(seed, 'twister');

    delta_set = struct('name', {}, 'delta', {});
    for idx = 1:numBpsk
        delta_set(end+1) = make_delta(sprintf('random_bpsk_%03d', idx), 2 * randi([0, 1], N, 1) - 1); %#ok<AGROW>
    end
    alphabet = [1, -1, 1i, -1i].';
    for idx = 1:numQpsk
        delta_set(end+1) = make_delta(sprintf('random_qpsk_%03d', idx), alphabet(randi(numel(alphabet), N, 1))); %#ok<AGROW>
    end
end

function row_set = evaluate_delta_set(delta_set, case_def, pattern_name, d_gps, H1_base, H2_base, H1_gps, H2_gps, c2_base, N)
    c2_values.base = c2_base;
    row_template = struct( ...
        'case_name', '', ...
        'l2', 0, ...
        'alpha2', 0, ...
        'L', 0, ...
        'path_phase', 0, ...
        'gps_pattern_name', '', ...
        'delta_type', '', ...
        'min_E_base', 0, ...
        'mean_E_base', 0, ...
        'min_E_GPS', 0, ...
        'mean_E_GPS', 0, ...
        'improvement_min', 0, ...
        'improvement_mean', 0, ...
        'rank_base', 0, ...
        'rank_GPS', 0, ...
        'sigma_ratio_base', 0, ...
        'sigma_ratio_GPS', 0, ...
        'col_corr_base', 0, ...
        'col_corr_GPS', 0);
    row_set = repmat(row_template, 1, numel(delta_set));
    for delta_idx = 1:numel(delta_set)
        row = row_template;
        delta = delta_set(delta_idx).delta;
        metrics = compute_bemani_equation_error(delta, d_gps, case_def.L, case_def.l2, N, c2_values);
        base_phi = evaluate_phi_metrics(H1_base, H2_base, delta);
        gps_phi = evaluate_phi_metrics(H1_gps, H2_gps, delta);

        row.case_name = case_def.name;
        row.l2 = case_def.l2;
        row.alpha2 = case_def.alpha2;
        row.L = case_def.L;
        row.path_phase = case_def.path_phase;
        row.gps_pattern_name = pattern_name;
        row.delta_type = delta_set(delta_idx).name;
        row.min_E_base = metrics.min_E_base;
        row.mean_E_base = metrics.mean_E_base;
        row.min_E_GPS = metrics.min_E_GPS;
        row.mean_E_GPS = metrics.mean_E_GPS;
        row.improvement_min = metrics.min_E_base / (metrics.min_E_GPS + eps);
        row.improvement_mean = metrics.mean_E_base / (metrics.mean_E_GPS + eps);
        row.rank_base = base_phi.rank;
        row.rank_GPS = gps_phi.rank;
        row.sigma_ratio_base = base_phi.sigma_ratio;
        row.sigma_ratio_GPS = gps_phi.sigma_ratio;
        row.col_corr_base = base_phi.col_corr;
        row.col_corr_GPS = gps_phi.col_corr;
        row_set(delta_idx) = row;
    end
end

function tf = has_strong_case(rows)
    if isempty(rows)
        tf = false;
        return;
    end
    sigma_base = [rows.sigma_ratio_base];
    sigma_gps = [rows.sigma_ratio_GPS];
    corr_base = [rows.col_corr_base];
    corr_gps = [rows.col_corr_GPS];
    rank_base = [rows.rank_base];
    rank_gps = [rows.rank_GPS];
    mean_base = [rows.mean_E_base];
    mean_gps = [rows.mean_E_GPS];
    tf = any(rank_gps < 2 & rank_base == 2) || ...
        any(sigma_gps < 0.1 & sigma_base > 0.7) || ...
        any(corr_gps > 0.95 & corr_base < 0.5) || ...
        any(mean_base ./ (mean_gps + eps) > 10);
end

function item = make_delta(name, delta)
    item.name = name;
    item.delta = delta(:);
end

function delta = group_delta(N, M, values)
    delta = zeros(N, 1);
    for group_idx = 1:numel(values)
        delta((group_idx - 1) * M + 1:group_idx * M) = values(group_idx);
    end
end

function print_top_tables(T)
    print_sorted(T, 'mean_E_GPS', 'ascend', 'Top 20 by mean_E_GPS');
    print_sorted(T, 'improvement_mean', 'descend', 'Top 20 by improvement_mean');
    print_sorted(T, 'sigma_ratio_gps_over_base', 'ascend', 'Top 20 by sigma_ratio_GPS/sigma_ratio_base');
    print_sorted(T, 'col_corr_gain', 'descend', 'Top 20 by col_corr_GPS-col_corr_base');
end

function focused_best_cases = select_focused_best_cases(T, N, V, M, alpha_max, c1, c2_base)
    % 保存两个代表性 case：
    % 1) mismatch_case：key-equation mismatch 改善最大；
    % 2) phi_case：Phi(delta) 的 sigma_ratio 降低、列相关升高最明显。
    [~, mismatch_idx] = max(T.improvement_mean);
    [~, phi_idx] = min(T.sigma_ratio_gps_over_base);

    focused_best_cases.mismatch_case = table_row_to_best_case(T(mismatch_idx, :), N, V, M, alpha_max, c1, c2_base);
    focused_best_cases.phi_case = table_row_to_best_case(T(phi_idx, :), N, V, M, alpha_max, c1, c2_base);

    fprintf('\nFocused mismatch best case: Case %s L=%d pattern=%s delta=%s improvement=%.4g\n', ...
        focused_best_cases.mismatch_case.case_name, focused_best_cases.mismatch_case.L, ...
        focused_best_cases.mismatch_case.gps_pattern, focused_best_cases.mismatch_case.delta_type, ...
        focused_best_cases.mismatch_case.improvement_mean);
    fprintf('Focused Phi best case: Case %s L=%d pattern=%s delta=%s sigma GPS/base=%.4g corr gain=%.4g\n', ...
        focused_best_cases.phi_case.case_name, focused_best_cases.phi_case.L, ...
        focused_best_cases.phi_case.gps_pattern, focused_best_cases.phi_case.delta_type, ...
        focused_best_cases.phi_case.sigma_ratio_GPS / max(focused_best_cases.phi_case.sigma_ratio_base, eps), ...
        focused_best_cases.phi_case.col_corr_GPS - focused_best_cases.phi_case.col_corr_base);
end

function best_case = table_row_to_best_case(row, N, V, M, alpha_max, c1, c2_base)
    % run_bestcase_ber_snr 需要的字段在这里集中补齐。
    best_case = struct();
    best_case.N = N;
    best_case.V = V;
    best_case.M = M;
    best_case.alpha_max = alpha_max;
    best_case.c1 = c1;
    best_case.c2_base = c2_base;
    best_case.l1 = 0;
    best_case.alpha1 = 0;
    best_case.l2 = row.l2;
    best_case.alpha2 = row.alpha2;
    best_case.L = row.L;
    best_case.case_name = char(row.case_name);
    best_case.gps_pattern = char(row.gps_pattern_name);
    best_case.proposed_pattern = char(row.gps_pattern_name);
    best_case.delta_type = char(row.delta_type);
    best_case.mean_E_base = row.mean_E_base;
    best_case.mean_E_GPS = row.mean_E_GPS;
    best_case.improvement_mean = row.improvement_mean;
    best_case.rank_base = row.rank_base;
    best_case.rank_GPS = row.rank_GPS;
    best_case.sigma_ratio_base = row.sigma_ratio_base;
    best_case.sigma_ratio_GPS = row.sigma_ratio_GPS;
    best_case.col_corr_base = row.col_corr_base;
    best_case.col_corr_GPS = row.col_corr_GPS;
end

function print_sorted(T, field_name, direction, title_text)
    [~, order] = sort(T.(field_name), direction);
    top = T(order(1:min(20, height(T))), :);
    fprintf('\n%s:\n', title_text);
    disp(top(:, {'case_name','l2','alpha2','L','gps_pattern_name','delta_type', ...
        'mean_E_base','mean_E_GPS','improvement_mean','rank_base','rank_GPS', ...
        'sigma_ratio_base','sigma_ratio_GPS','col_corr_base','col_corr_GPS'}));
end

function plot_focused_results(T, outputDir)
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    summary = best_row_per_case(T);
    labels = strcat("Case ", summary.case_name, " L=", string(summary.L));

    figure('Name', 'Focused mean equation error', 'Color', 'w');
    bar([summary.mean_E_base, summary.mean_E_GPS]);
    grid on;
    xticklabels(labels);
    xlabel('Focused two-path case');
    ylabel('Mean key-equation mismatch');
    title('Baseline vs GPS key-equation mismatch');
    legend({'Baseline', 'GPS'}, 'Location', 'northwest');
    saveas(gcf, fullfile(outputDir, ['focused_mean_error_' timestamp '.png']));

    figure('Name', 'Focused Phi sigma ratio', 'Color', 'w');
    bar([summary.sigma_ratio_base, summary.sigma_ratio_GPS]);
    grid on;
    xticklabels(labels);
    xlabel('Focused two-path case');
    ylabel('\sigma_{min}(\Phi)/\sigma_{max}(\Phi)');
    title('Baseline vs GPS Phi conditioning');
    legend({'Baseline', 'GPS'}, 'Location', 'southwest');
    saveas(gcf, fullfile(outputDir, ['focused_phi_sigma_ratio_' timestamp '.png']));

    figure('Name', 'Focused Phi column correlation', 'Color', 'w');
    bar([summary.col_corr_base, summary.col_corr_GPS]);
    grid on;
    ylim([0, 1]);
    xticklabels(labels);
    xlabel('Focused two-path case');
    ylabel('Column correlation');
    title('Baseline vs GPS Phi column correlation');
    legend({'Baseline', 'GPS'}, 'Location', 'northwest');
    saveas(gcf, fullfile(outputDir, ['focused_phi_col_corr_' timestamp '.png']));
end

function summary = best_row_per_case(T)
    case_names = unique(T.case_name, 'stable');
    summary = table();
    for idx = 1:numel(case_names)
        rows = T(T.case_name == case_names(idx), :);
        [~, best] = max(rows.improvement_mean);
        summary = [summary; rows(best, :)]; %#ok<AGROW>
    end
end
