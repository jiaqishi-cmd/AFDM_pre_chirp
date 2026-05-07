% MAIN_BEMANI_GPS_KEY_EQUATION_SEARCH
% Search Bemani two-path key-equation matches for Yuan GPS-AFDM c2,m.

rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
setup_paths(rootDir);

N = 64;
alpha_max = 3;
c1 = (2 * alpha_max + 1) / (2 * N);
V = 4;
W = 2; %#ok<NASGU>
M = N / V;

c2_base = sqrt(2) / (10 * N);
c2_values.base = c2_base;
c2_values.small = 1 / (100 * N);
c2_values.zero = 0;
c2_values.ocdm = 1 / (2 * N);

patterns = { ...
    'pattern1', [1 1 1 1]; ...
    'pattern2', [1 2 1 2]; ...
    'pattern3', [1 1 2 2]; ...
    'pattern4', [1 2 2 1]; ...
    'pattern5', [2 1 2 1]; ...
    'pattern6', [2 2 1 1] ...
    };

proposed_patterns = build_all_patterns(3, V);
proposed_delta = c2_base / 16;

delta_set = generate_delta_set_for_equation_search(N, V, 1000, 1000, 20260427);
path_cases = build_path_cases(alpha_max, N);

numRows = numel(path_cases) * size(patterns, 1) * numel(delta_set);
fprintf('========== Bemani/GPS key-equation search ==========\n');
fprintf('N=%d, alpha_max=%d, c1=%.10g, V=%d, M=%d\n', N, alpha_max, c1, V, M);
fprintf('Path cases=%d, GPS patterns=%d, deltas=%d, rows=%d\n\n', ...
    numel(path_cases), size(patterns, 1), numel(delta_set), numRows);

rows = initialize_rows(numRows);
row = 0;

gps_cache = cell(size(patterns, 1), 2);
for pattern_idx = 1:size(patterns, 1)
    [gps_cache{pattern_idx, 1}, gps_cache{pattern_idx, 2}] = ...
        build_c2m_gps_pattern(N, V, patterns{pattern_idx, 2});
end

proposed_cache = cell(size(proposed_patterns, 1), 2);
for proposed_idx = 1:size(proposed_patterns, 1)
    [proposed_cache{proposed_idx, 1}, proposed_cache{proposed_idx, 2}] = ...
        build_c2m_proposed_pattern(N, V, proposed_patterns(proposed_idx, :), ...
        c2_base, proposed_delta);
end

for case_idx = 1:numel(path_cases)
    caseDef = path_cases(case_idx);
    if caseDef.L == 0
        continue;
    end

    proposed_by_delta = repmat(struct( ...
        'min_E_proposed', 0, ...
        'mean_E_proposed', 0, ...
        'max_E_proposed', 0, ...
        'pattern', zeros(1, V), ...
        'pattern_name', ''), numel(delta_set), 1);
    for delta_idx = 1:numel(delta_set)
        proposed_by_delta(delta_idx) = best_proposed_equation_error( ...
            delta_set(delta_idx).delta, proposed_cache, proposed_patterns, caseDef.L, ...
            caseDef.l2 - caseDef.l1, N, c2_values);
    end

    for pattern_idx = 1:size(patterns, 1)
        d_gps = gps_cache{pattern_idx, 2};

        for delta_idx = 1:numel(delta_set)
            delta = delta_set(delta_idx).delta;
            proposed_metrics = proposed_by_delta(delta_idx);
            metrics = compute_bemani_equation_error( ...
                delta, d_gps, caseDef.L, caseDef.l2 - caseDef.l1, N, c2_values);

            row = row + 1;
            rows.case_id(row) = string(caseDef.name);
            rows.l1(row) = caseDef.l1;
            rows.alpha1(row) = caseDef.alpha1;
            rows.l2(row) = caseDef.l2;
            rows.alpha2(row) = caseDef.alpha2;
            rows.loc1(row) = caseDef.loc1;
            rows.loc2(row) = caseDef.loc2;
            rows.L(row) = caseDef.L;
            rows.gps_pattern_name(row) = string(patterns{pattern_idx, 1});
            rows.delta_type(row) = string(delta_set(delta_idx).name);

            rows.min_E_base(row) = metrics.min_E_base;
            rows.mean_E_base(row) = metrics.mean_E_base;
            rows.max_E_base(row) = metrics.max_E_base;
            rows.min_E_GPS(row) = metrics.min_E_GPS;
            rows.mean_E_GPS(row) = metrics.mean_E_GPS;
            rows.max_E_GPS(row) = metrics.max_E_GPS;
            rows.min_E_proposed(row) = proposed_metrics.min_E_proposed;
            rows.mean_E_proposed(row) = proposed_metrics.mean_E_proposed;
            rows.max_E_proposed(row) = proposed_metrics.max_E_proposed;
            rows.proposed_pattern_name(row) = string(proposed_metrics.pattern_name);
            rows.improvement_min(row) = metrics.min_E_base / (metrics.min_E_GPS + eps);
            rows.improvement_mean(row) = metrics.mean_E_base / (metrics.mean_E_GPS + eps);
            rows.improvement_mean_proposed(row) = ...
                metrics.mean_E_base / (proposed_metrics.mean_E_proposed + eps);
            rows.min_E_small(row) = metrics.min_E_small;
            rows.mean_E_small(row) = metrics.mean_E_small;
            rows.min_E_zero(row) = metrics.min_E_zero;
            rows.mean_E_zero(row) = metrics.mean_E_zero;
            rows.min_E_ocdm(row) = metrics.min_E_ocdm;
            rows.mean_E_ocdm(row) = metrics.mean_E_ocdm;
        end
    end

    if mod(case_idx, 10) == 0 || case_idx == numel(path_cases)
        fprintf('Path case %3d/%d complete\n', case_idx, numel(path_cases));
    end
end

results_table = struct_to_table(rows, row);
suspect = results_table.min_E_GPS < 1e-8 | ...
    results_table.min_E_GPS < 1e-6 | ...
    results_table.improvement_mean > 10 | ...
    results_table.improvement_min > 100 | ...
    startsWith(results_table.case_id, "Case_");

fprintf('\nSuspicious rows for Phi(delta) verification: %d\n', nnz(suspect));
results_table = verify_suspect_phi(results_table, suspect, patterns, delta_set, ...
    N, V, c1, c2_base, proposed_delta);

outputDir = fullfile(fileparts(rootDir), 'results');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
matPath = fullfile(outputDir, 'bemani_gps_key_equation_search.mat');
csvPath = fullfile(outputDir, 'bemani_gps_key_equation_search.csv');
save(matPath, 'results_table');
writetable(results_table, csvPath);

print_ranked_tables(results_table);
phi_plot = select_phi_plot(results_table, patterns, delta_set, N, V, c1, c2_base);
plot_bemani_gps_search_results(results_table, outputDir, timestamp, phi_plot);

fprintf('\nSaved MAT to %s\n', matPath);
fprintf('Saved CSV to %s\n', csvPath);

function cases = build_path_cases(alpha_max, N)
    l1 = 0;
    alpha1 = 0;
    loc1 = 0;
    factor = 2 * alpha_max + 1;

    cases = struct('name', {}, 'l1', {}, 'alpha1', {}, 'l2', {}, 'alpha2', {}, ...
        'loc1', {}, 'loc2', {}, 'L_raw', {}, 'L', {});

    for l2 = 0:8
        for alpha2 = -alpha_max:alpha_max
            if l2 == l1 && alpha2 == alpha1
                continue;
            end
            loc2 = alpha2 + factor * l2;
            cases(end+1) = make_case(sprintf('enum_l%d_a%+d', l2, alpha2), ...
                l1, alpha1, l2, alpha2, loc1, loc2, N); %#ok<AGROW>
        end
    end

    focus = { ...
        'Case_A', 2, 2; ...
        'Case_B', 5, -3; ...
        'Case_C', 2, 1; ...
        'Case_D', 2, 3 ...
        };
    for idx = 1:size(focus, 1)
        l2 = focus{idx, 2};
        alpha2 = focus{idx, 3};
        loc2 = alpha2 + factor * l2;
        cases(end+1) = make_case(focus{idx, 1}, ...
            l1, alpha1, l2, alpha2, loc1, loc2, N); %#ok<AGROW>
    end
end

function patterns = build_all_patterns(base, width)
    count = base^width;
    patterns = zeros(count, width);
    for row = 1:count
        value = row - 1;
        for idx = 1:width
            patterns(row, idx) = mod(value, base) + 1;
            value = floor(value / base);
        end
    end
end

function proposed_metrics = best_proposed_equation_error(delta, proposed_cache, proposed_patterns, L, ldiff, N, c2_values)
    best_mean = inf;
    best_metrics = [];
    best_idx = 1;

    for proposed_idx = 1:size(proposed_patterns, 1)
        d_proposed = proposed_cache{proposed_idx, 2};
        metrics = compute_bemani_equation_error(delta, d_proposed, L, ldiff, N, c2_values);
        if metrics.mean_E_GPS < best_mean
            best_mean = metrics.mean_E_GPS;
            best_metrics = metrics;
            best_idx = proposed_idx;
        end
    end

    proposed_metrics.min_E_proposed = best_metrics.min_E_GPS;
    proposed_metrics.mean_E_proposed = best_metrics.mean_E_GPS;
    proposed_metrics.max_E_proposed = best_metrics.max_E_GPS;
    proposed_metrics.pattern = proposed_patterns(best_idx, :);
    proposed_metrics.pattern_name = sprintf('proposed_%s', strrep(mat2str(proposed_metrics.pattern), ' ', ''));
end

function caseDef = make_case(name, l1, alpha1, l2, alpha2, loc1, loc2, N)
    caseDef.name = name;
    caseDef.l1 = l1;
    caseDef.alpha1 = alpha1;
    caseDef.l2 = l2;
    caseDef.alpha2 = alpha2;
    caseDef.loc1 = loc1;
    caseDef.loc2 = loc2;
    caseDef.L_raw = loc2 - loc1;
    caseDef.L = mod(caseDef.L_raw, N);
end

function rows = initialize_rows(numRows)
    rows.case_id = strings(numRows, 1);
    rows.l1 = zeros(numRows, 1);
    rows.alpha1 = zeros(numRows, 1);
    rows.l2 = zeros(numRows, 1);
    rows.alpha2 = zeros(numRows, 1);
    rows.loc1 = zeros(numRows, 1);
    rows.loc2 = zeros(numRows, 1);
    rows.L = zeros(numRows, 1);
    rows.gps_pattern_name = strings(numRows, 1);
    rows.proposed_pattern_name = strings(numRows, 1);
    rows.delta_type = strings(numRows, 1);
    metricNames = {'min_E_base','mean_E_base','max_E_base','min_E_GPS','mean_E_GPS', ...
        'max_E_GPS','min_E_proposed','mean_E_proposed','max_E_proposed', ...
        'improvement_min','improvement_mean','improvement_mean_proposed', ...
        'min_E_small','mean_E_small','min_E_zero','mean_E_zero','min_E_ocdm','mean_E_ocdm'};
    for idx = 1:numel(metricNames)
        rows.(metricNames{idx}) = zeros(numRows, 1);
    end
end

function T = struct_to_table(rows, rowCount)
    names = fieldnames(rows);
    args = cell(1, 2 * numel(names));
    for idx = 1:numel(names)
        args{2 * idx - 1} = names{idx};
        args{2 * idx} = rows.(names{idx})(1:rowCount);
    end
    T = table(args{2:2:end}, 'VariableNames', args(1:2:end));

    phiNames = {'rank_base','sigma_min_base','sigma_ratio_base','col_corr_base', ...
        'rank_GPS','sigma_min_GPS','sigma_ratio_GPS','col_corr_GPS', ...
        'rank_proposed','sigma_min_proposed','sigma_ratio_proposed','col_corr_proposed'};
    for idx = 1:numel(phiNames)
        T.(phiNames{idx}) = NaN(height(T), 1);
    end
end

function T = verify_suspect_phi(T, suspect, patterns, delta_set, N, V, c1, c2_base, proposed_delta)
    suspect_idx = find(suspect);
    H_cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
    delta_names = strings(numel(delta_set), 1);
    for delta_idx = 1:numel(delta_set)
        delta_names(delta_idx) = string(delta_set(delta_idx).name);
    end

    for idx = 1:numel(suspect_idx)
        rowIdx = suspect_idx(idx);
        pattern = pattern_from_name(patterns, T.gps_pattern_name(rowIdx));
        [c2_gps, ~] = build_c2m_gps_pattern(N, V, pattern);
        proposed_pattern = parse_proposed_pattern(T.proposed_pattern_name(rowIdx));
        [c2_proposed, ~] = build_c2m_proposed_pattern(N, V, proposed_pattern, c2_base, proposed_delta);
        delta = delta_from_name_fast(delta_set, delta_names, T.delta_type(rowIdx));

        H1_base = cached_H_path(H_cache, N, c1, c2_base, 'base', ...
            T.l1(rowIdx), T.alpha1(rowIdx));
        H2_base = cached_H_path(H_cache, N, c1, c2_base, 'base', ...
            T.l2(rowIdx), T.alpha2(rowIdx));
        H1_gps = cached_H_path(H_cache, N, c1, c2_gps, char(T.gps_pattern_name(rowIdx)), ...
            T.l1(rowIdx), T.alpha1(rowIdx));
        H2_gps = cached_H_path(H_cache, N, c1, c2_gps, char(T.gps_pattern_name(rowIdx)), ...
            T.l2(rowIdx), T.alpha2(rowIdx));
        H1_proposed = cached_H_path(H_cache, N, c1, c2_proposed, char(T.proposed_pattern_name(rowIdx)), ...
            T.l1(rowIdx), T.alpha1(rowIdx));
        H2_proposed = cached_H_path(H_cache, N, c1, c2_proposed, char(T.proposed_pattern_name(rowIdx)), ...
            T.l2(rowIdx), T.alpha2(rowIdx));

        baseMetrics = evaluate_phi_metrics(H1_base, H2_base, delta);
        gpsMetrics = evaluate_phi_metrics(H1_gps, H2_gps, delta);
        proposedMetrics = evaluate_phi_metrics(H1_proposed, H2_proposed, delta);

        T.rank_base(rowIdx) = baseMetrics.rank;
        T.sigma_min_base(rowIdx) = baseMetrics.sigma_min;
        T.sigma_ratio_base(rowIdx) = baseMetrics.sigma_ratio;
        T.col_corr_base(rowIdx) = baseMetrics.col_corr;
        T.rank_GPS(rowIdx) = gpsMetrics.rank;
        T.sigma_min_GPS(rowIdx) = gpsMetrics.sigma_min;
        T.sigma_ratio_GPS(rowIdx) = gpsMetrics.sigma_ratio;
        T.col_corr_GPS(rowIdx) = gpsMetrics.col_corr;
        T.rank_proposed(rowIdx) = proposedMetrics.rank;
        T.sigma_min_proposed(rowIdx) = proposedMetrics.sigma_min;
        T.sigma_ratio_proposed(rowIdx) = proposedMetrics.sigma_ratio;
        T.col_corr_proposed(rowIdx) = proposedMetrics.col_corr;

        if mod(idx, 1000) == 0 || idx == numel(suspect_idx)
            fprintf('Phi verification %6d/%d complete\n', idx, numel(suspect_idx));
        end
    end
end

function pattern = parse_proposed_pattern(name)
    text = char(name);
    tokens = regexp(text, '\[(\d)(\d)(\d)(\d)\]', 'tokens', 'once');
    if isempty(tokens)
        error('Unknown proposed pattern name: %s', text);
    end
    pattern = zeros(1, numel(tokens));
    for idx = 1:numel(tokens)
        pattern(idx) = str2double(tokens{idx});
    end
end

function H = cached_H_path(H_cache, N, c1, c2m, c2Key, l, alpha)
    key = sprintf('%s_l%d_a%d', c2Key, l, alpha);
    if isKey(H_cache, key)
        H = H_cache(key);
        return;
    end

    H = build_H_path_general_c2m(N, c1, c2m, l, alpha);
    H_cache(key) = H;
end

function pattern = pattern_from_name(patterns, name)
    for idx = 1:size(patterns, 1)
        if string(patterns{idx, 1}) == string(name)
            pattern = patterns{idx, 2};
            return;
        end
    end
    error('Unknown GPS pattern: %s', name);
end

function delta = delta_from_name(delta_set, name)
    for idx = 1:numel(delta_set)
        if string(delta_set(idx).name) == string(name)
            delta = delta_set(idx).delta;
            return;
        end
    end
    error('Unknown delta: %s', name);
end

function delta = delta_from_name_fast(delta_set, delta_names, name)
    idx = find(delta_names == string(name), 1);
    if isempty(idx)
        error('Unknown delta: %s', name);
    end
    delta = delta_set(idx).delta;
end

function print_ranked_tables(T)
    fprintf('\nTop 20 by min_E_GPS ascending:\n');
    print_top(T, 'min_E_GPS', 'ascend', 20);

    fprintf('\nTop 20 by improvement_mean descending:\n');
    print_top(T, 'improvement_mean', 'descend', 20);

    phiRows = T(~isnan(T.sigma_ratio_GPS) & ~isnan(T.sigma_ratio_base), :);
    if ~isempty(phiRows)
        phiRows.sigma_ratio_gps_over_base = ...
            phiRows.sigma_ratio_GPS ./ max(phiRows.sigma_ratio_base, eps);

        fprintf('\nTop 20 by sigma_ratio_GPS/sigma_ratio_base ascending:\n');
        print_top(phiRows, 'sigma_ratio_gps_over_base', 'ascend', 20);

        fprintf('\nTop 20 by col_corr_GPS descending:\n');
        print_top(phiRows, 'col_corr_GPS', 'descend', 20);
    end
end

function print_top(T, fieldName, direction, n)
    [~, order] = sort(T.(fieldName), direction);
    order = order(1:min(n, numel(order)));
    vars = {'case_id','l2','alpha2','L','gps_pattern_name','delta_type', ...
        'min_E_base','mean_E_base','min_E_GPS','mean_E_GPS','mean_E_proposed', ...
        'improvement_min','improvement_mean','improvement_mean_proposed', ...
        'sigma_ratio_base','sigma_ratio_GPS','sigma_ratio_proposed','col_corr_GPS','col_corr_proposed'};
    vars = vars(ismember(vars, T.Properties.VariableNames));
    disp(T(order, vars));
end

function phi_plot = select_phi_plot(T, patterns, delta_set, N, V, c1, c2_base)
    phi_plot = struct('phi1', [], 'phi2', [], 'title', '');
    phiRows = T(~isnan(T.col_corr_GPS), :);
    if isempty(phiRows)
        return;
    end
    [~, idx] = max(phiRows.col_corr_GPS);
    row = phiRows(idx, :);

    pattern = pattern_from_name(patterns, row.gps_pattern_name);
    [c2_gps, ~] = build_c2m_gps_pattern(N, V, pattern);
    delta = delta_from_name(delta_set, row.delta_type);
    H1_gps = build_H_path_general_c2m(N, c1, c2_gps, row.l1, row.alpha1);
    H2_gps = build_H_path_general_c2m(N, c1, c2_gps, row.l2, row.alpha2);
    metrics = evaluate_phi_metrics(H1_gps, H2_gps, delta);

    phi_plot.phi1 = metrics.phi1;
    phi_plot.phi2 = metrics.phi2;
    phi_plot.title = sprintf('%s %s %s colcorr %.6f', ...
        row.case_id, row.gps_pattern_name, row.delta_type, row.col_corr_GPS);
end
