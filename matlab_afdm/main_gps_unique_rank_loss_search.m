% MAIN_GPS_UNIQUE_RANK_LOSS_SEARCH
% 强化搜索：寻找 GPS 的 Phi(delta) 独有退秩/近退秩 case。
% 注意：本脚本不使用 total H_eff 的 rank 作为依据，只计算
% Phi(delta)=[H1*delta,H2*delta]，其中 H1/H2 是单位路径单径矩阵。

rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
setup_paths(rootDir);

cfgSearch.N_list = [16, 32, 64];
cfgSearch.V_list = [2, 4, 8];
cfgSearch.alpha_max_list = [1, 2, 3];
cfgSearch.max_exhaustive_bpsk = 4096;   % 全 2^16 可放开，但默认限量避免爆炸
cfgSearch.num_random_bpsk = 2000;
cfgSearch.num_random_qpsk = 2000;
cfgSearch.max_phi_rows = 200000;        % 只对 key-equation 强候选做 Phi 精算
cfgSearch.seed = 20260427;

fprintf('========== GPS unique Phi rank-loss search ==========\n');
fprintf('N=%s | V=%s | alpha_max=%s\n', ...
    mat2str(cfgSearch.N_list), mat2str(cfgSearch.V_list), mat2str(cfgSearch.alpha_max_list));
fprintf('Max Phi-verified rows: %d\n\n', cfgSearch.max_phi_rows);

rows = initialize_rows(cfgSearch.max_phi_rows);
rowCount = 0;

for N = cfgSearch.N_list
    for V = cfgSearch.V_list
        if mod(N, V) ~= 0
            continue;
        end

        M = N / V;
        gpsPatterns = build_all_patterns(2, V);
        proposedPatterns = build_all_patterns(3, V);

        for alpha_max = cfgSearch.alpha_max_list
            c1 = (2 * alpha_max + 1) / (2 * N);
            c2_base = sqrt(2) / (10 * N);
            proposed_delta = c2_base / 16;
            c2_values.base = c2_base;
            c2_values.small = 1 / (100 * N);
            c2_values.zero = 0;
            c2_values.ocdm = 1 / (2 * N);

            pathCases = build_path_cases(N, alpha_max);
            fprintf('N=%d V=%d M=%d alpha_max=%d | path cases=%d | GPS patterns=%d\n', ...
                N, V, M, alpha_max, numel(pathCases), size(gpsPatterns, 1));

            for caseIdx = 1:numel(pathCases)
                caseDef = pathCases(caseIdx);
                if caseDef.L == 0
                    continue;
                end

                proposedCache = build_proposed_cache(N, V, proposedPatterns, c2_base, proposed_delta);

                for gpsIdx = 1:size(gpsPatterns, 1)
                    gpsPattern = gpsPatterns(gpsIdx, :);
                    [c2_gps, d_gps] = build_c2m_gps_pattern(N, V, gpsPattern);
                    path_phase = exp(1i * 2 * pi / N * (caseDef.l2 - caseDef.l1) * caseDef.L);
                    deltaOptions = struct( ...
                        'max_exhaustive_bpsk', cfgSearch.max_exhaustive_bpsk, ...
                        'num_random_bpsk', cfgSearch.num_random_bpsk, ...
                        'num_random_qpsk', cfgSearch.num_random_qpsk, ...
                        'seed', cfgSearch.seed + 1000 * gpsIdx + caseIdx);
                    deltaSet = generate_delta_set_extended(N, V, caseDef.L, d_gps, path_phase, deltaOptions);

                    H_cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
                    H1_base = cached_H(H_cache, N, c1, c2_base, 'base', caseDef.l1, caseDef.alpha1);
                    H2_base = cached_H(H_cache, N, c1, c2_base, 'base', caseDef.l2, caseDef.alpha2);
                    H1_gps = cached_H(H_cache, N, c1, c2_gps, pattern_key('gps', gpsPattern), caseDef.l1, caseDef.alpha1);
                    H2_gps = cached_H(H_cache, N, c1, c2_gps, pattern_key('gps', gpsPattern), caseDef.l2, caseDef.alpha2);

                    for deltaIdx = 1:numel(deltaSet)
                        delta = deltaSet(deltaIdx).delta;
                        gpsErr = compute_bemani_equation_error(delta, d_gps, caseDef.L, ...
                            caseDef.l2 - caseDef.l1, N, c2_values);
                        proposedErr = best_proposed_error(delta, proposedCache, proposedPatterns, ...
                            caseDef.L, caseDef.l2 - caseDef.l1, N, c2_values);

                        if ~is_equation_candidate(gpsErr, proposedErr)
                            continue;
                        end
                        if rowCount >= cfgSearch.max_phi_rows
                            warning('Reached max_phi_rows=%d; remaining candidates are skipped.', cfgSearch.max_phi_rows);
                            break;
                        end

                        [c2_prop, ~] = build_c2m_proposed_pattern(N, V, proposedErr.pattern, c2_base, proposed_delta);
                        H1_prop = cached_H(H_cache, N, c1, c2_prop, pattern_key('prop', proposedErr.pattern), caseDef.l1, caseDef.alpha1);
                        H2_prop = cached_H(H_cache, N, c1, c2_prop, pattern_key('prop', proposedErr.pattern), caseDef.l2, caseDef.alpha2);

                        basePhi = evaluate_phi_metrics(H1_base, H2_base, delta);
                        gpsPhi = evaluate_phi_metrics(H1_gps, H2_gps, delta);
                        propPhi = evaluate_phi_metrics(H1_prop, H2_prop, delta);

                        rowCount = rowCount + 1;
                        rows = fill_row(rows, rowCount, N, V, M, alpha_max, c1, c2_base, ...
                            caseDef, gpsPattern, proposedErr.pattern, deltaSet(deltaIdx).name, ...
                            gpsErr, proposedErr, basePhi, gpsPhi, propPhi);
                    end

                    if rowCount >= cfgSearch.max_phi_rows
                        break;
                    end
                end

                if rowCount >= cfgSearch.max_phi_rows
                    break;
                end
            end
        end
    end
end

results_table = rows_to_table(rows, rowCount);
results_table = add_scores(results_table);
best_case = select_best_case(results_table);

outputDir = fullfile(fileparts(rootDir), 'results');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
save(fullfile(outputDir, 'gps_unique_rank_loss_search.mat'), 'results_table', 'best_case', 'cfgSearch');
writetable(results_table, fullfile(outputDir, 'gps_unique_rank_loss_search.csv'));
save(fullfile(outputDir, 'gps_unique_best_case.mat'), 'best_case');
plot_search_summary(results_table, outputDir, datestr(now, 'yyyymmdd_HHMMSS'));

print_top_cases(results_table);
fprintf('\nSaved search table and best_case to %s\n', outputDir);

function rows = initialize_rows(maxRows)
    names = {'N','V','M','alpha_max','c1','c2_base','l1','alpha1','l2','alpha2','loc1','loc2','L', ...
        'mean_E_base','mean_E_GPS','mean_E_proposed','min_E_base','min_E_GPS','min_E_proposed', ...
        'rank_base','rank_GPS','rank_prop','sigma_ratio_base','sigma_ratio_GPS','sigma_ratio_prop', ...
        'col_corr_base','col_corr_GPS','col_corr_prop'};
    for idx = 1:numel(names)
        rows.(names{idx}) = zeros(maxRows, 1);
    end
    rows.gps_pattern = strings(maxRows, 1);
    rows.proposed_pattern = strings(maxRows, 1);
    rows.delta_type = strings(maxRows, 1);
end

function pathCases = build_path_cases(N, alpha_max)
    pathCases = struct('l1', {}, 'alpha1', {}, 'l2', {}, 'alpha2', {}, 'loc1', {}, 'loc2', {}, 'L', {});
    factor = 2 * alpha_max + 1;
    lmax_search = min(8, floor(N / 4));
    for l2 = 0:lmax_search
        for alpha2 = -alpha_max:alpha_max
            if l2 == 0 && alpha2 == 0
                continue;
            end
            caseDef.l1 = 0;
            caseDef.alpha1 = 0;
            caseDef.l2 = l2;
            caseDef.alpha2 = alpha2;
            caseDef.loc1 = 0;
            caseDef.loc2 = alpha2 + factor * l2;
            caseDef.L = mod(caseDef.loc2, N);
            if caseDef.L ~= 0
                pathCases(end+1) = caseDef; %#ok<AGROW>
            end
        end
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

function cache = build_proposed_cache(N, V, patterns, c2_base, delta)
    cache = cell(size(patterns, 1), 2);
    for idx = 1:size(patterns, 1)
        [cache{idx, 1}, cache{idx, 2}] = build_c2m_proposed_pattern(N, V, patterns(idx, :), c2_base, delta);
    end
end

function tf = is_equation_candidate(gpsErr, proposedErr)
    impBase = gpsErr.mean_E_base / (gpsErr.mean_E_GPS + eps);
    impProp = proposedErr.mean_E_proposed / (gpsErr.mean_E_GPS + eps);
    tf = gpsErr.mean_E_GPS < 0.2 || gpsErr.min_E_GPS < 1e-8 || impBase > 5 || impProp > 3;
end

function proposedErr = best_proposed_error(delta, proposedCache, proposedPatterns, L, ldiff, N, c2_values)
    bestMean = inf;
    bestMetrics = [];
    bestIdx = 1;
    for idx = 1:size(proposedPatterns, 1)
        metrics = compute_bemani_equation_error(delta, proposedCache{idx, 2}, L, ldiff, N, c2_values);
        if metrics.mean_E_GPS < bestMean
            bestMean = metrics.mean_E_GPS;
            bestMetrics = metrics;
            bestIdx = idx;
        end
    end
    proposedErr.min_E_proposed = bestMetrics.min_E_GPS;
    proposedErr.mean_E_proposed = bestMetrics.mean_E_GPS;
    proposedErr.pattern = proposedPatterns(bestIdx, :);
end

function H = cached_H(H_cache, N, c1, c2m, keyPrefix, l, alpha)
    key = sprintf('%s_l%d_a%d', keyPrefix, l, alpha);
    if isKey(H_cache, key)
        H = H_cache(key);
        return;
    end
    H = build_H_path_general_c2m(N, c1, c2m, l, alpha);
    H_cache(key) = H;
end

function key = pattern_key(prefix, pattern)
    key = sprintf('%s_%s', prefix, sprintf('%d', pattern));
end

function rows = fill_row(rows, idx, N, V, M, alpha_max, c1, c2_base, caseDef, gpsPattern, propPattern, deltaName, gpsErr, propErr, basePhi, gpsPhi, propPhi)
    rows.N(idx) = N;
    rows.V(idx) = V;
    rows.M(idx) = M;
    rows.alpha_max(idx) = alpha_max;
    rows.c1(idx) = c1;
    rows.c2_base(idx) = c2_base;
    rows.l1(idx) = caseDef.l1;
    rows.alpha1(idx) = caseDef.alpha1;
    rows.l2(idx) = caseDef.l2;
    rows.alpha2(idx) = caseDef.alpha2;
    rows.loc1(idx) = caseDef.loc1;
    rows.loc2(idx) = caseDef.loc2;
    rows.L(idx) = caseDef.L;
    rows.gps_pattern(idx) = string(mat2str(gpsPattern));
    rows.proposed_pattern(idx) = string(mat2str(propPattern));
    rows.delta_type(idx) = string(deltaName);
    rows.mean_E_base(idx) = gpsErr.mean_E_base;
    rows.mean_E_GPS(idx) = gpsErr.mean_E_GPS;
    rows.mean_E_proposed(idx) = propErr.mean_E_proposed;
    rows.min_E_base(idx) = gpsErr.min_E_base;
    rows.min_E_GPS(idx) = gpsErr.min_E_GPS;
    rows.min_E_proposed(idx) = propErr.min_E_proposed;
    rows.rank_base(idx) = basePhi.rank;
    rows.rank_GPS(idx) = gpsPhi.rank;
    rows.rank_prop(idx) = propPhi.rank;
    rows.sigma_ratio_base(idx) = basePhi.sigma_ratio;
    rows.sigma_ratio_GPS(idx) = gpsPhi.sigma_ratio;
    rows.sigma_ratio_prop(idx) = propPhi.sigma_ratio;
    rows.col_corr_base(idx) = basePhi.col_corr;
    rows.col_corr_GPS(idx) = gpsPhi.col_corr;
    rows.col_corr_prop(idx) = propPhi.col_corr;
end

function T = rows_to_table(rows, rowCount)
    fields = fieldnames(rows);
    vars = cell(1, numel(fields));
    for idx = 1:numel(fields)
        vars{idx} = rows.(fields{idx})(1:rowCount);
    end
    T = table(vars{:}, 'VariableNames', fields);
end

function T = add_scores(T)
    T.improvement_GPS_over_base_mean = T.mean_E_base ./ (T.mean_E_GPS + eps);
    T.improvement_GPS_over_prop_mean = T.mean_E_proposed ./ (T.mean_E_GPS + eps);
    T.score1 = (T.sigma_ratio_base + T.sigma_ratio_prop) ./ (2 * (T.sigma_ratio_GPS + eps));
    T.score2 = T.col_corr_GPS - max(T.col_corr_base, T.col_corr_prop);
    T.score3 = T.improvement_GPS_over_base_mean;
    T.score4 = T.improvement_GPS_over_prop_mean;
end

function best_case = select_best_case(T)
    strict = T.rank_GPS < 2 & T.rank_base == 2 & T.rank_prop == 2;
    near = T.sigma_ratio_GPS < 1e-3 & T.sigma_ratio_base > 0.5 & T.sigma_ratio_prop > 0.5;
    corr = T.col_corr_GPS > 0.999 & T.col_corr_base < 0.5 & T.col_corr_prop < 0.5;
    if any(strict)
        C = T(strict, :);
        [~, idx] = min(C.sigma_ratio_GPS);
    elseif any(near)
        C = T(near, :);
        [~, idx] = min(C.sigma_ratio_GPS);
    elseif any(corr)
        C = T(corr, :);
        [~, idx] = max(C.col_corr_GPS);
    else
        C = T;
        [~, idx] = max(C.score1);
    end
    best_case = table2struct(C(idx, :));
end

function print_top_cases(T)
    print_table(T(T.rank_GPS < 2 & T.rank_base == 2 & T.rank_prop == 2, :), 'Strict GPS-only rank loss');
    print_sorted(T, 'sigma_ratio_GPS', 'ascend', 'Smallest GPS sigma ratio');
    print_sorted(T, 'score1', 'descend', 'Largest sigma-ratio separation score');
    C = T(T.col_corr_base < 0.5 & T.col_corr_prop < 0.5, :);
    print_sorted(C, 'col_corr_GPS', 'descend', 'Largest GPS col corr with safe baseline/proposed');
    print_sorted(T, 'score3', 'descend', 'Largest key-equation GPS/base improvement');
end

function print_sorted(T, fieldName, direction, titleText)
    if isempty(T)
        fprintf('\n%s: none\n', titleText);
        return;
    end
    [~, order] = sort(T.(fieldName), direction);
    print_table(T(order(1:min(50, numel(order))), :), titleText);
end

function print_table(T, titleText)
    fprintf('\n%s (%d rows):\n', titleText, height(T));
    if isempty(T)
        return;
    end
    vars = {'N','V','M','alpha_max','l2','alpha2','L','gps_pattern','delta_type', ...
        'mean_E_base','mean_E_GPS','mean_E_proposed','rank_base','rank_GPS','rank_prop', ...
        'sigma_ratio_base','sigma_ratio_GPS','sigma_ratio_prop','col_corr_base','col_corr_GPS','col_corr_prop'};
    disp(T(1:min(50, height(T)), vars));
end
