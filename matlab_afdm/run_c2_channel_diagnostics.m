function results = run_c2_channel_diagnostics(maxCasesPerScenario)
%RUN_C2_CHANNEL_DIAGNOSTICS Search channels that expose c2-selection tradeoffs.
%   results = run_c2_channel_diagnostics(maxCasesPerScenario)

    rootDir = fileparts(mfilename('fullpath'));
    addpath(rootDir);
    setup_paths(rootDir);

    if nargin < 1 || isempty(maxCasesPerScenario)
        maxCasesPerScenario = inf;
    end

    base_config = afdm_config();
    base_config.channel.multipath = true;
    base_config.channel.add_noise = false;
    base_config.channel.snr_db = 100;

    rng(base_config.simulation.random_seed, 'twister');
    [symbols, bits] = random_data_generator( ...
        base_config.waveform.NumSubcarriers, ...
        base_config.modulation.M_mod, ...
        base_config.modulation.modType);

    schemes = {'baseline', 'paper_grouping', 'proposed_grouping'};
    labels = {'Baseline', 'GPS', 'Proposed'};
    profiles = build_scheme_profiles(base_config, symbols, schemes);

    scenarios = { ...
        build_two_path_cases(base_config), ...
        build_equal_grid_cases(base_config), ...
        build_fractional_doppler_cases(base_config) ...
        };

    results.schemes = schemes;
    results.labels = labels;
    results.symbols = symbols;
    results.bits = bits;
    results.scenarios = cell(1, numel(scenarios));

    fprintf('========== c2/channel diagnostics ==========\n');
    for scenario_idx = 1:numel(scenarios)
        scenario = scenarios{scenario_idx};
        if isfinite(maxCasesPerScenario)
            scenario.cases = scenario.cases(1:min(maxCasesPerScenario, numel(scenario.cases)));
        end

        fprintf('\nScenario: %s | cases: %d\n', scenario.name, numel(scenario.cases));
        scenario_results = evaluate_scenario(base_config, profiles, scenario);
        print_top_cases(scenario_results, labels, 5);
        results.scenarios{scenario_idx} = scenario_results;
    end

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    results.output_dir = outputDir;
    results.timestamp = timestamp;

    plot_diagnostic_results(results);

    matPath = fullfile(outputDir, ['c2_channel_diagnostics_' timestamp '.mat']);
    save(matPath, 'results');
    fprintf('\nSaved diagnostics to %s\n', matPath);
end

function profiles = build_scheme_profiles(base_config, symbols, schemes)
    profiles = cell(1, numel(schemes));

    for scheme_idx = 1:numel(schemes)
        cfg = apply_pre_chirp_scheme(base_config, schemes{scheme_idx});
        profiles{scheme_idx} = select_pre_chirp_for_symbols(symbols, cfg);
    end
end

function scenario = build_two_path_cases(config)
    maxDelay = max(1, min(config.waveform.CPPLength, 3));
    dopplers = -3:0.5:3;
    phases = 0:pi/8:(2*pi - pi/8);

    cases = {};
    for delay = 1:maxDelay
        for doppler = dopplers
            for phase = phases
                caseDef.name = sprintf('two_path_l%d_nu%+.1f_phi%.2f', delay, doppler, phase);
                caseDef.delay_taps = [0, delay];
                caseDef.doppler_taps = [0, doppler];
                caseDef.chan_coef = [1, exp(1i * phase)] / sqrt(2);
                cases{end+1} = caseDef; %#ok<AGROW>
            end
        end
    end

    scenario.name = 'near-collision two-path phase sweep';
    scenario.cases = cases;
end

function scenario = build_equal_grid_cases(config)
    delays = 0:config.channel.max_integer_delay;
    dopplers = -config.channel.max_integer_doppler:config.channel.max_integer_doppler;
    [delay_grid, doppler_grid] = ndgrid(delays, dopplers);
    delay_taps = delay_grid(:).';
    doppler_taps = doppler_grid(:).';
    paths = numel(delay_taps);
    phase_slopes = 0:pi/12:(2*pi - pi/12);

    cases = {};
    path_idx = 0:paths-1;
    for slope = phase_slopes
        phases = slope * path_idx;
        caseDef.name = sprintf('equal_grid_slope%.2f', slope);
        caseDef.delay_taps = delay_taps;
        caseDef.doppler_taps = doppler_taps;
        caseDef.chan_coef = exp(1i * phases) / sqrt(paths);
        cases{end+1} = caseDef; %#ok<AGROW>
    end

    scenario.name = 'equal-power 21-path phase-ramp grid';
    scenario.cases = cases;
end

function scenario = build_fractional_doppler_cases(~)
    fractional_sets = { ...
        [-2.5, 0.5, 2.5], ...
        [-2.25, 0.25, 2.75], ...
        [-1.5, 0.5, 1.5], ...
        [-2.75, -0.25, 2.25] ...
        };
    phases = 0:pi/8:(2*pi - pi/8);

    cases = {};
    for set_idx = 1:numel(fractional_sets)
        for phase1 = phases
            for phase2 = phases
                caseDef.name = sprintf('frac_set%d_phi%.2f_%.2f', set_idx, phase1, phase2);
                caseDef.delay_taps = [0, 1, 2];
                caseDef.doppler_taps = fractional_sets{set_idx};
                caseDef.chan_coef = [1, exp(1i * phase1), exp(1i * phase2)] / sqrt(3);
                cases{end+1} = caseDef; %#ok<AGROW>
            end
        end
    end

    scenario.name = 'fractional-Doppler three-path phase sweep';
    scenario.cases = cases;
end

function scenario_results = evaluate_scenario(base_config, profiles, scenario)
    num_cases = numel(scenario.cases);
    num_schemes = numel(profiles);

    min_sv = zeros(num_cases, num_schemes);
    cond_num = zeros(num_cases, num_schemes);
    offdiag_ratio = zeros(num_cases, num_schemes);
    single_symbol_distance = zeros(num_cases, num_schemes);
    two_symbol_distance = zeros(num_cases, num_schemes);

    for case_idx = 1:num_cases
        cfg = apply_case_to_config(base_config, scenario.cases{case_idx});

        for scheme_idx = 1:num_schemes
            H_eff = estimate_effective_channel( ...
                cfg.waveform.NumSubcarriers, ...
                cfg.waveform.c1, ...
                profiles{scheme_idx}.c2, ...
                cfg.channel.chan_coef, ...
                cfg.channel.delay_taps, ...
                cfg.channel.doppler_freq);

            svals = svd(H_eff);
            min_sv(case_idx, scheme_idx) = min(svals);
            cond_num(case_idx, scheme_idx) = max(svals) / max(min(svals), eps);
            offdiag_ratio(case_idx, scheme_idx) = norm(H_eff - diag(diag(H_eff)), 'fro') / norm(H_eff, 'fro');

            distance_metrics = finite_constellation_distance_metrics( ...
                H_eff, ...
                cfg.modulation.M_mod, ...
                cfg.modulation.modType);
            single_symbol_distance(case_idx, scheme_idx) = distance_metrics.single_symbol_distance;
            two_symbol_distance(case_idx, scheme_idx) = distance_metrics.two_symbol_distance;
        end
    end

    gps_idx = 2;
    proposed_idx = 3;
    score = log10((single_symbol_distance(:, proposed_idx) + eps) ./ ...
        (single_symbol_distance(:, gps_idx) + eps)) + ...
        log10((two_symbol_distance(:, proposed_idx) + eps) ./ ...
        (two_symbol_distance(:, gps_idx) + eps));

    scenario_results.name = scenario.name;
    scenario_results.cases = scenario.cases;
    scenario_results.min_sv = min_sv;
    scenario_results.cond_num = cond_num;
    scenario_results.offdiag_ratio = offdiag_ratio;
    scenario_results.single_symbol_distance = single_symbol_distance;
    scenario_results.two_symbol_distance = two_symbol_distance;
    scenario_results.gps_exposure_score = score;
end

function cfg = apply_case_to_config(base_config, caseDef)
    cfg = base_config;
    cfg.channel.profile = caseDef.name;
    cfg.channel.delay_taps = caseDef.delay_taps;
    cfg.channel.doppler_taps = caseDef.doppler_taps;
    cfg.channel.doppler_freq = caseDef.doppler_taps / cfg.waveform.NumSubcarriers;
    cfg.channel.chan_coef = caseDef.chan_coef / norm(caseDef.chan_coef);
end

function print_top_cases(scenario_results, labels, topK)
    [~, order] = sort(scenario_results.gps_exposure_score, 'descend');
    topK = min(topK, numel(order));

    for rank_idx = 1:topK
        case_idx = order(rank_idx);
        caseDef = scenario_results.cases{case_idx};
        fprintf('  #%d %s | score %.3f\n', ...
            rank_idx, caseDef.name, scenario_results.gps_exposure_score(case_idx));
        fprintf('     min_sv: %s\n', format_metric_row(scenario_results.min_sv(case_idx, :), labels));
        fprintf('     cond  : %s\n', format_metric_row(scenario_results.cond_num(case_idx, :), labels));
        fprintf('     offdiag: %s\n', format_metric_row(scenario_results.offdiag_ratio(case_idx, :), labels));
        fprintf('     d1    : %s\n', format_metric_row(scenario_results.single_symbol_distance(case_idx, :), labels));
        fprintf('     d2    : %s\n', format_metric_row(scenario_results.two_symbol_distance(case_idx, :), labels));
    end
end

function text = format_metric_row(values, labels)
    parts = cell(1, numel(values));
    for idx = 1:numel(values)
        parts{idx} = sprintf('%s=%.3g', labels{idx}, values(idx));
    end
    text = strjoin(parts, ', ');
end

function plot_diagnostic_results(results)
    for scenario_idx = 1:numel(results.scenarios)
        scenario = results.scenarios{scenario_idx};
        safeName = regexprep(lower(scenario.name), '[^a-z0-9]+', '_');
        safeName = regexprep(safeName, '^_|_$', '');

        plot_score_ranking(scenario, results.output_dir, results.timestamp, safeName);
        plot_distance_comparison(scenario, results.labels, results.output_dir, results.timestamp, safeName);
        plot_distance_ratio_histogram(scenario, results.output_dir, results.timestamp, safeName);
    end
end

function plot_score_ranking(scenario, outputDir, timestamp, safeName)
    [sortedScore, order] = sort(scenario.gps_exposure_score, 'descend');
    topK = min(30, numel(sortedScore));
    topOrder = order(1:topK);

    figure('Name', ['GPS exposure score - ' scenario.name], 'Color', 'w');
    bar(sortedScore(1:topK));
    grid on;
    xlabel('Ranked channel case');
    ylabel('GPS exposure score');
    title(['Top GPS exposure cases: ' scenario.name], 'Interpreter', 'none');
    xticks(1:topK);
    xticklabels(short_case_names(scenario.cases(topOrder)));
    xtickangle(45);
    saveas(gcf, fullfile(outputDir, ['c2_score_' safeName '_' timestamp '.png']));
end

function plot_distance_comparison(scenario, labels, outputDir, timestamp, safeName)
    [~, order] = sort(scenario.gps_exposure_score, 'descend');
    topK = min(20, numel(order));
    topOrder = order(1:topK);

    figure('Name', ['Two-symbol distance - ' scenario.name], 'Color', 'w');
    bar(scenario.two_symbol_distance(topOrder, :));
    grid on;
    xlabel('Ranked channel case');
    ylabel('Two-symbol nearest distance');
    title(['Two-symbol distance comparison: ' scenario.name], 'Interpreter', 'none');
    legend(labels, 'Location', 'best');
    xticks(1:topK);
    xticklabels(short_case_names(scenario.cases(topOrder)));
    xtickangle(45);
    saveas(gcf, fullfile(outputDir, ['c2_two_symbol_distance_' safeName '_' timestamp '.png']));
end

function plot_distance_ratio_histogram(scenario, outputDir, timestamp, safeName)
    gpsDistance = scenario.two_symbol_distance(:, 2);
    proposedDistance = scenario.two_symbol_distance(:, 3);
    ratioDb = 10 * log10((gpsDistance + eps) ./ (proposedDistance + eps));

    figure('Name', ['GPS/Proposed distance ratio - ' scenario.name], 'Color', 'w');
    histogram(ratioDb, 30);
    grid on;
    xlabel('10log10(d2 GPS / d2 Proposed) (dB)');
    ylabel('Channel cases');
    title(['GPS vs proposed distance ratio: ' scenario.name], 'Interpreter', 'none');
    xline(0, '--', 'Equal');
    saveas(gcf, fullfile(outputDir, ['c2_distance_ratio_' safeName '_' timestamp '.png']));
end

function labels = short_case_names(cases)
    labels = cell(1, numel(cases));
    for idx = 1:numel(cases)
        labels{idx} = cases{idx}.name;
        if strlength(labels{idx}) > 24
            labels{idx} = extractBefore(labels{idx}, 25);
        end
    end
end

function metrics = finite_constellation_distance_metrics(H_eff, M_mod, modType)
    nearest_deltas = nearest_constellation_deltas(M_mod, modType);
    min_delta_power = min(abs(nearest_deltas).^2);

    column_norm_power = sum(abs(H_eff).^2, 1);
    metrics.single_symbol_distance = min(column_norm_power) * min_delta_power;
    metrics.two_symbol_distance = two_symbol_min_distance(H_eff, nearest_deltas);
end

function deltas = nearest_constellation_deltas(M_mod, modType)
    switch lower(modType)
        case 'qam'
            constellation = qammod((0:M_mod-1).', M_mod, 'UnitAveragePower', true);
        case 'psk'
            constellation = pskmod((0:M_mod-1).', M_mod, pi/M_mod);
        otherwise
            error('Unsupported modulation type: %s', modType);
    end

    diff_matrix = constellation - constellation.';
    all_deltas = diff_matrix(abs(diff_matrix) > 0);
    min_abs = min(abs(all_deltas));
    deltas = all_deltas(abs(abs(all_deltas) - min_abs) < 1e-10);
    deltas = unique(round(deltas * 1e12) / 1e12);
end

function min_distance = two_symbol_min_distance(H_eff, deltas)
    num_columns = size(H_eff, 2);
    min_distance = inf;

    for first_col = 1:num_columns-1
        h_first = H_eff(:, first_col);
        for second_col = first_col+1:num_columns
            h_second = H_eff(:, second_col);
            for first_delta_idx = 1:numel(deltas)
                first_term = h_first * deltas(first_delta_idx);
                for second_delta_idx = 1:numel(deltas)
                    diff_vec = first_term + h_second * deltas(second_delta_idx);
                    distance = sum(abs(diff_vec).^2);
                    if distance < min_distance
                        min_distance = distance;
                    end
                end
            end
        end
    end
end
