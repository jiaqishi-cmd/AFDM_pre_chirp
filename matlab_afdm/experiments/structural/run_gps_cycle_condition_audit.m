function results = run_gps_cycle_condition_audit(options)
%RUN_GPS_CYCLE_CONDITION_AUDIT Explain GPS rank loss through path cycles.

    if nargin < 1
        options = struct();
    end

    rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    setup_paths(rootDir);

    configurations = get_option(options, 'configurations', default_configurations());
    pattern = get_option(options, 'pattern', [2 2 1 1]);
    differenceAlphabet = get_option(options, ...
        'difference_alphabet', [-2i; 2i]);

    schemes = ["baseline", "GPS", "proposed"];
    rows = repmat(empty_row(), numel(configurations) * numel(schemes), 1);
    details = cell(size(rows));
    cursor = 0;

    for configIdx = 1:numel(configurations)
        cfg = configurations(configIdx);
        c1 = (2 * cfg.alpha_max + 1) / (2 * cfg.N);
        c2Base = sqrt(2) / (10 * cfg.N);
        c2Values = { ...
            c2Base, ...
            afdm.chirp.build_gps_pattern(cfg.N, cfg.V, pattern), ...
            afdm.chirp.build_proposed_pattern( ...
                cfg.N, cfg.V, pattern, c2Base, c2Base / 16)};

        for schemeIdx = 1:numel(schemes)
            h1 = afdm.analysis.build_path_matrix( ...
                cfg.N, c1, c2Values{schemeIdx}, 0, 0);
            h2 = afdm.analysis.build_path_matrix( ...
                cfg.N, c1, c2Values{schemeIdx}, cfg.delay2, cfg.doppler2);
            relativeOperator = h1' * h2;
            cycleAnalysis = afdm.analysis.monomial_cycle_analysis( ...
                relativeOperator, differenceAlphabet);

            cursor = cursor + 1;
            details{cursor} = cycleAnalysis;
            row = empty_row();
            row.N = cfg.N;
            row.V = cfg.V;
            row.alpha_max = cfg.alpha_max;
            row.delay2 = cfg.delay2;
            row.doppler2 = cfg.doppler2;
            row.scheme = schemes(schemeIdx);
            row.effective_shift = cycleAnalysis.constant_shift;
            row.gcd_shift_N = gcd(round(cycleAnalysis.constant_shift), cfg.N);
            row.predicted_cycle_length = cfg.N / row.gcd_shift_N;
            row.observed_min_cycle_length = min( ...
                cycleAnalysis.cycle_table.cycle_length);
            row.compatible_eigenvectors = ...
                cycleAnalysis.num_compatible_eigenvectors;
            row.min_compatible_weight = cycleAnalysis.min_compatible_weight;
            row.monomial_leakage = ...
                cycleAnalysis.max_monomial_leakage_power;
            rows(cursor) = row;
        end
    end

    summary = struct2table(rows);
    results.summary = summary;
    results.details = details;
    results.configurations = configurations;
    results.pattern = pattern;
    results.difference_alphabet = differenceAlphabet;

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    save(fullfile(outputDir, 'gps_cycle_condition_audit.mat'), 'results');
    writetable(summary, fullfile(outputDir, 'gps_cycle_condition_audit.csv'));
    disp(summary);

    verify_compatible_events(results, schemes, configurations);
end

function verify_compatible_events(results, schemes, configurations)
    for rowIdx = 1:height(results.summary)
        if results.summary.compatible_eigenvectors(rowIdx) == 0
            continue;
        end
        cfg = configurations(ceil(rowIdx / numel(schemes)));
        c1 = (2 * cfg.alpha_max + 1) / (2 * cfg.N);
        c2Base = sqrt(2) / (10 * cfg.N);
        switch results.summary.scheme(rowIdx)
            case "baseline"
                c2 = c2Base;
            case "GPS"
                c2 = afdm.chirp.build_gps_pattern(cfg.N, cfg.V, results.pattern);
            case "proposed"
                c2 = afdm.chirp.build_proposed_pattern( ...
                    cfg.N, cfg.V, results.pattern, c2Base, c2Base / 16);
        end
        h1 = afdm.analysis.build_path_matrix(cfg.N, c1, c2, 0, 0);
        h2 = afdm.analysis.build_path_matrix( ...
            cfg.N, c1, c2, cfg.delay2, cfg.doppler2);
        delta = results.details{rowIdx}.compatible_deltas{1};
        psi = afdm.analysis.build_pairwise_path_response({h1, h2}, delta);
        metrics = afdm.analysis.pairwise_diversity_metrics( ...
            psi, 1e-8 * norm(psi, 2));
        if metrics.rank ~= 1
            error('Compatible cycle event did not produce rank-one Psi.');
        end
    end
end

function configurations = default_configurations()
    configurations = struct( ...
        'N', {8, 16, 64}, ...
        'V', {4, 4, 4}, ...
        'alpha_max', {2, 2, 3}, ...
        'delay2', {2, 2, 2}, ...
        'doppler2', {2, 2, 2});
end

function row = empty_row()
    row = struct( ...
        'N', NaN, ...
        'V', NaN, ...
        'alpha_max', NaN, ...
        'delay2', NaN, ...
        'doppler2', NaN, ...
        'scheme', "", ...
        'effective_shift', NaN, ...
        'gcd_shift_N', NaN, ...
        'predicted_cycle_length', NaN, ...
        'observed_min_cycle_length', NaN, ...
        'compatible_eigenvectors', NaN, ...
        'min_compatible_weight', NaN, ...
        'monomial_leakage', NaN);
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = defaultValue;
    end
end
