function results = run_caseA_projection_floor_audit(options)
%RUN_CASEA_PROJECTION_FLOOR_AUDIT Explain the fixed-channel MMSE error floor.

    if nargin < 1
        options = struct();
    end

    rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    setup_paths(rootDir);

    N = get_option(options, 'N', 64);
    V = get_option(options, 'V', 4);
    alphaMax = get_option(options, 'alpha_max', 3);
    c1 = get_option(options, 'c1', (2 * alphaMax + 1) / (2 * N));
    c2Base = get_option(options, 'c2_base', sqrt(2) / (10 * N));
    delta = get_option(options, 'proposed_delta', c2Base / 16);
    pattern = get_option(options, 'pattern', [2 2 1 1]);
    theta = get_option(options, 'theta', pi);
    numFrames = get_option(options, 'num_frames', 10000);
    seed = get_option(options, 'seed', 20260907);
    delays = get_option(options, 'delays', [0, 2]);
    dopplers = get_option(options, 'dopplers', [0, 2]);

    c2Gps = afdm.chirp.build_gps_pattern(N, V, pattern);
    c2Prop = afdm.chirp.build_proposed_pattern( ...
        N, V, pattern, c2Base, delta);
    schemes = ["baseline", "GPS", "proposed"];
    c2Values = {c2Base, c2Gps, c2Prop};
    gains = [1, exp(1i * theta)] / sqrt(2);

    rows = repmat(empty_row(), numel(schemes), 1);
    details = cell(numel(schemes), 1);
    for schemeIdx = 1:numel(schemes)
        hEff = zeros(N, N);
        for pathIdx = 1:numel(delays)
            hEff = hEff + gains(pathIdx) * afdm.analysis.build_path_matrix( ...
                N, c1, c2Values{schemeIdx}, ...
                delays(pathIdx), dopplers(pathIdx));
        end
        metrics = afdm.analysis.linear_detector_floor_metrics( ...
            hEff, 2, 'psk', numFrames, seed);
        details{schemeIdx} = metrics;

        row = empty_row();
        row.scheme = schemes(schemeIdx);
        row.rank_h = metrics.rank_h;
        row.min_singular_h = metrics.min_singular_h;
        row.projector_error_fro = metrics.projector_error_fro;
        row.min_row_sir_db = metrics.min_row_sir_db;
        row.median_row_sir_db = metrics.median_row_sir_db;
        row.noiseless_ber = metrics.noiseless_ber;
        row.min_decision_margin = metrics.min_decision_margin;
        row.p001_decision_margin = metrics.p001_decision_margin;
        row.nonpositive_margin_fraction = metrics.nonpositive_margin_fraction;
        row.bit_errors = metrics.noiseless_bit_errors;
        row.total_bits = metrics.total_bits;
        rows(schemeIdx) = row;
    end

    summary = struct2table(rows);
    results.config = struct( ...
        'N', N, 'V', V, 'alpha_max', alphaMax, 'c1', c1, ...
        'c2_base', c2Base, 'delta', delta, 'pattern', pattern, ...
        'theta', theta, 'delays', delays, 'dopplers', dopplers, ...
        'num_frames', numFrames, 'seed', seed);
    results.summary = summary;
    results.details = details;

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    save(fullfile(outputDir, 'caseA_projection_floor_audit.mat'), 'results');
    writetable(summary, fullfile(outputDir, 'caseA_projection_floor_audit.csv'));
    disp(summary);
end

function row = empty_row()
    row = struct( ...
        'scheme', "", ...
        'rank_h', NaN, ...
        'min_singular_h', NaN, ...
        'projector_error_fro', NaN, ...
        'min_row_sir_db', NaN, ...
        'median_row_sir_db', NaN, ...
        'noiseless_ber', NaN, ...
        'min_decision_margin', NaN, ...
        'p001_decision_margin', NaN, ...
        'nonpositive_margin_fraction', NaN, ...
        'bit_errors', NaN, ...
        'total_bits', NaN);
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = defaultValue;
    end
end
