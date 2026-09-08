function results = run_gps_pep_bound_comparison(options)
%RUN_GPS_PEP_BOUND_COMPARISON Compare Rayleigh PEP for one danger event.

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
    localDelta = get_option(options, 'proposed_delta', c2Base / 16);
    pattern = get_option(options, 'pattern', [2 2 1 1]);
    delay2 = get_option(options, 'delay2', 2);
    doppler2 = get_option(options, 'doppler2', 2);
    snrDb = get_option(options, 'snr_db', 0:2:120);

    [variants, c2Patterns] = build_variants( ...
        N, V, pattern, c2Base, localDelta);
    numVariants = numel(variants);

    c2Rational = c2Patterns{2};
    h1Rational = afdm.analysis.build_path_matrix(N, c1, c2Rational, 0, 0);
    h2Rational = afdm.analysis.build_path_matrix( ...
        N, c1, c2Rational, delay2, doppler2);
    cycleAnalysis = afdm.analysis.monomial_cycle_analysis( ...
        h1Rational' * h2Rational, [-2i; 2i], 1e-12);
    if isempty(cycleAnalysis.compatible_deltas)
        error('No rational-GPS compatible error event found.');
    end
    deltaX = cycleAnalysis.compatible_deltas{1};
    deltaX = deltaX / norm(deltaX);

    rows = repmat(empty_row(), numVariants, 1);
    pepBound = zeros(numel(snrDb), numVariants);
    P = 2;
    snrLinear = 10.^(snrDb(:) / 10);

    for variantIdx = 1:numVariants
        h1 = afdm.analysis.build_path_matrix( ...
            N, c1, c2Patterns{variantIdx}, 0, 0);
        h2 = afdm.analysis.build_path_matrix( ...
            N, c1, c2Patterns{variantIdx}, delay2, doppler2);
        psi = afdm.analysis.build_pairwise_path_response({h1, h2}, deltaX);
        metrics = afdm.analysis.pairwise_diversity_metrics(psi, 1e-10);
        eigenvalues = metrics.gram_eigenvalues;
        activeEigenvalues = eigenvalues(eigenvalues > 1e-20);

        bound = ones(numel(snrLinear), 1);
        for eigenIdx = 1:numel(activeEigenvalues)
            bound = bound ./ (1 + ...
                activeEigenvalues(eigenIdx) * snrLinear / (4 * P));
        end
        pepBound(:, variantIdx) = bound;

        row = empty_row();
        row.variant = variants(variantIdx);
        row.rank_psi = metrics.rank;
        row.lambda_max = eigenvalues(1);
        row.lambda_min = eigenvalues(end);
        row.min_sigma = metrics.min_singular_value;
        row.gram_det = metrics.gram_determinant;
        row.snr_for_pep_1e_3 = first_snr_below(snrDb, bound, 1e-3);
        row.snr_for_pep_1e_5 = first_snr_below(snrDb, bound, 1e-5);
        rows(variantIdx) = row;
    end

    summary = struct2table(rows);
    results.config = struct( ...
        'N', N, 'V', V, 'alpha_max', alphaMax, 'c1', c1, ...
        'pattern', pattern, 'delay2', delay2, 'doppler2', doppler2);
    results.delta_x = deltaX;
    results.snr_db = snrDb;
    results.pep_bound = pepBound;
    results.summary = summary;

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    save(fullfile(outputDir, 'gps_pep_bound_comparison.mat'), 'results');
    writetable(summary, fullfile(outputDir, 'gps_pep_bound_comparison.csv'));
    plot_results(results, variants, outputDir);
    disp(summary);
end

function [variants, c2Patterns] = build_variants( ...
        N, V, pattern, c2Base, localDelta)
    variants = ["baseline", "GPS-rational", "GPS-irrational-k1", ...
        "GPS-irrational-k2", "GPS-q8-k2", "GPS-q12-k2", ...
        "proposed-2", "proposed-3"];
    c2Patterns = cell(numel(variants), 1);
    c2Patterns{1} = c2Base;
    c2Patterns{2} = afdm.chirp.build_gps_pattern_variant( ...
        N, V, pattern, 'rational');
    c2Patterns{3} = afdm.chirp.build_gps_pattern_variant( ...
        N, V, pattern, 'irrational', struct('precision_k', 1));
    c2Patterns{4} = afdm.chirp.build_gps_pattern_variant( ...
        N, V, pattern, 'irrational', struct('precision_k', 2));
    c2Patterns{5} = afdm.chirp.build_gps_pattern_variant( ...
        N, V, pattern, 'phase_quantized', ...
        struct('precision_k', 2, 'phase_bits', 8));
    c2Patterns{6} = afdm.chirp.build_gps_pattern_variant( ...
        N, V, pattern, 'phase_quantized', ...
        struct('precision_k', 2, 'phase_bits', 12));

    groupIndex = afdm.chirp.group_index(N, V, 'contiguous');
    candidate2 = c2Base + repmat([-localDelta, localDelta], N, 1);
    candidate3 = c2Base + repmat([0, -localDelta, localDelta], N, 1);
    c2Patterns{7} = select_pattern(candidate2, groupIndex, pattern);
    c2Patterns{8} = select_pattern(candidate3, groupIndex, pattern);
end

function c2 = select_pattern(candidateSet, groupIndex, pattern)
    c2 = zeros(size(groupIndex));
    for groupIdx = 1:numel(pattern)
        selected = groupIndex == groupIdx;
        c2(selected) = candidateSet(selected, pattern(groupIdx));
    end
end

function snrValue = first_snr_below(snrDb, values, threshold)
    index = find(values <= threshold, 1);
    if isempty(index)
        snrValue = NaN;
    else
        snrValue = snrDb(index);
    end
end

function plot_results(results, variants, outputDir)
    figure('Name', 'GPS danger-event PEP bound', 'Color', 'w');
    semilogy(results.snr_db, results.pep_bound, 'LineWidth', 1.8);
    grid on;
    xlabel('SNR (dB)');
    ylabel('Rayleigh PEP upper bound');
    title('PEP bound for the rational-GPS danger event');
    legend(cellstr(variants), 'Location', 'southwest');
    saveas(gcf, fullfile(outputDir, 'gps_pep_bound_comparison.png'));
end

function row = empty_row()
    row = struct( ...
        'variant', "", ...
        'rank_psi', NaN, ...
        'lambda_max', NaN, ...
        'lambda_min', NaN, ...
        'min_sigma', NaN, ...
        'gram_det', NaN, ...
        'snr_for_pep_1e_3', NaN, ...
        'snr_for_pep_1e_5', NaN);
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = defaultValue;
    end
end
