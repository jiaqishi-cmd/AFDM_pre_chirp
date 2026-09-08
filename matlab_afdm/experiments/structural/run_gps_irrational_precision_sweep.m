function results = run_gps_irrational_precision_sweep(options)
%RUN_GPS_IRRATIONAL_PRECISION_SWEEP Sweep the GPS irrational precision k.

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
    pattern = get_option(options, 'pattern', [2 2 1 1]);
    kList = get_option(options, 'k_list', 1:8);
    delayList = get_option(options, 'delay_list', 1:8);
    dopplerList = get_option(options, 'doppler_list', -alphaMax:alphaMax);
    numPaprFrames = get_option(options, 'num_papr_frames', 500);
    oversamplingFactor = get_option(options, 'oversampling_factor', 4);
    seed = get_option(options, 'seed', 20260908);
    differenceAlphabet = [-2i; 2i];
    groupIndex = afdm.chirp.group_index(N, V, 'contiguous');

    rows = repmat(empty_row(), numel(kList), 1);
    candidateSets = cell(numel(kList), 1);
    c2Patterns = cell(numel(kList), 1);

    for kIdx = 1:numel(kList)
        gpsOptions.precision_k = kList(kIdx);
        [candidateSet, metadata] = afdm.chirp.gps_candidate_set_variant( ...
            N, 'irrational', gpsOptions);
        c2Pattern = select_pattern(candidateSet, groupIndex, pattern);
        candidateSets{kIdx} = candidateSet;
        c2Patterns{kIdx} = c2Pattern;

        h1 = afdm.analysis.build_path_matrix(N, c1, c2Pattern, 0, 0);
        minMismatch = Inf;
        minSigma = Inf;
        riskySupports = 0;
        for delay2 = delayList
            for doppler2 = dopplerList
                h2 = afdm.analysis.build_path_matrix( ...
                    N, c1, c2Pattern, delay2, doppler2);
                cycleAnalysis = afdm.analysis.monomial_cycle_analysis( ...
                    h1' * h2, differenceAlphabet, 1e-12);
                riskySupports = riskySupports + ...
                    (cycleAnalysis.num_compatible_eigenvectors > 0);
                minMismatch = min(minMismatch, ...
                    cycleAnalysis.min_difference_alphabet_mismatch);
                delta = cycleAnalysis.closest_delta;
                psi = afdm.analysis.build_pairwise_path_response( ...
                    {h1, h2}, delta / norm(delta));
                metrics = afdm.analysis.pairwise_diversity_metrics(psi, 1e-12);
                minSigma = min(minSigma, metrics.min_singular_value);
            end
        end

        row = empty_row();
        row.precision_k = kList(kIdx);
        row.irrational_scale_minus_one = metadata.irrational_scale - 1;
        row.phase_error_rad = metadata.phase_error_from_pi_over_2;
        row.risky_supports = riskySupports;
        row.min_alphabet_mismatch = minMismatch;
        row.min_closest_sigma = minSigma;
        rows(kIdx) = row;
    end

    paprSamples = zeros(numPaprFrames, numel(kList));
    oldRngState = rng;
    cleanup = onCleanup(@() rng(oldRngState)); %#ok<NASGU>
    for frameIdx = 1:numPaprFrames
        rng(seed + frameIdx, 'twister');
        [symbols, ~] = afdm.tx.random_data(N, 16, 'qam');
        for kIdx = 1:numel(kList)
            selection = afdm.search.greedy_group_papr_selection( ...
                symbols, N, c1, candidateSets{kIdx}, groupIndex, ...
                oversamplingFactor);
            paprSamples(frameIdx, kIdx) = selection.papr;
        end
    end

    summary = struct2table(rows);
    summary.mean_papr = mean(paprSamples, 1).';
    summary.p99_papr = percentile_by_sort(paprSamples, 0.99).';

    results.config = struct( ...
        'N', N, 'V', V, 'alpha_max', alphaMax, 'c1', c1, ...
        'pattern', pattern, 'k_list', kList, ...
        'delay_list', delayList, 'doppler_list', dopplerList, ...
        'num_papr_frames', numPaprFrames, ...
        'oversampling_factor', oversamplingFactor);
    results.summary = summary;
    results.papr_samples = paprSamples;
    results.c2_patterns = c2Patterns;

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    save(fullfile(outputDir, 'gps_irrational_precision_sweep.mat'), 'results');
    writetable(summary, ...
        fullfile(outputDir, 'gps_irrational_precision_sweep.csv'));
    plot_results(summary, outputDir);
    disp(summary);
end

function c2 = select_pattern(candidateSet, groupIndex, pattern)
    c2 = zeros(size(groupIndex));
    for groupIdx = 1:numel(pattern)
        selected = groupIndex == groupIdx;
        c2(selected) = candidateSet(selected, pattern(groupIdx));
    end
end

function plot_results(summary, outputDir)
    figure('Name', 'GPS irrational precision safety', 'Color', 'w');
    semilogy(summary.precision_k, summary.min_closest_sigma, ...
        'o-', 'LineWidth', 2);
    grid on;
    xlabel('GPS irrational precision k');
    ylabel('Minimum closest-event singular value');
    title('Full rank does not imply a stable coding gain');
    saveas(gcf, fullfile(outputDir, 'gps_irrational_precision_safety.png'));
end

function values = percentile_by_sort(samples, probability)
    sorted = sort(samples, 1);
    index = max(1, min(size(sorted, 1), ...
        ceil(probability * size(sorted, 1))));
    values = sorted(index, :);
end

function row = empty_row()
    row = struct( ...
        'precision_k', NaN, ...
        'irrational_scale_minus_one', NaN, ...
        'phase_error_rad', NaN, ...
        'risky_supports', NaN, ...
        'min_alphabet_mismatch', NaN, ...
        'min_closest_sigma', NaN);
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = defaultValue;
    end
end
