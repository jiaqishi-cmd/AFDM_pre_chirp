function results = run_prechirp_search_runtime_comparison(numRepeats, options)
%RUN_PRECHIRP_SEARCH_RUNTIME_COMPARISON Compare proposed c2 search runtime.
%
% Patent support point:
% This experiment validates that partial waveform reuse reduces the number of
% full IDAFT waveform constructions from K^V in naive exhaustive search to
% V*K precomputed partial waveforms, and that group-wise greedy search reduces
% candidate evaluations to V*K. It assumes equal-length groups, which is the
% intended practical grouping for the proposed pre-chirp design. The proposed
% pre-chirp c2 structure therefore admits lower-complexity runtime
% implementations for PAPR-aware c2 selection.

    rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    setup_paths(rootDir);

    if nargin < 1 || isempty(numRepeats)
        numRepeats = 10;
    end
    if nargin < 2
        options = struct();
    end

    N = get_option(options, 'N', 64);
    K = get_option(options, 'K', 3);
    VList = get_option(options, 'V_list', [2, 4, 8]);
    rho = get_option(options, 'rho', 0.2);
    os = get_option(options, 'oversampling', 4);
    seed = get_option(options, 'seed', 20260513);

    cfg = afdm_config();
    baseC2 = get_option(options, 'base_c2', sqrt(2) / (10 * N));
    offsets = get_option(options, 'offsets', linspace(-rho, rho, K) * baseC2);
    offsets = offsets(:).';

    numV = numel(VList);
    runtimeNaive = zeros(numV, 1);
    runtimePartialExhaustive = zeros(numV, 1);
    runtimeGreedy = zeros(numV, 1);
    naiveCandidates = zeros(numV, 1);
    partialPrecomputeCandidates = zeros(numV, 1);
    greedyCandidates = zeros(numV, 1);
    bestPaprNaive = zeros(numV, 1);
    bestPaprPartialExhaustive = zeros(numV, 1);
    bestPaprGreedy = zeros(numV, 1);

    fprintf('========== proposed pre-chirp search runtime comparison ==========\n');
    fprintf('N=%d, K=%d, rho=%.3g, repeats=%d, oversampling=%d\n', N, K, rho, numRepeats, os);

    oldRngState = rng;
    cleanup = onCleanup(@() rng(oldRngState)); %#ok<NASGU>

    for vIdx = 1:numV
        V = VList(vIdx);
        if mod(N, V) ~= 0
            error('V=%d does not divide N=%d. Use equal-length groups for this experiment.', V, N);
        end
        groupIndex = repelem((1:V).', N / V);
        patterns = enumerate_patterns(K, V);

        naiveCandidates(vIdx) = K^V;
        partialPrecomputeCandidates(vIdx) = V * K;
        greedyCandidates(vIdx) = V * K;

        tNaive = zeros(numRepeats, 1);
        tPartial = zeros(numRepeats, 1);
        tGreedy = zeros(numRepeats, 1);
        paprNaive = zeros(numRepeats, 1);
        paprPartial = zeros(numRepeats, 1);
        paprGreedy = zeros(numRepeats, 1);

        fprintf('\nV=%d | naive candidates=%d | partial precompute=%d | greedy eval=%d\n', ...
            V, naiveCandidates(vIdx), partialPrecomputeCandidates(vIdx), greedyCandidates(vIdx));

        for repeatIdx = 1:numRepeats
            rng(seed + 1000 * V + repeatIdx, 'twister');
            symbols = random_qam_symbols(N, cfg.modulation.M_mod);

            tic;
            paprNaive(repeatIdx) = naive_exhaustive_search(symbols, baseC2, offsets, groupIndex, patterns, os);
            tNaive(repeatIdx) = toc;

            tic;
            paprPartial(repeatIdx) = partial_reuse_exhaustive_search(symbols, baseC2, offsets, groupIndex, patterns, os);
            tPartial(repeatIdx) = toc;

            tic;
            paprGreedy(repeatIdx) = partial_reuse_greedy_search(symbols, baseC2, offsets, groupIndex, os);
            tGreedy(repeatIdx) = toc;
        end

        runtimeNaive(vIdx) = mean(tNaive);
        runtimePartialExhaustive(vIdx) = mean(tPartial);
        runtimeGreedy(vIdx) = mean(tGreedy);
        bestPaprNaive(vIdx) = mean(paprNaive);
        bestPaprPartialExhaustive(vIdx) = mean(paprPartial);
        bestPaprGreedy(vIdx) = mean(paprGreedy);

        fprintf('  runtime naive/partial/greedy = %.4g / %.4g / %.4g s\n', ...
            runtimeNaive(vIdx), runtimePartialExhaustive(vIdx), runtimeGreedy(vIdx));
    end

    results = table( ...
        VList(:), ...
        repmat(K, numV, 1), ...
        naiveCandidates, ...
        partialPrecomputeCandidates, ...
        greedyCandidates, ...
        runtimeNaive, ...
        runtimePartialExhaustive, ...
        runtimeGreedy, ...
        bestPaprNaive, ...
        bestPaprPartialExhaustive, ...
        bestPaprGreedy, ...
        'VariableNames', {'V', 'K', 'naive_candidates', 'partial_precompute_candidates', ...
        'greedy_candidates', 'runtime_naive', 'runtime_partial_exhaustive', ...
        'runtime_greedy', 'best_papr_naive', 'best_papr_partial_exhaustive', ...
        'best_papr_greedy'});

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    csvPath = fullfile(outputDir, ['prechirp_search_runtime_comparison_' timestamp '.csv']);
    pngPath = fullfile(outputDir, ['prechirp_search_runtime_comparison_' timestamp '.png']);
    matPath = fullfile(outputDir, ['prechirp_search_runtime_comparison_' timestamp '.mat']);

    writetable(results, csvPath);
    save(matPath, 'results', 'offsets', 'N', 'K', 'rho', 'numRepeats', 'os');
    plot_results(results, pngPath);

    fprintf('\nSaved runtime comparison to %s\n', csvPath);
    disp(results);
end

function bestPapr = naive_exhaustive_search(symbols, baseC2, offsets, groupIndex, patterns, os)
    bestPapr = inf;
    for idx = 1:size(patterns, 1)
        signal = afdm.search.direct_full_waveform(symbols, baseC2, offsets, groupIndex, patterns(idx, :), os);
        papr = afdm.tx.compute_papr(signal);
        if papr < bestPapr
            bestPapr = papr;
        end
    end
end

function bestPapr = partial_reuse_exhaustive_search(symbols, baseC2, offsets, groupIndex, patterns, os)
    sPart = afdm.search.precompute_partial_waveforms(symbols, baseC2, offsets, groupIndex, os);
    bestPapr = inf;
    for idx = 1:size(patterns, 1)
        signal = afdm.search.combine_partial_waveform(sPart, patterns(idx, :));
        papr = afdm.tx.compute_papr(signal);
        if papr < bestPapr
            bestPapr = papr;
        end
    end
end

function bestPapr = partial_reuse_greedy_search(symbols, baseC2, offsets, groupIndex, os)
    V = max(groupIndex);
    K = numel(offsets);
    sPart = afdm.search.precompute_partial_waveforms(symbols, baseC2, offsets, groupIndex, os);
    pattern = ones(1, V);

    for groupId = 1:V
        groupBestPapr = inf;
        groupBestCandidate = pattern(groupId);
        for candidateId = 1:K
            trialPattern = pattern;
            trialPattern(groupId) = candidateId;
            signal = afdm.search.combine_partial_waveform(sPart, trialPattern);
            papr = afdm.tx.compute_papr(signal);
            if papr < groupBestPapr
                groupBestPapr = papr;
                groupBestCandidate = candidateId;
            end
        end
        pattern(groupId) = groupBestCandidate;
    end

    signal = afdm.search.combine_partial_waveform(sPart, pattern);
    bestPapr = afdm.tx.compute_papr(signal);
end

function patterns = enumerate_patterns(K, V)
    numPatterns = K^V;
    patterns = ones(numPatterns, V);
    for idx = 0:numPatterns-1
        value = idx;
        for groupId = 1:V
            patterns(idx + 1, groupId) = mod(value, K) + 1;
            value = floor(value / K);
        end
    end
end

function symbols = random_qam_symbols(N, MMod)
    bits = randi([0, 1], N * log2(MMod), 1);
    symbols = qammod(bits, MMod, 'InputType', 'bit', 'UnitAveragePower', true);
end

function plot_results(results, pngPath)
    fig = figure('Name', 'proposed pre-chirp search runtime comparison', 'Color', 'w');
    tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    semilogy(results.V, results.runtime_naive, 'o-', 'LineWidth', 2);
    hold on;
    semilogy(results.V, results.runtime_partial_exhaustive, 's-', 'LineWidth', 2);
    semilogy(results.V, results.runtime_greedy, '^-', 'LineWidth', 2);
    grid on;
    xlabel('Number of groups V');
    ylabel('Average runtime (s)');
    legend({'naive exhaustive', 'partial reuse exhaustive', 'partial reuse greedy'}, ...
        'Location', 'northwest');
    title('Runtime vs V');

    nexttile;
    semilogy(results.V, results.naive_candidates, 'o-', 'LineWidth', 2);
    hold on;
    semilogy(results.V, results.partial_precompute_candidates, 's-', 'LineWidth', 2);
    semilogy(results.V, results.greedy_candidates, '^-', 'LineWidth', 2);
    grid on;
    xlabel('Number of groups V');
    ylabel('Candidate count');
    legend({'K^V full waveforms', 'V*K partial precompute', 'V*K greedy eval'}, ...
        'Location', 'northwest');
    title('Theoretical complexity');

    exportgraphics(fig, pngPath, 'Resolution', 200);
    close(fig);
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = defaultValue;
    end
end
