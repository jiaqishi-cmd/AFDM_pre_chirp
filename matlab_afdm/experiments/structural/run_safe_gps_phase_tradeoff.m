function results = run_safe_gps_phase_tradeoff(options)
%RUN_SAFE_GPS_PHASE_TRADEOFF Test offline-safe GPS pattern filtering.

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
    epsilonList = get_option(options, 'epsilon_list', ...
        [0 0.01 0.02 0.05 0.1 0.2 0.3 0.4]);
    delayList = get_option(options, 'delay_list', 1:8);
    dopplerList = get_option(options, 'doppler_list', -alphaMax:alphaMax);
    fixedDelay = get_option(options, 'fixed_delay', 2);
    fixedDoppler = get_option(options, 'fixed_doppler', 2);
    paprFrames = get_option(options, 'papr_frames', 2000);
    berFrames = get_option(options, 'ber_frames', 5000);
    snrDb = get_option(options, 'snr_db', 28);
    oversamplingFactor = get_option(options, 'oversampling_factor', 4);
    seed = get_option(options, 'seed', 20260908);

    groupIndex = afdm.chirp.group_index(N, V, 'contiguous');
    patterns = enumerate_patterns(2, V);
    qam16Delta = afdm.delta.generate_sparse_difference_set(1, 16, 'qam', 1);
    differenceAlphabet = qam16Delta.difference_alphabet;
    rows = repmat(empty_row(), numel(epsilonList), 1);
    details = cell(numel(epsilonList), 1);

    for epsilonIdx = 1:numel(epsilonList)
        epsilon = epsilonList(epsilonIdx);
        candidateSet = afdm.chirp.gps_candidate_set_variant( ...
            N, 'phase_offset', struct('phase_offset', epsilon));
        patternData = analyze_patterns( ...
            candidateSet, patterns, groupIndex, N, c1, ...
            delayList, dopplerList, fixedDelay, fixedDoppler, ...
            differenceAlphabet, snrDb);

        safePatternIdx = find(~patternData.unsafe);
        [originalPapr, safePapr, originalBer, safeBer, ...
            originalRiskProbability] = simulate_selection( ...
            candidateSet, patterns, patternData, safePatternIdx, ...
            groupIndex, N, c1, paprFrames, berFrames, snrDb, ...
            oversamplingFactor, seed);

        row = empty_row();
        row.phase_offset = epsilon;
        row.unsafe_patterns = sum(patternData.unsafe);
        row.safe_patterns = numel(safePatternIdx);
        if isempty(safePatternIdx)
            row.min_safe_sigma = NaN;
        else
            row.min_safe_sigma = min( ...
                patternData.min_sigma_over_support(safePatternIdx));
        end
        row.original_mean_papr = mean(originalPapr);
        row.safe_mean_papr = mean(safePapr, 'omitnan');
        row.original_p99_papr = percentile_by_sort(originalPapr, 0.99);
        row.safe_p99_papr = percentile_by_sort(safePapr, 0.99);
        row.original_mmse_ber = originalBer;
        row.safe_mmse_ber = safeBer;
        row.original_risk_probability = originalRiskProbability;
        row.original_evaluations = V + 1;
        row.safe_evaluations = numel(safePatternIdx);
        rows(epsilonIdx) = row;

        details{epsilonIdx}.candidate_set = candidateSet;
        details{epsilonIdx}.pattern_data = patternData;
        details{epsilonIdx}.safe_pattern_indices = safePatternIdx;
    end

    summary = struct2table(rows);
    localReference = evaluate_local_references( ...
        N, V, c1, c2Base, localDelta, groupIndex, delayList, dopplerList, ...
        fixedDelay, fixedDoppler, differenceAlphabet, ...
        paprFrames, berFrames, snrDb, oversamplingFactor, seed);
    results.config = struct( ...
        'N', N, 'V', V, 'alpha_max', alphaMax, 'c1', c1, ...
        'c2_base', c2Base, 'local_delta', localDelta, ...
        'epsilon_list', epsilonList, 'delay_list', delayList, ...
        'doppler_list', dopplerList, 'fixed_delay', fixedDelay, ...
        'fixed_doppler', fixedDoppler, 'papr_frames', paprFrames, ...
        'ber_frames', berFrames, 'snr_db', snrDb, ...
        'oversampling_factor', oversamplingFactor, 'seed', seed);
    results.summary = summary;
    results.local_reference = localReference;
    results.details = details;

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    save(fullfile(outputDir, 'safe_gps_phase_tradeoff.mat'), 'results');
    writetable(summary, fullfile(outputDir, 'safe_gps_phase_tradeoff.csv'));
    writetable(localReference, ...
        fullfile(outputDir, 'safe_gps_local_reference.csv'));
    plot_results(summary, outputDir);
    disp(summary);
    disp(localReference);
end

function summary = evaluate_local_references( ...
        N, V, c1, c2Base, localDelta, groupIndex, ...
        delayList, dopplerList, fixedDelay, fixedDoppler, ...
        differenceAlphabet, paprFrames, berFrames, snrDb, ...
        oversamplingFactor, seed)
    names = ["proposed-2", "proposed-3"];
    candidateSets = { ...
        c2Base + repmat([-localDelta, localDelta], N, 1), ...
        c2Base + repmat([0, -localDelta, localDelta], N, 1)};
    rows = repmat(struct( ...
        'scheme', "", 'num_candidates', NaN, 'unsafe_patterns', NaN, ...
        'mean_papr', NaN, 'p99_papr', NaN, 'mmse_ber', NaN, ...
        'risk_probability', NaN, 'online_evaluations', NaN), 2, 1);

    for schemeIdx = 1:2
        candidateSet = candidateSets{schemeIdx};
        numCandidates = size(candidateSet, 2);
        patterns = enumerate_patterns(numCandidates, V);
        patternData = analyze_patterns( ...
            candidateSet, patterns, groupIndex, N, c1, ...
            delayList, dopplerList, fixedDelay, fixedDoppler, ...
            differenceAlphabet, snrDb);
        [papr, ~, ber, ~, riskProbability] = simulate_selection( ...
            candidateSet, patterns, patternData, [], groupIndex, ...
            N, c1, paprFrames, berFrames, snrDb, ...
            oversamplingFactor, seed);
        row = rows(schemeIdx);
        row.scheme = names(schemeIdx);
        row.num_candidates = numCandidates;
        row.unsafe_patterns = sum(patternData.unsafe);
        row.mean_papr = mean(papr);
        row.p99_papr = percentile_by_sort(papr, 0.99);
        row.mmse_ber = ber;
        row.risk_probability = riskProbability;
        row.online_evaluations = 1 + V * (numCandidates - 1);
        rows(schemeIdx) = row;
    end
    summary = struct2table(rows);
end

function data = analyze_patterns( ...
        candidateSet, patterns, groupIndex, N, c1, ...
        delayList, dopplerList, fixedDelay, fixedDoppler, ...
        differenceAlphabet, snrDb)
    numPatterns = size(patterns, 1);
    unsafe = false(numPatterns, 1);
    minSigma = Inf(numPatterns, 1);
    hFixed = cell(numPatterns, 1);
    equalizer = cell(numPatterns, 1);
    noiseVar = 10^(-snrDb / 10);

    for patternIdx = 1:numPatterns
        c2 = select_c2(candidateSet, groupIndex, patterns(patternIdx, :));
        h1 = eye(N);
        for delay2 = delayList
            for doppler2 = dopplerList
                path = afdm.analysis.integer_path_monomial_parameters( ...
                    N, c1, c2, delay2, doppler2);
                h2 = build_monomial_matrix(path, N);
                cycle = afdm.analysis.monomial_cycle_analysis( ...
                    h1' * h2, differenceAlphabet, 1e-12);
                unsafe(patternIdx) = unsafe(patternIdx) || ...
                    cycle.num_compatible_eigenvectors > 0;
                delta = cycle.closest_delta;
                psi = afdm.analysis.build_pairwise_path_response( ...
                    {h1, h2}, delta / norm(delta));
                metrics = afdm.analysis.pairwise_diversity_metrics(psi, 1e-12);
                minSigma(patternIdx) = min( ...
                    minSigma(patternIdx), metrics.min_singular_value);
            end
        end

        fixedPath = afdm.analysis.integer_path_monomial_parameters( ...
            N, c1, c2, fixedDelay, fixedDoppler);
        h2Fixed = build_monomial_matrix(fixedPath, N);
        hFixed{patternIdx} = (h1 - h2Fixed) / sqrt(2);
        equalizer{patternIdx} = ...
            (hFixed{patternIdx}' * hFixed{patternIdx} + ...
            noiseVar * eye(N)) \ hFixed{patternIdx}';
    end

    data.unsafe = unsafe;
    data.min_sigma_over_support = minSigma;
    data.h_fixed = hFixed;
    data.equalizer = equalizer;
end

function matrix = build_monomial_matrix(path, N)
    matrix = zeros(N, N);
    linearIdx = sub2ind([N, N], path.permutation, (1:N).');
    matrix(linearIdx) = path.coefficient;
end

function [originalPapr, safePapr, originalBer, safeBer, riskProbability] = ...
        simulate_selection(candidateSet, patterns, patternData, safePatternIdx, ...
        groupIndex, N, c1, paprFrames, berFrames, snrDb, ...
        oversamplingFactor, seed)
    originalPapr = zeros(paprFrames, 1);
    safePapr = NaN(paprFrames, 1);
    rng(seed, 'twister');
    for frameIdx = 1:paprFrames
        [symbols, ~] = afdm.tx.random_data(N, 16, 'qam');
        selection = afdm.search.greedy_group_papr_selection( ...
            symbols, N, c1, candidateSet, groupIndex, oversamplingFactor);
        originalPapr(frameIdx) = selection.papr;
        if ~isempty(safePatternIdx)
            safePapr(frameIdx) = best_pattern_papr( ...
                symbols, candidateSet, patterns, safePatternIdx, ...
                groupIndex, oversamplingFactor);
        end
    end

    originalErrors = 0;
    safeErrors = 0;
    selectedUnsafe = 0;
    noiseVar = 10^(-snrDb / 10);
    rng(seed + 50000000, 'twister');
    for frameIdx = 1:berFrames
        bits = randi([0, 1], N, 1);
        symbols = 1i * (1 - 2 * bits);
        noise = sqrt(noiseVar / 2) * ...
            (randn(N, 1) + 1i * randn(N, 1));

        originalSelection = afdm.search.greedy_group_papr_selection( ...
            symbols, N, c1, candidateSet, groupIndex, oversamplingFactor);
        originalIdx = encode_pattern( ...
            originalSelection.selected_candidate_index, 2);
        originalEstimate = patternData.equalizer{originalIdx} * ...
            (patternData.h_fixed{originalIdx} * symbols + noise);
        originalBits = pskdemod( ...
            originalEstimate, 2, pi / 2, 'OutputType', 'bit');
        originalErrors = originalErrors + sum(originalBits ~= bits);
        selectedUnsafe = selectedUnsafe + patternData.unsafe(originalIdx);

        if ~isempty(safePatternIdx)
            safeIdx = best_pattern_index( ...
                symbols, candidateSet, patterns, safePatternIdx, ...
                groupIndex, oversamplingFactor);
            safeEstimate = patternData.equalizer{safeIdx} * ...
                (patternData.h_fixed{safeIdx} * symbols + noise);
            safeBits = pskdemod( ...
                safeEstimate, 2, pi / 2, 'OutputType', 'bit');
            safeErrors = safeErrors + sum(safeBits ~= bits);
        end
    end

    originalBer = originalErrors / (berFrames * N);
    if isempty(safePatternIdx)
        safeBer = NaN;
    else
        safeBer = safeErrors / (berFrames * N);
    end
    riskProbability = selectedUnsafe / berFrames;
end

function value = best_pattern_papr( ...
        symbols, candidateSet, patterns, allowedIdx, ...
        groupIndex, oversamplingFactor)
    [~, value] = best_pattern_index( ...
        symbols, candidateSet, patterns, allowedIdx, ...
        groupIndex, oversamplingFactor);
end

function [bestIdx, bestPapr] = best_pattern_index( ...
        symbols, candidateSet, patterns, allowedIdx, ...
        groupIndex, oversamplingFactor)
    bestPapr = Inf;
    bestIdx = allowedIdx(1);
    for listIdx = 1:numel(allowedIdx)
        patternIdx = allowedIdx(listIdx);
        c2 = select_c2(candidateSet, groupIndex, patterns(patternIdx, :));
        waveform = afdm.search.full_waveform( ...
            symbols, c2, oversamplingFactor);
        papr = afdm.tx.compute_papr(waveform);
        if papr < bestPapr
            bestPapr = papr;
            bestIdx = patternIdx;
        end
    end
end

function patternIdx = encode_pattern(pattern, numCandidates)
    powers = numCandidates .^ (0:numel(pattern)-1);
    patternIdx = 1 + sum((pattern(:).' - 1) .* powers);
end

function patterns = enumerate_patterns(numCandidates, numGroups)
    patterns = ones(numCandidates^numGroups, numGroups);
    for patternIdx = 0:size(patterns, 1)-1
        value = patternIdx;
        for groupIdx = 1:numGroups
            patterns(patternIdx + 1, groupIdx) = mod(value, numCandidates) + 1;
            value = floor(value / numCandidates);
        end
    end
end

function c2 = select_c2(candidateSet, groupIndex, pattern)
    c2 = zeros(size(groupIndex));
    for groupIdx = 1:numel(pattern)
        selected = groupIndex == groupIdx;
        c2(selected) = candidateSet(selected, pattern(groupIdx));
    end
end

function value = percentile_by_sort(samples, probability)
    samples = samples(~isnan(samples));
    if isempty(samples)
        value = NaN;
        return;
    end
    samples = sort(samples(:));
    index = max(1, min(numel(samples), ceil(probability * numel(samples))));
    value = samples(index);
end

function plot_results(summary, outputDir)
    figure('Name', 'Safe GPS phase tradeoff', 'Color', 'w');
    yyaxis left;
    plot(summary.phase_offset, summary.safe_mean_papr, ...
        'o-', 'LineWidth', 1.8);
    ylabel('Safe-selection mean PAPR (dB)');
    yyaxis right;
    semilogy(summary.phase_offset, summary.safe_mmse_ber, ...
        's-', 'LineWidth', 1.8);
    ylabel('Safe-selection MMSE BER');
    grid on;
    xlabel('GPS phase offset from pi/2 (rad)');
    title('Safety-filtered GPS PAPR and BER tradeoff');
    saveas(gcf, fullfile(outputDir, 'safe_gps_phase_tradeoff.png'));
end

function row = empty_row()
    row = struct( ...
        'phase_offset', NaN, ...
        'unsafe_patterns', NaN, ...
        'safe_patterns', NaN, ...
        'min_safe_sigma', NaN, ...
        'original_mean_papr', NaN, ...
        'safe_mean_papr', NaN, ...
        'original_p99_papr', NaN, ...
        'safe_p99_papr', NaN, ...
        'original_mmse_ber', NaN, ...
        'safe_mmse_ber', NaN, ...
        'original_risk_probability', NaN, ...
        'original_evaluations', NaN, ...
        'safe_evaluations', NaN);
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = defaultValue;
    end
end
