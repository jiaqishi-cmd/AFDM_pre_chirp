function results = run_gps_variant_comparison(options)
%RUN_GPS_VARIANT_COMPARISON Compare GPS number models and quantization.

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
    delayList = get_option(options, 'delay_list', 1:8);
    dopplerList = get_option(options, 'doppler_list', -alphaMax:alphaMax);
    phaseBitsList = get_option(options, 'phase_bits_list', [4 6 8 10 12]);
    precisionK = get_option(options, 'precision_k', 2);
    numPaprFrames = get_option(options, 'num_papr_frames', 2000);
    numMmseFrames = get_option(options, 'num_mmse_frames', 5000);
    mmseSnrDb = get_option(options, 'mmse_snr_db', 28);
    oversamplingFactor = get_option(options, 'oversampling_factor', 4);
    seed = get_option(options, 'seed', 20260908);
    differenceAlphabet = [-2i; 2i];

    [variants, c2Patterns, candidateSets, phaseErrors] = build_variants( ...
        N, V, pattern, c2Base, delta, phaseBitsList, precisionK);
    numVariants = numel(variants);
    groupIndex = afdm.chirp.group_index(N, V, 'contiguous');

    structuralRows = repmat(empty_structural_row(), numVariants, 1);
    for variantIdx = 1:numVariants
        h1 = afdm.analysis.build_path_matrix( ...
            N, c1, c2Patterns{variantIdx}, 0, 0);
        riskySupports = 0;
        minMismatch = Inf;
        minClosestSigma = Inf;

        for delay2 = delayList
            for doppler2 = dopplerList
                h2 = afdm.analysis.build_path_matrix( ...
                    N, c1, c2Patterns{variantIdx}, delay2, doppler2);
                relativeOperator = h1' * h2;
                cycleAnalysis = afdm.analysis.monomial_cycle_analysis( ...
                    relativeOperator, differenceAlphabet);
                riskySupports = riskySupports + ...
                    (cycleAnalysis.num_compatible_eigenvectors > 0);
                minMismatch = min(minMismatch, ...
                    cycleAnalysis.min_difference_alphabet_mismatch);

                closestDelta = cycleAnalysis.closest_delta;
                psi = afdm.analysis.build_pairwise_path_response( ...
                    {h1, h2}, closestDelta / norm(closestDelta));
                metrics = afdm.analysis.pairwise_diversity_metrics( ...
                    psi, 1e-8);
                minClosestSigma = min(minClosestSigma, ...
                    metrics.min_singular_value);
            end
        end

        row = empty_structural_row();
        row.variant = variants(variantIdx);
        row.phase_error_rad = phaseErrors(variantIdx);
        row.risky_supports = riskySupports;
        row.total_supports = numel(delayList) * numel(dopplerList);
        row.min_alphabet_mismatch = minMismatch;
        row.min_closest_sigma = minClosestSigma;
        row.num_candidates = size(candidateSets{variantIdx}, 2);
        structuralRows(variantIdx) = row;
    end

    paprSamples = run_papr_trials(candidateSets, variants, groupIndex, ...
        c1, numPaprFrames, oversamplingFactor, seed);
    structuralSummary = struct2table(structuralRows);
    structuralSummary.mean_papr = mean(paprSamples, 1).';
    structuralSummary.p90_papr = percentile_by_sort(paprSamples, 0.90).';
    structuralSummary.p99_papr = percentile_by_sort(paprSamples, 0.99).';
    [mmseBer, zeroNoiseBer] = run_fixed_mmse_trials( ...
        c2Patterns, N, c1, numMmseFrames, mmseSnrDb, seed + 100000);
    structuralSummary.mmse_ber = mmseBer;
    structuralSummary.zero_noise_ber = zeroNoiseBer;

    results.config = struct( ...
        'N', N, 'V', V, 'alpha_max', alphaMax, 'c1', c1, ...
        'c2_base', c2Base, 'delta', delta, 'pattern', pattern, ...
        'delay_list', delayList, 'doppler_list', dopplerList, ...
        'phase_bits_list', phaseBitsList, 'precision_k', precisionK, ...
        'num_papr_frames', numPaprFrames, ...
        'num_mmse_frames', numMmseFrames, 'mmse_snr_db', mmseSnrDb, ...
        'oversampling_factor', oversamplingFactor, 'seed', seed);
    results.summary = structuralSummary;
    results.papr_samples = paprSamples;
    results.c2_patterns = c2Patterns;
    results.candidate_sets = candidateSets;

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    outputTag = sprintf('gps_variant_comparison_k%d_snr%d', ...
        precisionK, round(mmseSnrDb));
    save(fullfile(outputDir, [outputTag '.mat']), 'results');
    writetable(structuralSummary, ...
        fullfile(outputDir, [outputTag '.csv']));
    disp(structuralSummary);
end

function [ber, zeroNoiseBer] = run_fixed_mmse_trials( ...
        c2Patterns, N, c1, numFrames, snrDb, seed)
    numVariants = numel(c2Patterns);
    rng(seed, 'twister');
    bits = randi([0, 1], N, numFrames);
    symbols = 1i * (1 - 2 * bits);
    unitNoise = (randn(N, numFrames) + 1i * randn(N, numFrames)) / sqrt(2);
    snrLinear = 10^(snrDb / 10);
    ber = zeros(numVariants, 1);
    zeroNoiseBer = zeros(numVariants, 1);

    for variantIdx = 1:numVariants
        h1 = afdm.analysis.build_path_matrix( ...
            N, c1, c2Patterns{variantIdx}, 0, 0);
        h2 = afdm.analysis.build_path_matrix( ...
            N, c1, c2Patterns{variantIdx}, 2, 2);
        hEff = (h1 - h2) / sqrt(2);
        cleanRx = hEff * symbols;
        noiseVar = mean(abs(cleanRx).^2, 'all') / snrLinear;
        equalizer = (hEff' * hEff + noiseVar * eye(N)) \ hEff';
        estimate = equalizer * (cleanRx + sqrt(noiseVar) * unitNoise);
        decisions = pskdemod(estimate(:), 2, pi / 2, 'OutputType', 'bit');
        ber(variantIdx) = mean(decisions ~= bits(:));

        zeroEstimate = pinv(hEff) * cleanRx;
        zeroDecisions = pskdemod( ...
            zeroEstimate(:), 2, pi / 2, 'OutputType', 'bit');
        zeroNoiseBer(variantIdx) = mean(zeroDecisions ~= bits(:));
    end
end

function paprSamples = run_papr_trials(candidateSets, variants, groupIndex, ...
        c1, numFrames, oversamplingFactor, seed)
    numVariants = numel(variants);
    N = numel(groupIndex);
    paprSamples = zeros(numFrames, numVariants);

    oldRngState = rng;
    cleanup = onCleanup(@() rng(oldRngState)); %#ok<NASGU>
    for frameIdx = 1:numFrames
        rng(seed + frameIdx, 'twister');
        [symbols, ~] = afdm.tx.random_data(N, 16, 'qam');
        for variantIdx = 1:numVariants
            candidateSet = candidateSets{variantIdx};
            if size(candidateSet, 2) == 1
                waveform = afdm.search.full_waveform( ...
                    symbols, candidateSet(:, 1), oversamplingFactor);
                paprSamples(frameIdx, variantIdx) = ...
                    afdm.tx.compute_papr(waveform);
            else
                selection = afdm.search.greedy_group_papr_selection( ...
                    symbols, N, c1, candidateSet, groupIndex, ...
                    oversamplingFactor);
                paprSamples(frameIdx, variantIdx) = selection.papr;
            end
        end
    end
end

function [variants, c2Patterns, candidateSets, phaseErrors] = ...
        build_variants(N, V, pattern, c2Base, delta, phaseBitsList, precisionK)
    variants = ["baseline", "GPS-rational", "GPS-irrational"];
    for bits = phaseBitsList
        variants(end + 1) = "GPS-q" + string(bits); %#ok<AGROW>
    end
    variants(end + 1) = "proposed-2";
    variants(end + 1) = "proposed-3";

    numVariants = numel(variants);
    c2Patterns = cell(numVariants, 1);
    candidateSets = cell(numVariants, 1);
    phaseErrors = NaN(numVariants, 1);

    c2Patterns{1} = c2Base;
    candidateSets{1} = c2Base * ones(N, 1);

    options.precision_k = precisionK;
    [candidateSets{2}, metadata] = afdm.chirp.gps_candidate_set_variant( ...
        N, 'rational', options);
    c2Patterns{2} = select_pattern(candidateSets{2}, N, V, pattern);
    phaseErrors(2) = metadata.phase_error_from_pi_over_2;

    [candidateSets{3}, metadata] = afdm.chirp.gps_candidate_set_variant( ...
        N, 'irrational', options);
    c2Patterns{3} = select_pattern(candidateSets{3}, N, V, pattern);
    phaseErrors(3) = metadata.phase_error_from_pi_over_2;

    for bitsIdx = 1:numel(phaseBitsList)
        variantIdx = 3 + bitsIdx;
        options.phase_bits = phaseBitsList(bitsIdx);
        [candidateSets{variantIdx}, metadata] = ...
            afdm.chirp.gps_candidate_set_variant( ...
            N, 'phase_quantized', options);
        c2Patterns{variantIdx} = select_pattern( ...
            candidateSets{variantIdx}, N, V, pattern);
        phaseErrors(variantIdx) = metadata.phase_error_from_pi_over_2;
    end

    proposed2Idx = numVariants - 1;
    candidateSets{proposed2Idx} = c2Base + ...
        repmat([-delta, delta], N, 1);
    c2Patterns{proposed2Idx} = select_pattern( ...
        candidateSets{proposed2Idx}, N, V, pattern);

    proposedIdx = numVariants;
    proposedOffsets = [0, -delta, delta];
    candidateSets{proposedIdx} = c2Base + ...
        repmat(proposedOffsets, N, 1);
    c2Patterns{proposedIdx} = select_pattern( ...
        candidateSets{proposedIdx}, N, V, pattern);
end

function c2 = select_pattern(candidateSet, N, V, pattern)
    groupIndex = afdm.chirp.group_index(N, V, 'contiguous');
    c2 = zeros(N, 1);
    for groupIdx = 1:V
        selected = groupIndex == groupIdx;
        c2(selected) = candidateSet(selected, pattern(groupIdx));
    end
end

function values = percentile_by_sort(samples, probability)
    sorted = sort(samples, 1);
    index = max(1, min(size(sorted, 1), ...
        ceil(probability * size(sorted, 1))));
    values = sorted(index, :);
end

function row = empty_structural_row()
    row = struct( ...
        'variant', "", ...
        'phase_error_rad', NaN, ...
        'risky_supports', NaN, ...
        'total_supports', NaN, ...
        'min_alphabet_mismatch', NaN, ...
        'min_closest_sigma', NaN, ...
        'num_candidates', NaN);
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = defaultValue;
    end
end
