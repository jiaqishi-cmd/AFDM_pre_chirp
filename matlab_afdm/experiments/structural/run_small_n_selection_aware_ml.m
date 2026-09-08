function results = run_small_n_selection_aware_ml(options)
%RUN_SMALL_N_SELECTION_AWARE_ML Exact ML with data-driven PAPR patterns.

    if nargin < 1
        options = struct();
    end

    rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    setup_paths(rootDir);

    N = get_option(options, 'N', 8);
    V = get_option(options, 'V', 4);
    alphaMax = get_option(options, 'alpha_max', 2);
    c1 = get_option(options, 'c1', (2 * alphaMax + 1) / (2 * N));
    c2Base = get_option(options, 'c2_base', sqrt(2) / (10 * N));
    localDelta = get_option(options, 'proposed_delta', c2Base / 16);
    delay2 = get_option(options, 'delay2', 2);
    doppler2 = get_option(options, 'doppler2', 2);
    snrDb = get_option(options, 'snr_db', 12:4:28);
    numFrames = get_option(options, 'num_frames', 20000);
    oversamplingFactor = get_option(options, 'oversampling_factor', 4);
    seed = get_option(options, 'seed', 20260908);

    schemes = build_schemes(N, c2Base, localDelta);
    numSchemes = numel(schemes);
    groupIndex = afdm.chirp.group_index(N, V, 'contiguous');
    numCodewords = 2^N;
    codewordBits = de2bi(0:numCodewords-1, N, 'left-msb').';
    codewordSymbols = 1i * (1 - 2 * codewordBits);
    patternCache = build_pattern_cache( ...
        schemes, groupIndex, V, N, c1, delay2, doppler2);

    errorBits = zeros(numel(snrDb), numSchemes);
    selectedRisk = zeros(numel(snrDb), numSchemes);
    totalBits = numFrames * N * ones(numel(snrDb), numSchemes);
    oldRngState = rng;
    cleanup = onCleanup(@() rng(oldRngState)); %#ok<NASGU>

    for snrIdx = 1:numel(snrDb)
        rng(seed + 100000 * snrIdx, 'twister');
        noiseStd = sqrt(10^(-snrDb(snrIdx) / 10) / 2);
        txIndices = randi(numCodewords, numFrames, 1);
        gains = (randn(2, numFrames) + 1i * randn(2, numFrames)) / 2;
        unitNoise = randn(N, numFrames) + 1i * randn(N, numFrames);

        for frameIdx = 1:numFrames
            txIndex = txIndices(frameIdx);
            txBits = codewordBits(:, txIndex);
            txSymbols = codewordSymbols(:, txIndex);
            noise = noiseStd * unitNoise(:, frameIdx);

            for schemeIdx = 1:numSchemes
                selection = select_pattern( ...
                    txSymbols, schemes(schemeIdx).candidate_set, ...
                    groupIndex, c1, oversamplingFactor);
                patternIdx = encode_pattern(selection.pattern, ...
                    size(schemes(schemeIdx).candidate_set, 2));
                hEff = gains(1, frameIdx) * ...
                    patternCache{schemeIdx}.h1{patternIdx} + ...
                    gains(2, frameIdx) * ...
                    patternCache{schemeIdx}.h2{patternIdx};
                received = hEff * txSymbols + noise;
                candidateRx = hEff * codewordSymbols;
                distances = sum(abs(candidateRx - received).^2, 1);
                [~, detectedIndex] = min(distances);
                errorBits(snrIdx, schemeIdx) = ...
                    errorBits(snrIdx, schemeIdx) + ...
                    sum(codewordBits(:, detectedIndex) ~= txBits);
                selectedRisk(snrIdx, schemeIdx) = ...
                    selectedRisk(snrIdx, schemeIdx) + ...
                    patternCache{schemeIdx}.exact_risk(patternIdx);
            end
        end
        fprintf('Selection-aware ML SNR %g dB complete\n', snrDb(snrIdx));
    end

    ber = errorBits ./ totalBits;
    riskProbability = selectedRisk / numFrames;
    diversityEstimate = estimate_slopes(snrDb, ber, errorBits);

    results.config = struct( ...
        'N', N, 'V', V, 'alpha_max', alphaMax, 'c1', c1, ...
        'delay2', delay2, 'doppler2', doppler2, ...
        'snr_db', snrDb, 'num_frames', numFrames, ...
        'oversampling_factor', oversamplingFactor, 'seed', seed);
    results.schemes = string({schemes.name});
    results.ber = ber;
    results.error_bits = errorBits;
    results.total_bits = totalBits;
    results.selected_exact_risk_probability = riskProbability;
    results.diversity_estimate = diversityEstimate;

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    save(fullfile(outputDir, 'small_n_selection_aware_ml.mat'), 'results');
    write_summary(results, snrDb, outputDir);
    plot_results(results, snrDb, outputDir);
end

function schemes = build_schemes(N, c2Base, localDelta)
    schemes = repmat(struct('name', "", 'candidate_set', []), 6, 1);
    schemes(1).name = "baseline";
    schemes(1).candidate_set = c2Base * ones(N, 1);
    schemes(2).name = "GPS-rational";
    schemes(2).candidate_set = afdm.chirp.gps_candidate_set_variant( ...
        N, 'rational');
    schemes(3).name = "GPS-irrational-k1";
    schemes(3).candidate_set = afdm.chirp.gps_candidate_set_variant( ...
        N, 'irrational', struct('precision_k', 1));
    schemes(4).name = "GPS-irrational-k2";
    schemes(4).candidate_set = afdm.chirp.gps_candidate_set_variant( ...
        N, 'irrational', struct('precision_k', 2));
    schemes(5).name = "proposed-2";
    schemes(5).candidate_set = c2Base + ...
        repmat([-localDelta, localDelta], N, 1);
    schemes(6).name = "proposed-3";
    schemes(6).candidate_set = c2Base + ...
        repmat([0, -localDelta, localDelta], N, 1);
end

function cache = build_pattern_cache( ...
        schemes, groupIndex, V, N, c1, delay2, doppler2)
    cache = cell(numel(schemes), 1);
    for schemeIdx = 1:numel(schemes)
        numCandidates = size(schemes(schemeIdx).candidate_set, 2);
        patterns = enumerate_patterns(numCandidates, V);
        h1 = cell(size(patterns, 1), 1);
        h2 = cell(size(patterns, 1), 1);
        exactRisk = false(size(patterns, 1), 1);
        for patternIdx = 1:size(patterns, 1)
            c2 = select_c2( ...
                schemes(schemeIdx).candidate_set, groupIndex, ...
                patterns(patternIdx, :));
            h1{patternIdx} = afdm.analysis.build_path_matrix( ...
                N, c1, c2, 0, 0);
            h2{patternIdx} = afdm.analysis.build_path_matrix( ...
                N, c1, c2, delay2, doppler2);
            cycle = afdm.analysis.monomial_cycle_analysis( ...
                h1{patternIdx}' * h2{patternIdx}, [-2i; 2i], 1e-12);
            exactRisk(patternIdx) = cycle.num_compatible_eigenvectors > 0;
        end
        cache{schemeIdx}.patterns = patterns;
        cache{schemeIdx}.h1 = h1;
        cache{schemeIdx}.h2 = h2;
        cache{schemeIdx}.exact_risk = exactRisk;
    end
end

function selection = select_pattern( ...
        symbols, candidateSet, groupIndex, c1, oversamplingFactor)
    if size(candidateSet, 2) == 1
        selection.pattern = ones(max(groupIndex), 1);
        return;
    end
    selected = afdm.search.greedy_group_papr_selection( ...
        symbols, numel(symbols), c1, candidateSet, groupIndex, ...
        oversamplingFactor);
    selection.pattern = selected.selected_candidate_index;
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

function slopes = estimate_slopes(snrDb, ber, errorBits)
    snrLinear = 10.^(snrDb(:) / 10);
    slopes = NaN(1, size(ber, 2));
    for schemeIdx = 1:size(ber, 2)
        valid = snrDb(:) >= 16 & ber(:, schemeIdx) > 0 & ...
            errorBits(:, schemeIdx) >= 20;
        if nnz(valid) >= 3
            fit = polyfit(log10(snrLinear(valid)), ...
                log10(ber(valid, schemeIdx)), 1);
            slopes(schemeIdx) = -fit(1);
        end
    end
end

function write_summary(results, snrDb, outputDir)
    summary = table(results.schemes(:), results.diversity_estimate(:), ...
        'VariableNames', {'scheme', 'estimated_diversity'});
    for snrIdx = 1:numel(snrDb)
        summary.(sprintf('ber_%ddB', round(snrDb(snrIdx)))) = ...
            results.ber(snrIdx, :).';
        summary.(sprintf('risk_%ddB', round(snrDb(snrIdx)))) = ...
            results.selected_exact_risk_probability(snrIdx, :).';
    end
    writetable(summary, ...
        fullfile(outputDir, 'small_n_selection_aware_ml.csv'));
    disp(summary);
end

function plot_results(results, snrDb, outputDir)
    berPlot = max(results.ber, 0.5 ./ results.total_bits);
    figure('Name', 'Selection-aware exact ML', 'Color', 'w');
    semilogy(snrDb, berPlot, 'LineWidth', 1.8);
    grid on;
    xlabel('SNR (dB)');
    ylabel('BER');
    title('Exact ML with data-driven PAPR pattern selection');
    legend(cellstr(results.schemes), 'Location', 'southwest');
    saveas(gcf, fullfile(outputDir, 'small_n_selection_aware_ml.png'));
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = defaultValue;
    end
end
