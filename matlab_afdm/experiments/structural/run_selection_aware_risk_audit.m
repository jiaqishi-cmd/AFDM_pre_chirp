function results = run_selection_aware_risk_audit(options)
%RUN_SELECTION_AWARE_RISK_AUDIT Connect PAPR selection to receiver risk.

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
    delay2 = get_option(options, 'delay2', 2);
    doppler2 = get_option(options, 'doppler2', 2);
    selectionFrames = get_option(options, 'selection_frames', 5000);
    berFrames = get_option(options, 'ber_frames', 20000);
    snrDb = get_option(options, 'snr_db', 28);
    oversamplingFactor = get_option(options, 'oversampling_factor', 4);
    vulnerabilityThreshold = get_option(options, ...
        'vulnerability_threshold', 0.05);
    seed = get_option(options, 'seed', 20260908);

    schemes = build_schemes(N, c2Base, localDelta);
    groupIndex = afdm.chirp.group_index(N, V, 'contiguous');
    modulations = struct( ...
        'name', {"BPSK", "QPSK", "16QAM"}, ...
        'M_mod', {2, 4, 16}, ...
        'mod_type', {'psk', 'psk', 'qam'});

    selectionRows = repmat(empty_selection_row(), ...
        numel(schemes) * numel(modulations), 1);
    patternDetails = cell(numel(schemes), numel(modulations));
    rowCursor = 0;

    for modulationIdx = 1:numel(modulations)
        modulation = modulations(modulationIdx);
        deltaSet = afdm.delta.generate_sparse_difference_set( ...
            1, modulation.M_mod, modulation.mod_type, 1);
        structural = precompute_pattern_structure( ...
            schemes, groupIndex, V, N, c1, delay2, doppler2, ...
            deltaSet.difference_alphabet, vulnerabilityThreshold);
        counts = cell(numel(schemes), 1);
        paprSamples = zeros(selectionFrames, numel(schemes));
        for schemeIdx = 1:numel(schemes)
            counts{schemeIdx} = zeros(size(structural{schemeIdx}.patterns, 1), 1);
        end

        for frameIdx = 1:selectionFrames
            rng(seed + 1000000 * modulationIdx + frameIdx, 'twister');
            [symbols, ~] = afdm.tx.random_data( ...
                N, modulation.M_mod, modulation.mod_type);
            for schemeIdx = 1:numel(schemes)
                [pattern, papr] = select_pattern_for_symbols( ...
                    symbols, schemes(schemeIdx).candidate_set, groupIndex, ...
                    c1, oversamplingFactor);
                patternIdx = encode_pattern(pattern, ...
                    size(schemes(schemeIdx).candidate_set, 2));
                counts{schemeIdx}(patternIdx) = counts{schemeIdx}(patternIdx) + 1;
                paprSamples(frameIdx, schemeIdx) = papr;
            end
        end

        for schemeIdx = 1:numel(schemes)
            rowCursor = rowCursor + 1;
            frequencies = counts{schemeIdx} / selectionFrames;
            detail = structural{schemeIdx};
            row = empty_selection_row();
            row.modulation = modulation.name;
            row.scheme = schemes(schemeIdx).name;
            row.pattern_count = numel(frequencies);
            row.exact_risky_patterns = sum(detail.exact_risk);
            row.vulnerable_patterns = sum(detail.vulnerable);
            row.selected_exact_risk_probability = ...
                sum(frequencies(detail.exact_risk));
            row.selected_vulnerable_probability = ...
                sum(frequencies(detail.vulnerable));
            row.selected_mean_sigma = sum(frequencies .* detail.min_sigma);
            row.mean_papr = mean(paprSamples(:, schemeIdx));
            row.p99_papr = percentile_by_sort( ...
                paprSamples(:, schemeIdx), 0.99);
            selectionRows(rowCursor) = row;
            detail.selection_frequency = frequencies;
            patternDetails{schemeIdx, modulationIdx} = detail;
        end
    end

    selectionSummary = struct2table(selectionRows);
    berSummary = run_selection_aware_ber( ...
        schemes, groupIndex, V, N, c1, delay2, doppler2, ...
        berFrames, snrDb, oversamplingFactor, ...
        vulnerabilityThreshold, seed + 9000000);

    results.config = struct( ...
        'N', N, 'V', V, 'alpha_max', alphaMax, 'c1', c1, ...
        'delay2', delay2, 'doppler2', doppler2, ...
        'selection_frames', selectionFrames, 'ber_frames', berFrames, ...
        'snr_db', snrDb, 'oversampling_factor', oversamplingFactor, ...
        'vulnerability_threshold', vulnerabilityThreshold, 'seed', seed);
    results.selection_summary = selectionSummary;
    results.ber_summary = berSummary;
    results.pattern_details = patternDetails;

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    save(fullfile(outputDir, 'selection_aware_risk_audit.mat'), 'results');
    writetable(selectionSummary, ...
        fullfile(outputDir, 'selection_aware_pattern_risk.csv'));
    writetable(berSummary, ...
        fullfile(outputDir, 'selection_aware_mmse_ber.csv'));
    disp(selectionSummary);
    disp(berSummary);
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

function structural = precompute_pattern_structure( ...
        schemes, groupIndex, V, N, c1, delay2, doppler2, ...
        differenceAlphabet, vulnerabilityThreshold)
    structural = cell(numel(schemes), 1);
    for schemeIdx = 1:numel(schemes)
        numCandidates = size(schemes(schemeIdx).candidate_set, 2);
        patterns = enumerate_patterns(numCandidates, V);
        exactRisk = false(size(patterns, 1), 1);
        minSigma = zeros(size(patterns, 1), 1);
        hEff = cell(size(patterns, 1), 1);
        for patternIdx = 1:size(patterns, 1)
            c2 = select_c2( ...
                schemes(schemeIdx).candidate_set, groupIndex, ...
                patterns(patternIdx, :));
            h1 = afdm.analysis.build_path_matrix(N, c1, c2, 0, 0);
            h2 = afdm.analysis.build_path_matrix( ...
                N, c1, c2, delay2, doppler2);
            cycle = afdm.analysis.monomial_cycle_analysis( ...
                h1' * h2, differenceAlphabet, 1e-12);
            exactRisk(patternIdx) = cycle.num_compatible_eigenvectors > 0;
            delta = cycle.closest_delta;
            psi = afdm.analysis.build_pairwise_path_response( ...
                {h1, h2}, delta / norm(delta));
            metrics = afdm.analysis.pairwise_diversity_metrics(psi, 1e-12);
            minSigma(patternIdx) = metrics.min_singular_value;
            hEff{patternIdx} = (h1 - h2) / sqrt(2);
        end
        detail.patterns = patterns;
        detail.exact_risk = exactRisk;
        detail.min_sigma = minSigma;
        detail.vulnerable = minSigma < vulnerabilityThreshold;
        detail.fixed_channel_h = hEff;
        structural{schemeIdx} = detail;
    end
end

function summary = run_selection_aware_ber( ...
        schemes, groupIndex, V, N, c1, delay2, doppler2, ...
        numFrames, snrDb, oversamplingFactor, vulnerabilityThreshold, seed)
    structure = precompute_pattern_structure( ...
        schemes, groupIndex, V, N, c1, delay2, doppler2, ...
        [-2i; 2i], vulnerabilityThreshold);
    noiseVar = 10^(-snrDb / 10);
    equalizers = cell(numel(schemes), 1);
    for schemeIdx = 1:numel(schemes)
        equalizers{schemeIdx} = cell(numel(structure{schemeIdx}.fixed_channel_h), 1);
        for patternIdx = 1:numel(structure{schemeIdx}.fixed_channel_h)
            hEff = structure{schemeIdx}.fixed_channel_h{patternIdx};
            equalizers{schemeIdx}{patternIdx} = ...
                (hEff' * hEff + noiseVar * eye(N)) \ hEff';
        end
    end

    errors = zeros(numel(schemes), 1);
    selectedExact = zeros(numel(schemes), 1);
    selectedVulnerable = zeros(numel(schemes), 1);
    rng(seed, 'twister');
    for frameIdx = 1:numFrames
        bits = randi([0, 1], N, 1);
        symbols = 1i * (1 - 2 * bits);
        noise = sqrt(noiseVar / 2) * ...
            (randn(N, 1) + 1i * randn(N, 1));
        for schemeIdx = 1:numel(schemes)
            pattern = select_pattern_for_symbols( ...
                symbols, schemes(schemeIdx).candidate_set, groupIndex, ...
                c1, oversamplingFactor);
            patternIdx = encode_pattern(pattern, ...
                size(schemes(schemeIdx).candidate_set, 2));
            hEff = structure{schemeIdx}.fixed_channel_h{patternIdx};
            estimate = equalizers{schemeIdx}{patternIdx} * ...
                (hEff * symbols + noise);
            decisions = pskdemod(estimate, 2, pi / 2, 'OutputType', 'bit');
            errors(schemeIdx) = errors(schemeIdx) + sum(decisions ~= bits);
            selectedExact(schemeIdx) = selectedExact(schemeIdx) + ...
                structure{schemeIdx}.exact_risk(patternIdx);
            selectedVulnerable(schemeIdx) = selectedVulnerable(schemeIdx) + ...
                structure{schemeIdx}.vulnerable(patternIdx);
        end
    end

    scheme = string({schemes.name}).';
    ber = errors / (numFrames * N);
    bit_errors = errors;
    total_bits = numFrames * N * ones(numel(schemes), 1);
    selected_exact_risk_probability = selectedExact / numFrames;
    selected_vulnerable_probability = selectedVulnerable / numFrames;
    summary = table(scheme, ber, bit_errors, total_bits, ...
        selected_exact_risk_probability, selected_vulnerable_probability);
end

function [pattern, papr] = select_pattern_for_symbols( ...
        symbols, candidateSet, groupIndex, c1, oversamplingFactor)
    if size(candidateSet, 2) == 1
        pattern = ones(max(groupIndex), 1);
        waveform = afdm.search.full_waveform( ...
            symbols, candidateSet(:, 1), oversamplingFactor);
        papr = afdm.tx.compute_papr(waveform);
        return;
    end
    selection = afdm.search.greedy_group_papr_selection( ...
        symbols, numel(symbols), c1, candidateSet, groupIndex, ...
        oversamplingFactor);
    pattern = selection.selected_candidate_index;
    papr = selection.papr;
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
    samples = sort(samples(:));
    index = max(1, min(numel(samples), ceil(probability * numel(samples))));
    value = samples(index);
end

function row = empty_selection_row()
    row = struct( ...
        'modulation', "", ...
        'scheme', "", ...
        'pattern_count', NaN, ...
        'exact_risky_patterns', NaN, ...
        'vulnerable_patterns', NaN, ...
        'selected_exact_risk_probability', NaN, ...
        'selected_vulnerable_probability', NaN, ...
        'selected_mean_sigma', NaN, ...
        'mean_papr', NaN, ...
        'p99_papr', NaN);
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = defaultValue;
    end
end
