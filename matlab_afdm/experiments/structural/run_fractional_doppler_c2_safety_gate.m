function results = run_fractional_doppler_c2_safety_gate(options)
%RUN_FRACTIONAL_DOPPLER_C2_SAFETY_GATE Track c2-codeword risk off grid.
%
% This is a direction gate, not an exhaustive diversity proof. The candidate
% error set contains all BPSK events of weight one, a deterministic sample of
% weight-two events, and the integer-Doppler cycle events of every scheme.

    rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    setup_paths(rootDir);

    if nargin < 1
        options = struct();
    end

    N = get_option(options, 'N', 64);
    V = get_option(options, 'V', 4);
    alphaMax = get_option(options, 'alpha_max', 3);
    betaValues = get_option(options, 'beta_values', -0.5:0.05:0.5);
    maxSampledWeightTwo = get_option(options, 'max_sampled_weight_two', 1500);
    seed = get_option(options, 'seed', 20260909);
    snrDb = get_option(options, 'pep_snr_db', 30);
    c1 = get_option(options, 'c1', (2 * alphaMax + 1) / (2 * N));
    c2Base = get_option(options, 'c2_base', sqrt(2) / (10 * N));
    localDelta = get_option(options, 'proposed_delta', c2Base / 16);
    delays = get_option(options, 'delays', [0, 2]);
    integerDopplers = get_option(options, 'integer_dopplers', [0, 2]);

    if numel(delays) ~= 2 || numel(integerDopplers) ~= 2
        error('The first gate currently requires exactly two paths.');
    end
    if mod(N, V) ~= 0
        error('N must be divisible by V.');
    end

    schemes = build_schemes(N, V, c2Base, localDelta);
    differenceAlphabet = afdm.delta.generate_sparse_difference_set( ...
        1, 2, 'psk', 1).difference_alphabet;
    eventSet = build_candidate_event_set( ...
        N, schemes, c1, delays, integerDopplers, differenceAlphabet, ...
        maxSampledWeightTwo, seed);
    integerWorstIndices = find_integer_worst_events( ...
        schemes, eventSet, N, c1, delays, integerDopplers);

    numSchemes = numel(schemes);
    numBeta = numel(betaValues);
    numRows = numSchemes * numBeta;
    rows = repmat(empty_row(), numRows, 1);
    leakageRows = repmat(empty_leakage_row(), numBeta, 1);
    rowIdx = 0;

    for betaIdx = 1:numBeta
        dopplers = integerDopplers;
        dopplers(2) = dopplers(2) + betaValues(betaIdx);
        magnitudeReferences = cell(numSchemes, 1);

        for schemeIdx = 1:numSchemes
            pathMatrices = build_path_matrices( ...
                N, c1, schemes(schemeIdx).c2, delays, dopplers);
            metrics = evaluate_events(pathMatrices, eventSet, snrDb);
            trackedIdx = integerWorstIndices(schemeIdx);

            rowIdx = rowIdx + 1;
            row = empty_row();
            row.scheme = schemes(schemeIdx).name;
            row.fractional_doppler = betaValues(betaIdx);
            row.min_sigma = metrics.sigma(metrics.worst_index);
            row.tracked_integer_event_sigma = metrics.sigma(trackedIdx);
            row.min_gram_determinant = metrics.gram_determinant(metrics.worst_index);
            row.worst_pep_surrogate = metrics.pep_surrogate(metrics.worst_index);
            row.near_rank_loss_events = sum(metrics.sigma < 1e-8);
            row.worst_event_index = metrics.worst_index;
            row.integer_worst_event_index = trackedIdx;
            rows(rowIdx) = row;
            magnitudeReferences{schemeIdx} = abs(pathMatrices{2});
        end

        leakage = leakage_metrics(magnitudeReferences{1});
        maxMagnitudeDifference = 0;
        for schemeIdx = 2:numSchemes
            maxMagnitudeDifference = max(maxMagnitudeDifference, ...
                max(abs(magnitudeReferences{schemeIdx}(:) - ...
                magnitudeReferences{1}(:))));
        end
        leakageRow = empty_leakage_row();
        leakageRow.fractional_doppler = betaValues(betaIdx);
        leakageRow.mean_top1_energy = leakage.top1;
        leakageRow.mean_top3_energy = leakage.top3;
        leakageRow.mean_entries_for_95pct = leakage.entries95;
        leakageRow.max_c2_magnitude_difference = maxMagnitudeDifference;
        leakageRows(betaIdx) = leakageRow;
    end

    auditTable = struct2table(rows);
    leakageTable = struct2table(leakageRows);
    results.config = struct( ...
        'N', N, 'V', V, 'alpha_max', alphaMax, 'c1', c1, ...
        'c2_base', c2Base, 'local_delta', localDelta, ...
        'delays', delays, 'integer_dopplers', integerDopplers, ...
        'beta_values', betaValues, 'event_count', size(eventSet, 2), ...
        'pep_snr_db', snrDb, 'seed', seed);
    results.schemes = schemes;
    results.event_set = eventSet;
    results.integer_worst_indices = integerWorstIndices;
    results.audit_table = auditTable;
    results.leakage_table = leakageTable;

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    outputTag = sprintf('fractional_doppler_c2_safety_gate_d%d_a%d', ...
        delays(2), integerDopplers(2));
    save(fullfile(outputDir, [outputTag '.mat']), 'results');
    writetable(auditTable, ...
        fullfile(outputDir, [outputTag '.csv']));
    writetable(leakageTable, ...
        fullfile(outputDir, [outputTag '_leakage.csv']));
    plot_results(auditTable, leakageTable, schemes, outputDir, outputTag);

    fprintf('Fractional-Doppler c2 safety gate: %d candidate error events\n', ...
        size(eventSet, 2));
    print_key_rows(auditTable);
    disp(leakageTable);
end

function schemes = build_schemes(N, V, c2Base, localDelta)
    pattern = [2, 2, 1, 1];
    if V ~= numel(pattern)
        error('The default gate currently uses V=4 and pattern [2 2 1 1].');
    end

    gpsRational = afdm.chirp.gps_candidate_set_variant(N, 'rational');
    gpsSafeCandidates = afdm.chirp.gps_candidate_set_variant( ...
        N, 'phase_offset', struct('phase_offset', 0.3));

    schemes = repmat(struct('name', "", 'c2', []), 4, 1);
    schemes(1).name = "baseline";
    schemes(1).c2 = c2Base * ones(N, 1);
    schemes(2).name = "GPS-rational";
    schemes(2).c2 = select_grouped_c2(gpsRational, N, V, pattern);
    schemes(3).name = "local-proposed";
    schemes(3).c2 = afdm.chirp.build_proposed_pattern( ...
        N, V, pattern, c2Base, localDelta);
    schemes(4).name = "integer-safe-GPS";
    schemes(4).c2 = select_grouped_c2(gpsSafeCandidates, N, V, pattern);
end

function c2 = select_grouped_c2(candidateSet, N, V, pattern)
    groupIndex = afdm.chirp.group_index(N, V, 'contiguous');
    c2 = zeros(N, 1);
    for groupIdx = 1:V
        selected = groupIndex == groupIdx;
        c2(selected) = candidateSet(selected, pattern(groupIdx));
    end
end

function eventSet = build_candidate_event_set( ...
        N, schemes, c1, delays, dopplers, differenceAlphabet, ...
        maxSampledWeightTwo, seed)
    sparseSet = afdm.delta.generate_sparse_difference_set(N, 2, 'psk', 2);
    weightOne = sparseSet.delta(:, sparseSet.weight == 1);
    weightTwo = sparseSet.delta(:, sparseSet.weight == 2);

    rng(seed, 'twister');
    sampleCount = min(maxSampledWeightTwo, size(weightTwo, 2));
    sampledIndices = randperm(size(weightTwo, 2), sampleCount);
    eventSet = [weightOne, weightTwo(:, sampledIndices)];

    for schemeIdx = 1:numel(schemes)
        pathMatrices = build_path_matrices( ...
            N, c1, schemes(schemeIdx).c2, delays, dopplers);
        cycle = afdm.analysis.monomial_cycle_analysis( ...
            pathMatrices{1}' * pathMatrices{2}, differenceAlphabet, 1e-12);
        eventSet(:, end + 1) = cycle.closest_delta; %#ok<AGROW>
    end

    norms = vecnorm(eventSet, 2, 1);
    eventSet = eventSet(:, norms > 0);
    eventSet = eventSet ./ vecnorm(eventSet, 2, 1);
end

function indices = find_integer_worst_events( ...
        schemes, eventSet, N, c1, delays, dopplers)
    indices = zeros(numel(schemes), 1);
    for schemeIdx = 1:numel(schemes)
        pathMatrices = build_path_matrices( ...
            N, c1, schemes(schemeIdx).c2, delays, dopplers);
        metrics = evaluate_events(pathMatrices, eventSet, 30);
        indices(schemeIdx) = metrics.worst_index;
    end
end

function pathMatrices = build_path_matrices(N, c1, c2, delays, dopplers)
    pathMatrices = cell(1, 2);
    for pathIdx = 1:2
        pathMatrices{pathIdx} = afdm.analysis.build_path_matrix( ...
            N, c1, c2, delays(pathIdx), dopplers(pathIdx));
    end
end

function metrics = evaluate_events(pathMatrices, eventSet, snrDb)
    firstResponse = pathMatrices{1} * eventSet;
    secondResponse = pathMatrices{2} * eventSet;
    firstEnergy = sum(abs(firstResponse).^2, 1);
    secondEnergy = sum(abs(secondResponse).^2, 1);
    crossTerm = sum(conj(firstResponse) .* secondResponse, 1);

    discriminant = sqrt(max( ...
        (firstEnergy - secondEnergy).^2 + 4 * abs(crossTerm).^2, 0));
    lambdaMin = max((firstEnergy + secondEnergy - discriminant) / 2, 0);
    lambdaMax = max((firstEnergy + secondEnergy + discriminant) / 2, 0);
    sigma = sqrt(lambdaMin);
    gramDeterminant = max(firstEnergy .* secondEnergy - abs(crossTerm).^2, 0);
    snrLinear = 10^(snrDb / 10);
    pepSurrogate = 1 ./ ((1 + snrLinear * lambdaMin / 4) .* ...
        (1 + snrLinear * lambdaMax / 4));

    [~, worstIndex] = min(sigma);
    metrics.sigma = sigma;
    metrics.gram_determinant = gramDeterminant;
    metrics.pep_surrogate = pepSurrogate;
    metrics.worst_index = worstIndex;
end

function metrics = leakage_metrics(magnitudeMatrix)
    energy = magnitudeMatrix.^2;
    energy = energy ./ max(sum(energy, 2), eps);
    sortedEnergy = sort(energy, 2, 'descend');
    cumulativeEnergy = cumsum(sortedEnergy, 2);
    metrics.top1 = mean(cumulativeEnergy(:, 1));
    metrics.top3 = mean(cumulativeEnergy(:, min(3, size(energy, 2))));
    entries95 = zeros(size(energy, 1), 1);
    for rowIdx = 1:size(energy, 1)
        firstIndex = find(cumulativeEnergy(rowIdx, :) >= 0.95, 1, 'first');
        entries95(rowIdx) = firstIndex;
    end
    metrics.entries95 = mean(entries95);
end

function plot_results(auditTable, leakageTable, schemes, outputDir, outputTag)
    fig = figure('Name', 'Fractional Doppler c2 safety gate', 'Color', 'w');
    tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    colors = lines(numel(schemes));

    nexttile;
    hold on;
    for schemeIdx = 1:numel(schemes)
        rows = auditTable.scheme == schemes(schemeIdx).name;
        semilogy(auditTable.fractional_doppler(rows), ...
            max(auditTable.min_sigma(rows), eps), 'o-', ...
            'Color', colors(schemeIdx, :), 'LineWidth', 1.4);
    end
    set(gca, 'YScale', 'log');
    grid on;
    xlabel('fractional Doppler \beta');
    ylabel('minimum \sigma_{min}');
    title('Worst event in candidate set');

    nexttile;
    hold on;
    for schemeIdx = 1:numel(schemes)
        rows = auditTable.scheme == schemes(schemeIdx).name;
        semilogy(auditTable.fractional_doppler(rows), ...
            max(auditTable.tracked_integer_event_sigma(rows), eps), 'o-', ...
            'Color', colors(schemeIdx, :), 'LineWidth', 1.4);
    end
    set(gca, 'YScale', 'log');
    grid on;
    xlabel('fractional Doppler \beta');
    ylabel('tracked \sigma_{min}');
    title('Integer worst event');

    nexttile;
    hold on;
    for schemeIdx = 1:numel(schemes)
        rows = auditTable.scheme == schemes(schemeIdx).name;
        semilogy(auditTable.fractional_doppler(rows), ...
            max(auditTable.worst_pep_surrogate(rows), eps), 'o-', ...
            'Color', colors(schemeIdx, :), 'LineWidth', 1.4);
    end
    set(gca, 'YScale', 'log');
    grid on;
    xlabel('fractional Doppler \beta');
    ylabel('PEP surrogate');
    title('30 dB pairwise bound indicator');

    nexttile;
    plot(leakageTable.fractional_doppler, leakageTable.mean_top1_energy, ...
        'o-', 'LineWidth', 1.4);
    hold on;
    plot(leakageTable.fractional_doppler, leakageTable.mean_top3_energy, ...
        's-', 'LineWidth', 1.4);
    grid on;
    ylim([0, 1.02]);
    xlabel('fractional Doppler \beta');
    ylabel('captured path energy');
    legend({'top 1', 'top 3'}, 'Location', 'best');
    title('Single-path leakage');

    legend(nexttile(1), string({schemes.name}), 'Location', 'best');
    exportgraphics(fig, ...
        fullfile(outputDir, [outputTag '.png']), ...
        'Resolution', 200);
    close(fig);
end

function print_key_rows(auditTable)
    targetBeta = [0, 0.25, 0.5];
    selected = false(height(auditTable), 1);
    for idx = 1:numel(targetBeta)
        selected = selected | ...
            abs(auditTable.fractional_doppler - targetBeta(idx)) < 1e-12;
    end
    disp(auditTable(selected, {'scheme', 'fractional_doppler', ...
        'min_sigma', 'tracked_integer_event_sigma', ...
        'worst_pep_surrogate', 'near_rank_loss_events'}));
end

function row = empty_row()
    row = struct( ...
        'scheme', "", ...
        'fractional_doppler', NaN, ...
        'min_sigma', NaN, ...
        'tracked_integer_event_sigma', NaN, ...
        'min_gram_determinant', NaN, ...
        'worst_pep_surrogate', NaN, ...
        'near_rank_loss_events', NaN, ...
        'worst_event_index', NaN, ...
        'integer_worst_event_index', NaN);
end

function row = empty_leakage_row()
    row = struct( ...
        'fractional_doppler', NaN, ...
        'mean_top1_energy', NaN, ...
        'mean_top3_energy', NaN, ...
        'mean_entries_for_95pct', NaN, ...
        'max_c2_magnitude_difference', NaN);
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = defaultValue;
    end
end
