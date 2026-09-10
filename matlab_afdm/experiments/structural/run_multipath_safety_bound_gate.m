function results = run_multipath_safety_bound_gate(options)
%RUN_MULTIPATH_SAFETY_BOUND_GATE Test scalable P-path safety certificates.

    rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    setup_paths(rootDir);

    if nargin < 1
        options = struct();
    end

    N = get_option(options, 'N', 64);
    V = get_option(options, 'V', 4);
    alphaMax = get_option(options, 'alpha_max', 3);
    c1 = get_option(options, 'c1', (2 * alphaMax + 1) / (2 * N));
    c2Base = get_option(options, 'c2_base', sqrt(2) / (10 * N));
    localDelta = get_option(options, 'proposed_delta', c2Base / 16);
    maxWeightTwo = get_option(options, 'max_weight_two', 1000);
    numDenseEvents = get_option(options, 'num_dense_events', 500);
    seed = get_option(options, 'seed', 20260910);
    tolerance = get_option(options, 'tolerance', 1e-10);

    supportSets = { ...
        [0, 0; 2, 2; 5, -3], ...
        [0, 0; 2, 2; 8, 0], ...
        [0, 0; 2, 2; 5, -3; 8, 0]};
    schemes = build_schemes(N, V, c2Base, localDelta);
    [eventSet, eventClass] = build_event_set( ...
        N, c1, schemes, supportSets, maxWeightTwo, numDenseEvents, seed);

    numRows = numel(supportSets) * numel(schemes);
    rows = repmat(empty_row(), numRows, 1);
    details = cell(numRows, 1);
    rowIdx = 0;

    for supportIdx = 1:numel(supportSets)
        support = supportSets{supportIdx};
        for schemeIdx = 1:numel(schemes)
            rowIdx = rowIdx + 1;
            pathResponses = build_path_responses( ...
                N, c1, schemes(schemeIdx).c2, support, eventSet);
            metrics = evaluate_certificates(pathResponses, tolerance);

            row = empty_row();
            row.support_name = sprintf('P%d-set%d', size(support, 1), supportIdx);
            row.num_paths = size(support, 1);
            row.scheme = schemes(schemeIdx).name;
            row.num_events = size(eventSet, 2);
            row.min_exact_lambda = min(metrics.exact_lambda);
            row.min_gershgorin_bound = min(metrics.gershgorin_bound);
            row.min_coherence_bound = min(metrics.coherence_bound);
            row.positive_gershgorin_fraction = mean(metrics.gershgorin_bound > 0);
            row.positive_coherence_fraction = mean(metrics.coherence_bound > 0);
            row.exact_rank_loss_events = sum(metrics.exact_lambda < tolerance);
            row.pairwise_rank_loss_events = sum(metrics.pairwise_dependent);
            row.collective_only_rank_loss_events = sum( ...
                metrics.exact_lambda < tolerance & ~metrics.pairwise_dependent);
            row.safe_but_uncertified_fraction = mean( ...
                metrics.exact_lambda >= tolerance & metrics.gershgorin_bound <= 0);
            row.bound_exact_correlation = safe_correlation( ...
                metrics.gershgorin_bound, metrics.exact_lambda);
            row.mean_bound_gap = mean( ...
                metrics.exact_lambda - metrics.gershgorin_bound);
            row.bound_time_ms = 1000 * metrics.bound_time_seconds;
            row.eigen_time_ms = 1000 * metrics.eigen_time_seconds;
            rows(rowIdx) = row;

            details{rowIdx}.support = support;
            details{rowIdx}.scheme = schemes(schemeIdx).name;
            details{rowIdx}.event_class = eventClass;
            details{rowIdx}.metrics = metrics;
        end
    end

    summary = struct2table(rows);
    classSummary = summarize_event_classes(details, eventClass, tolerance);
    results.config = struct( ...
        'N', N, 'V', V, 'c1', c1, 'c2_base', c2Base, ...
        'local_delta', localDelta, 'max_weight_two', maxWeightTwo, ...
        'num_dense_events', numDenseEvents, 'seed', seed, ...
        'tolerance', tolerance);
    results.support_sets = supportSets;
    results.schemes = schemes;
    results.event_set = eventSet;
    results.event_class = eventClass;
    results.summary = summary;
    results.class_summary = classSummary;
    results.details = details;

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    save(fullfile(outputDir, 'multipath_safety_bound_gate.mat'), 'results');
    writetable(summary, fullfile(outputDir, 'multipath_safety_bound_gate.csv'));
    writetable(classSummary, ...
        fullfile(outputDir, 'multipath_safety_bound_by_event_class.csv'));
    plot_results(summary, details, outputDir);

    fprintf('Multipath safety bound gate: %d events\n', size(eventSet, 2));
    disp(summary);
    disp(classSummary);
end

function schemes = build_schemes(N, V, c2Base, localDelta)
    pattern = [2, 2, 1, 1];
    gpsRational = afdm.chirp.gps_candidate_set_variant(N, 'rational');
    gpsSafe = afdm.chirp.gps_candidate_set_variant( ...
        N, 'phase_offset', struct('phase_offset', 0.3));

    schemes = repmat(struct('name', "", 'c2', []), 4, 1);
    schemes(1).name = "baseline";
    schemes(1).c2 = c2Base * ones(N, 1);
    schemes(2).name = "GPS-rational";
    schemes(2).c2 = select_grouped_c2(gpsRational, N, V, pattern);
    schemes(3).name = "local-proposed";
    schemes(3).c2 = afdm.chirp.build_proposed_pattern( ...
        N, V, pattern, c2Base, localDelta);
    schemes(4).name = "pairwise-safe-GPS";
    schemes(4).c2 = select_grouped_c2(gpsSafe, N, V, pattern);
end

function c2 = select_grouped_c2(candidateSet, N, V, pattern)
    groupIndex = afdm.chirp.group_index(N, V, 'contiguous');
    c2 = zeros(N, 1);
    for groupIdx = 1:V
        selected = groupIndex == groupIdx;
        c2(selected) = candidateSet(selected, pattern(groupIdx));
    end
end

function [eventSet, eventClass] = build_event_set( ...
        N, c1, schemes, supportSets, maxWeightTwo, numDenseEvents, seed)
    sparseSet = afdm.delta.generate_sparse_difference_set(N, 2, 'psk', 2);
    weightOne = sparseSet.delta(:, sparseSet.weight == 1);
    weightTwo = sparseSet.delta(:, sparseSet.weight == 2);
    rng(seed, 'twister');
    sampledCount = min(maxWeightTwo, size(weightTwo, 2));
    sampledWeightTwo = weightTwo(:, randperm(size(weightTwo, 2), sampledCount));

    denseEvents = zeros(N, numDenseEvents);
    for eventIdx = 1:numDenseEvents
        firstBits = randi([0, 1], N, 1);
        secondBits = randi([0, 1], N, 1);
        while isequal(firstBits, secondBits)
            secondBits = randi([0, 1], N, 1);
        end
        denseEvents(:, eventIdx) = ...
            1i * (1 - 2 * firstBits) - 1i * (1 - 2 * secondBits);
    end

    eventSet = [weightOne, sampledWeightTwo, denseEvents];
    eventClass = [ ...
        repmat("weight-1", 1, size(weightOne, 2)), ...
        repmat("weight-2", 1, size(sampledWeightTwo, 2)), ...
        repmat("dense", 1, size(denseEvents, 2))];

    differenceAlphabet = afdm.delta.generate_sparse_difference_set( ...
        1, 2, 'psk', 1).difference_alphabet;
    allSupports = unique(vertcat(supportSets{:}), 'rows', 'stable');
    for schemeIdx = 1:numel(schemes)
        pathMatrices = cell(size(allSupports, 1), 1);
        for supportIdx = 1:size(allSupports, 1)
            pathMatrices{supportIdx} = afdm.analysis.build_path_matrix( ...
                N, c1, schemes(schemeIdx).c2, ...
                allSupports(supportIdx, 1), allSupports(supportIdx, 2));
        end
        for firstIdx = 1:numel(pathMatrices)-1
            for secondIdx = firstIdx+1:numel(pathMatrices)
                cycle = afdm.analysis.monomial_cycle_analysis( ...
                    pathMatrices{firstIdx}' * pathMatrices{secondIdx}, ...
                    differenceAlphabet, 1e-12);
                eventSet(:, end + 1) = cycle.closest_delta; %#ok<AGROW>
                eventClass(end + 1) = "pair-cycle"; %#ok<AGROW>
            end
        end
    end

    valid = vecnorm(eventSet, 2, 1) > 0;
    eventSet = eventSet(:, valid);
    eventClass = eventClass(valid).';
    eventSet = eventSet ./ vecnorm(eventSet, 2, 1);
end

function responses = build_path_responses(N, c1, c2, support, eventSet)
    numPaths = size(support, 1);
    responses = cell(numPaths, 1);
    for pathIdx = 1:numPaths
        pathMatrix = afdm.analysis.build_path_matrix( ...
            N, c1, c2, support(pathIdx, 1), support(pathIdx, 2));
        responses{pathIdx} = pathMatrix * eventSet;
    end
end

function metrics = evaluate_certificates(responses, tolerance)
    numPaths = numel(responses);
    numEvents = size(responses{1}, 2);
    gramMatrices = zeros(numPaths, numPaths, numEvents);

    for firstIdx = 1:numPaths
        for secondIdx = firstIdx:numPaths
            values = sum(conj(responses{firstIdx}) .* responses{secondIdx}, 1);
            gramMatrices(firstIdx, secondIdx, :) = values;
            gramMatrices(secondIdx, firstIdx, :) = conj(values);
        end
    end

    tic;
    gershgorinBound = zeros(numEvents, 1);
    coherenceBound = zeros(numEvents, 1);
    pairwiseDependent = false(numEvents, 1);
    for eventIdx = 1:numEvents
        gram = gramMatrices(:, :, eventIdx);
        diagonal = real(diag(gram));
        offDiagonal = abs(gram - diag(diag(gram)));
        gershgorinBound(eventIdx) = min(diagonal - sum(offDiagonal, 2));
        normalized = diag(1 ./ sqrt(max(diagonal, eps))) * gram * ...
            diag(1 ./ sqrt(max(diagonal, eps)));
        offNormalized = abs(normalized - eye(numPaths));
        coherence = max(offNormalized(:));
        coherenceBound(eventIdx) = 1 - (numPaths - 1) * coherence;
        pairwiseDependent(eventIdx) = coherence >= 1 - tolerance;
    end
    boundTime = toc;

    tic;
    exactLambda = zeros(numEvents, 1);
    for eventIdx = 1:numEvents
        eigenvalues = eig(gramMatrices(:, :, eventIdx), 'vector');
        exactLambda(eventIdx) = max(min(real(eigenvalues)), 0);
    end
    eigenTime = toc;

    metrics.exact_lambda = exactLambda;
    metrics.gershgorin_bound = gershgorinBound;
    metrics.coherence_bound = coherenceBound;
    metrics.pairwise_dependent = pairwiseDependent;
    metrics.bound_time_seconds = boundTime;
    metrics.eigen_time_seconds = eigenTime;
end

function summary = summarize_event_classes(details, eventClass, tolerance)
    classes = unique(eventClass, 'stable');
    rows = repmat(struct( ...
        'event_class', "", 'num_events', NaN, ...
        'positive_bound_fraction', NaN, 'exact_rank_loss_events', NaN, ...
        'collective_only_events', NaN), numel(classes), 1);

    for classIdx = 1:numel(classes)
        selected = eventClass == classes(classIdx);
        positive = 0;
        exactLoss = 0;
        collectiveOnly = 0;
        total = 0;
        for detailIdx = 1:numel(details)
            metrics = details{detailIdx}.metrics;
            positive = positive + sum(metrics.gershgorin_bound(selected) > 0);
            exactLoss = exactLoss + sum(metrics.exact_lambda(selected) < tolerance);
            collectiveOnly = collectiveOnly + sum( ...
                metrics.exact_lambda(selected) < tolerance & ...
                ~metrics.pairwise_dependent(selected));
            total = total + sum(selected);
        end
        rows(classIdx).event_class = classes(classIdx);
        rows(classIdx).num_events = sum(selected);
        rows(classIdx).positive_bound_fraction = positive / total;
        rows(classIdx).exact_rank_loss_events = exactLoss;
        rows(classIdx).collective_only_events = collectiveOnly;
    end
    summary = struct2table(rows);
end

function value = safe_correlation(first, second)
    if std(first) <= eps || std(second) <= eps
        value = NaN;
        return;
    end
    matrix = corrcoef(first, second);
    value = matrix(1, 2);
end

function plot_results(summary, details, outputDir)
    fig = figure('Name', 'Multipath safety bound gate', 'Color', 'w');
    tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    supportNames = unique(string(summary.support_name), 'stable');
    schemeNames = unique(summary.scheme, 'stable');
    exactMatrix = grouped_metric(summary, supportNames, schemeNames, ...
        'min_exact_lambda');
    rateMatrix = grouped_metric(summary, supportNames, schemeNames, ...
        'positive_gershgorin_fraction');
    costRatioMatrix = grouped_metric(summary, supportNames, schemeNames, ...
        'bound_time_ms') ./ grouped_metric(summary, supportNames, schemeNames, ...
        'eigen_time_ms');

    nexttile;
    bar(max(exactMatrix, eps), 'grouped');
    set(gca, 'YScale', 'log');
    set(gca, 'XTickLabel', supportNames);
    grid on;
    ylabel('minimum exact \lambda_{min}');
    title('Multipath conditioning');
    legend(schemeNames, 'Location', 'best');

    nexttile;
    bar(rateMatrix, 'grouped');
    set(gca, 'XTickLabel', supportNames);
    grid on;
    ylim([0, 1]);
    ylabel('positive-bound fraction');
    title('Gershgorin certification rate');

    nexttile;
    detail = details{end}.metrics;
    scatter(detail.gershgorin_bound, detail.exact_lambda, 12, 'filled');
    hold on;
    limits = [min(detail.gershgorin_bound), max(detail.exact_lambda)];
    plot(limits, limits, 'k--');
    grid on;
    xlabel('Gershgorin lower bound');
    ylabel('exact \lambda_{min}');
    title('Bound tightness, P=4 safe GPS');

    nexttile;
    bar(costRatioMatrix, 'grouped');
    set(gca, 'XTickLabel', supportNames);
    grid on;
    ylabel('bound time / eigen time');
    yline(1, 'k--');
    title('Measured certificate cost ratio');

    exportgraphics(fig, ...
        fullfile(outputDir, 'multipath_safety_bound_gate.png'), ...
        'Resolution', 200);
    close(fig);
end

function matrix = grouped_metric(summary, supportNames, schemeNames, variableName)
    matrix = zeros(numel(supportNames), numel(schemeNames));
    for supportIdx = 1:numel(supportNames)
        for schemeIdx = 1:numel(schemeNames)
            selected = string(summary.support_name) == supportNames(supportIdx) & ...
                summary.scheme == schemeNames(schemeIdx);
            matrix(supportIdx, schemeIdx) = summary.(variableName)(selected);
        end
    end
end

function row = empty_row()
    row = struct( ...
        'support_name', "", ...
        'num_paths', NaN, ...
        'scheme', "", ...
        'num_events', NaN, ...
        'min_exact_lambda', NaN, ...
        'min_gershgorin_bound', NaN, ...
        'min_coherence_bound', NaN, ...
        'positive_gershgorin_fraction', NaN, ...
        'positive_coherence_fraction', NaN, ...
        'exact_rank_loss_events', NaN, ...
        'pairwise_rank_loss_events', NaN, ...
        'collective_only_rank_loss_events', NaN, ...
        'safe_but_uncertified_fraction', NaN, ...
        'bound_exact_correlation', NaN, ...
        'mean_bound_gap', NaN, ...
        'bound_time_ms', NaN, ...
        'eigen_time_ms', NaN);
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = defaultValue;
    end
end
