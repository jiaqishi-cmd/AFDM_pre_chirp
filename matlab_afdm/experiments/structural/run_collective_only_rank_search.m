function results = run_collective_only_rank_search(options)
%RUN_COLLECTIVE_ONLY_RANK_SEARCH Exhaust small-system multipath failures.
%
% A collective-only event has a rank-deficient P=3 path-response matrix
% while every two-column submatrix remains full rank.

    rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    setup_paths(rootDir);

    if nargin < 1
        options = struct();
    end

    N = get_option(options, 'N', 8);
    V = get_option(options, 'V', 4);
    alphaMax = get_option(options, 'alpha_max', 1);
    maxDelay = get_option(options, 'max_delay', 2);
    c1 = get_option(options, 'c1', (2 * alphaMax + 1) / (2 * N));
    c2Base = get_option(options, 'c2_base', sqrt(2) / (10 * N));
    localDelta = get_option(options, 'proposed_delta', c2Base / 16);
    phaseOffset = get_option(options, 'phase_offset', 0.3);
    tolerance = get_option(options, 'tolerance', 1e-10);

    if mod(N, V) ~= 0
        error('N must be divisible by V.');
    end
    if N > 10
        error('This exhaustive search is intended for N <= 10.');
    end

    deltaSet = afdm.delta.generate_sparse_difference_set(N, 2, 'psk', N);
    eventSet = deltaSet.delta ./ sqrt(deltaSet.norm2);
    support = build_support(maxDelay, alphaMax);
    schemeFamilies = build_scheme_families( ...
        N, V, c2Base, localDelta, phaseOffset);
    triples = nchoosek(1:size(support, 1), 3);

    totalPatterns = sum(arrayfun(@(x) size(x.patterns, 1), schemeFamilies));
    rows = repmat(empty_row(), totalPatterns, 1);
    findings = repmat(empty_finding(), 0, 1);
    rowIdx = 0;

    fprintf('========== collective-only rank search ==========\n');
    fprintf('N=%d, events=%d, supports=%d, triples=%d, patterns=%d\n', ...
        N, size(eventSet, 2), size(support, 1), size(triples, 1), totalPatterns);

    for familyIdx = 1:numel(schemeFamilies)
        family = schemeFamilies(familyIdx);
        fprintf('Scanning %s: %d patterns\n', ...
            family.name, size(family.patterns, 1));

        for patternIdx = 1:size(family.patterns, 1)
            rowIdx = rowIdx + 1;
            pattern = family.patterns(patternIdx, :);
            c2 = family.build_c2(pattern);
            correlations = build_pair_correlations( ...
                N, c1, c2, support, eventSet);
            pairAudit = audit_pairs(correlations, tolerance);

            row = empty_row();
            row.scheme = family.name;
            row.pattern = join(string(pattern), '-');
            row.num_events = size(eventSet, 2);
            row.num_support_pairs = nchoosek(size(support, 1), 2);
            row.pairwise_rank_loss_events = pairAudit.rank_loss_count;
            row.max_pair_correlation = pairAudit.max_correlation;
            row.pairwise_safe = pairAudit.rank_loss_count == 0;

            if row.pairwise_safe
                collective = audit_triples( ...
                    correlations, triples, eventSet, support, ...
                    family.name, row.pattern, tolerance);
                row.collective_only_events = collective.count;
                row.min_collective_determinant = collective.min_determinant;
                row.min_collective_lambda = collective.min_lambda;
                if ~isempty(collective.findings)
                    findings = [findings; collective.findings]; %#ok<AGROW>
                end
            end
            rows(rowIdx) = row;
        end
    end

    summary = struct2table(rows);
    if isempty(findings)
        findingTable = struct2table(repmat(empty_finding(), 0, 1));
    else
        findingTable = struct2table(findings);
    end
    familySummary = summarize_families(summary);

    results.config = struct( ...
        'N', N, 'V', V, 'alpha_max', alphaMax, 'max_delay', maxDelay, ...
        'c1', c1, 'c2_base', c2Base, 'local_delta', localDelta, ...
        'phase_offset', phaseOffset, 'tolerance', tolerance);
    results.support = support;
    results.triples = triples;
    results.delta_set = deltaSet;
    results.summary = summary;
    results.family_summary = familySummary;
    results.findings = findingTable;

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    outputTag = sprintf('collective_only_rank_search_N%d_V%d_a%d_l%d', ...
        N, V, alphaMax, maxDelay);
    save(fullfile(outputDir, [outputTag '.mat']), 'results');
    writetable(summary, fullfile(outputDir, [outputTag '.csv']));
    writetable(familySummary, ...
        fullfile(outputDir, [outputTag '_families.csv']));
    if ~isempty(findingTable)
        writetable(findingTable, ...
            fullfile(outputDir, [outputTag '_findings.csv']));
    end
    plot_results(summary, familySummary, outputDir, outputTag);

    disp(familySummary);
    if isempty(findingTable)
        fprintf('No collective-only finite-alphabet rank-loss event found.\n');
    else
        fprintf('Found %d collective-only events.\n', height(findingTable));
        disp(findingTable(1:min(20, height(findingTable)), :));
    end
end

function support = build_support(maxDelay, alphaMax)
    delays = 0:maxDelay;
    dopplers = -alphaMax:alphaMax;
    [delayGrid, dopplerGrid] = ndgrid(delays, dopplers);
    support = [delayGrid(:), dopplerGrid(:)];
end

function families = build_scheme_families(N, V, c2Base, localDelta, phaseOffset)
    rationalCandidates = afdm.chirp.gps_candidate_set_variant(N, 'rational');
    offsetCandidates = afdm.chirp.gps_candidate_set_variant( ...
        N, 'phase_offset', struct('phase_offset', phaseOffset));
    localCandidates = c2Base + repmat([0, -localDelta, localDelta], N, 1);

    families = repmat(struct( ...
        'name', "", 'patterns', [], 'build_c2', []), 4, 1);
    families(1).name = "baseline";
    families(1).patterns = ones(1, V);
    families(1).build_c2 = @(pattern) c2Base * ones(N, 1);
    families(2).name = "GPS-rational";
    families(2).patterns = enumerate_patterns(2, V);
    families(2).build_c2 = @(pattern) select_grouped_c2( ...
        rationalCandidates, N, V, pattern);
    families(3).name = "GPS-offset";
    families(3).patterns = enumerate_patterns(2, V);
    families(3).build_c2 = @(pattern) select_grouped_c2( ...
        offsetCandidates, N, V, pattern);
    families(4).name = "local-proposed";
    families(4).patterns = enumerate_patterns(3, V);
    families(4).build_c2 = @(pattern) select_grouped_c2( ...
        localCandidates, N, V, pattern);
end

function correlations = build_pair_correlations(N, c1, c2, support, eventSet)
    numSupports = size(support, 1);
    responses = cell(numSupports, 1);
    for supportIdx = 1:numSupports
        pathMatrix = afdm.analysis.build_path_matrix( ...
            N, c1, c2, support(supportIdx, 1), support(supportIdx, 2));
        responses{supportIdx} = pathMatrix * eventSet;
    end

    correlations = cell(numSupports, numSupports);
    for firstIdx = 1:numSupports
        correlations{firstIdx, firstIdx} = ones(1, size(eventSet, 2));
        for secondIdx = firstIdx+1:numSupports
            values = sum(conj(responses{firstIdx}) .* responses{secondIdx}, 1);
            correlations{firstIdx, secondIdx} = values;
            correlations{secondIdx, firstIdx} = conj(values);
        end
    end
end

function audit = audit_pairs(correlations, tolerance)
    numSupports = size(correlations, 1);
    rankLossCount = 0;
    maxCorrelation = 0;
    for firstIdx = 1:numSupports-1
        for secondIdx = firstIdx+1:numSupports
            magnitude = abs(correlations{firstIdx, secondIdx});
            rankLossCount = rankLossCount + sum(magnitude >= 1 - tolerance);
            maxCorrelation = max(maxCorrelation, max(magnitude));
        end
    end
    audit.rank_loss_count = rankLossCount;
    audit.max_correlation = maxCorrelation;
end

function collective = audit_triples( ...
        correlations, triples, eventSet, support, scheme, pattern, tolerance)
    collective.count = 0;
    collective.min_determinant = Inf;
    collective.min_lambda = Inf;
    collective.findings = repmat(empty_finding(), 0, 1);

    for tripleIdx = 1:size(triples, 1)
        indices = triples(tripleIdx, :);
        c12 = correlations{indices(1), indices(2)};
        c13 = correlations{indices(1), indices(3)};
        c23 = correlations{indices(2), indices(3)};
        determinant = real(1 + 2 * c12 .* c23 .* conj(c13) ...
            - abs(c12).^2 - abs(c13).^2 - abs(c23).^2);
        pairwiseIndependent = abs(c12) < 1 - tolerance & ...
            abs(c13) < 1 - tolerance & abs(c23) < 1 - tolerance;
        eligible = find(pairwiseIndependent);
        if isempty(eligible)
            continue;
        end

        [minDeterminant, localIdx] = min(determinant(eligible));
        eventIdx = eligible(localIdx);
        gram = [1, c12(eventIdx), c13(eventIdx); ...
            conj(c12(eventIdx)), 1, c23(eventIdx); ...
            conj(c13(eventIdx)), conj(c23(eventIdx)), 1];
        minLambda = max(min(real(eig(gram))), 0);
        if minDeterminant < collective.min_determinant
            collective.min_determinant = minDeterminant;
            collective.min_lambda = minLambda;
        end

        foundEvents = find(pairwiseIndependent & determinant < tolerance);
        for findingIdx = 1:numel(foundEvents)
            selectedEvent = foundEvents(findingIdx);
            gram = [1, c12(selectedEvent), c13(selectedEvent); ...
                conj(c12(selectedEvent)), 1, c23(selectedEvent); ...
                conj(c13(selectedEvent)), conj(c23(selectedEvent)), 1];
            eigenvalues = sort(real(eig(gram)), 'ascend');
            finding = empty_finding();
            finding.scheme = scheme;
            finding.pattern = pattern;
            finding.support_1 = support_label(support(indices(1), :));
            finding.support_2 = support_label(support(indices(2), :));
            finding.support_3 = support_label(support(indices(3), :));
            finding.event_index = selectedEvent;
            finding.event_weight = sum(abs(eventSet(:, selectedEvent)) > tolerance);
            finding.gram_determinant = determinant(selectedEvent);
            finding.min_eigenvalue = max(eigenvalues(1), 0);
            finding.max_pair_correlation = max( ...
                [abs(c12(selectedEvent)), abs(c13(selectedEvent)), ...
                abs(c23(selectedEvent))]);
            collective.findings(end + 1, 1) = finding;
        end
        collective.count = collective.count + numel(foundEvents);
    end
end

function summary = summarize_families(patternSummary)
    schemes = unique(patternSummary.scheme, 'stable');
    rows = repmat(struct( ...
        'scheme', "", 'num_patterns', NaN, 'pairwise_safe_patterns', NaN, ...
        'patterns_with_collective_only_loss', NaN, ...
        'collective_only_events', NaN, 'best_min_collective_lambda', NaN, ...
        'worst_min_collective_lambda', NaN), numel(schemes), 1);

    for schemeIdx = 1:numel(schemes)
        selected = patternSummary.scheme == schemes(schemeIdx);
        safe = selected & patternSummary.pairwise_safe;
        rows(schemeIdx).scheme = schemes(schemeIdx);
        rows(schemeIdx).num_patterns = sum(selected);
        rows(schemeIdx).pairwise_safe_patterns = sum(safe);
        rows(schemeIdx).patterns_with_collective_only_loss = sum( ...
            safe & patternSummary.collective_only_events > 0);
        rows(schemeIdx).collective_only_events = sum( ...
            patternSummary.collective_only_events(safe), 'omitnan');
        if any(safe)
            rows(schemeIdx).best_min_collective_lambda = max( ...
                patternSummary.min_collective_lambda(safe));
            rows(schemeIdx).worst_min_collective_lambda = min( ...
                patternSummary.min_collective_lambda(safe));
        end
    end
    summary = struct2table(rows);
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

function c2 = select_grouped_c2(candidateSet, N, V, pattern)
    groupIndex = afdm.chirp.group_index(N, V, 'contiguous');
    c2 = zeros(N, 1);
    for groupIdx = 1:V
        selected = groupIndex == groupIdx;
        c2(selected) = candidateSet(selected, pattern(groupIdx));
    end
end

function label = support_label(support)
    label = sprintf('(%d,%d)', support(1), support(2));
end

function plot_results(summary, familySummary, outputDir, outputTag)
    fig = figure('Name', 'Collective-only rank search', 'Color', 'w');
    tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    bar(categorical(familySummary.scheme), ...
        [familySummary.num_patterns, familySummary.pairwise_safe_patterns]);
    grid on;
    ylabel('number of patterns');
    legend({'all', 'pairwise safe'}, 'Location', 'best');
    title('Pairwise screening');

    nexttile;
    safe = summary.pairwise_safe;
    values = summary.min_collective_lambda(safe);
    groups = categorical(summary.scheme(safe));
    boxchart(groups, max(values, eps));
    set(gca, 'YScale', 'log');
    grid on;
    ylabel('minimum collective \lambda_{min}');
    title('Pairwise-safe P=3 conditioning');

    exportgraphics(fig, ...
        fullfile(outputDir, [outputTag '.png']), ...
        'Resolution', 200);
    close(fig);
end

function row = empty_row()
    row = struct( ...
        'scheme', "", ...
        'pattern', "", ...
        'num_events', NaN, ...
        'num_support_pairs', NaN, ...
        'pairwise_rank_loss_events', NaN, ...
        'max_pair_correlation', NaN, ...
        'pairwise_safe', false, ...
        'collective_only_events', NaN, ...
        'min_collective_determinant', NaN, ...
        'min_collective_lambda', NaN);
end

function finding = empty_finding()
    finding = struct( ...
        'scheme', "", ...
        'pattern', "", ...
        'support_1', "", ...
        'support_2', "", ...
        'support_3', "", ...
        'event_index', NaN, ...
        'event_weight', NaN, ...
        'gram_determinant', NaN, ...
        'min_eigenvalue', NaN, ...
        'max_pair_correlation', NaN);
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = defaultValue;
    end
end
