function results = run_n64_collective_rank_construction(options)
%RUN_N64_COLLECTIVE_RANK_CONSTRUCTION Search four-cycle P=4 failures at N=64.

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
    phaseOffset = get_option(options, 'phase_offset', 0.3);
    tolerance = get_option(options, 'tolerance', 1e-10);
    support = get_option(options, 'support', [0, 0; 2, 2; 5, -3; 7, -1]);

    if N ~= 64 || V ~= 4
        error('The default construction currently requires N=64 and V=4.');
    end
    if size(support, 1) ~= 4
        error('support must contain four [delay,Doppler] rows.');
    end

    [eventSet, eventAmplitude, eventOffset] = build_four_block_events(N);
    families = build_scheme_families(N, V, c2Base, localDelta, phaseOffset);
    totalPatterns = sum(arrayfun(@(x) size(x.patterns, 1), families));
    rows = repmat(empty_row(), totalPatterns, 1);
    findings = repmat(empty_finding(), 0, 1);
    rowIdx = 0;

    fprintf('========== N=64 collective-rank construction ==========\n');
    fprintf('support=%s, events=%d, patterns=%d\n', ...
        mat2str(support), size(eventSet, 2), totalPatterns);

    for familyIdx = 1:numel(families)
        family = families(familyIdx);
        fprintf('Scanning %s: %d patterns\n', ...
            family.name, size(family.patterns, 1));
        for patternIdx = 1:size(family.patterns, 1)
            rowIdx = rowIdx + 1;
            pattern = family.patterns(patternIdx, :);
            c2 = family.build_c2(pattern);
            responses = build_responses(N, c1, c2, support, eventSet);
            audit = audit_four_path_events(responses, tolerance);

            row = empty_row();
            row.scheme = family.name;
            row.pattern = join(string(pattern), '-');
            row.num_events = size(eventSet, 2);
            row.pairwise_rank_loss_events = sum(audit.pairwise_dependent);
            row.collective_only_events = sum(audit.collective_only);
            row.min_exact_lambda = min(audit.min_lambda);
            if any(~audit.pairwise_dependent)
                safeEvents = find(~audit.pairwise_dependent);
                [row.min_pairwise_safe_lambda, localWorstIdx] = min( ...
                    audit.min_lambda(safeEvents));
                worstEventIdx = safeEvents(localWorstIdx);
                row.worst_pairwise_safe_event_index = worstEventIdx;
                row.worst_block_offset = eventOffset(worstEventIdx);
                row.worst_block_amplitudes = join( ...
                    string(eventAmplitude(worstEventIdx, :)), ',');
                row.max_pair_correlation_at_worst = ...
                    audit.max_pair_correlation(worstEventIdx);
            end
            if any(audit.collective_only)
                row.max_pair_correlation_at_collective = max( ...
                    audit.max_pair_correlation(audit.collective_only));
            end
            rows(rowIdx) = row;

            selectedEvents = find(audit.collective_only);
            for selectedIdx = 1:numel(selectedEvents)
                eventIdx = selectedEvents(selectedIdx);
                finding = empty_finding();
                finding.scheme = family.name;
                finding.pattern = row.pattern;
                finding.event_index = eventIdx;
                finding.block_offset = eventOffset(eventIdx);
                finding.block_amplitudes = join(string(eventAmplitude(eventIdx, :)), ',');
                finding.event_weight = sum(abs(eventSet(:, eventIdx)) > tolerance);
                finding.min_eigenvalue = audit.min_lambda(eventIdx);
                finding.max_pair_correlation = audit.max_pair_correlation(eventIdx);
                findings(end + 1, 1) = finding; %#ok<AGROW>
            end
        end
    end

    summary = struct2table(rows);
    findingTable = struct2table(findings);
    familySummary = summarize_families(summary);
    results.config = struct( ...
        'N', N, 'V', V, 'alpha_max', alphaMax, 'c1', c1, ...
        'c2_base', c2Base, 'local_delta', localDelta, ...
        'phase_offset', phaseOffset, 'tolerance', tolerance);
    results.support = support;
    results.event_set = eventSet;
    results.event_amplitude = eventAmplitude;
    results.event_offset = eventOffset;
    results.summary = summary;
    results.family_summary = familySummary;
    results.findings = findingTable;

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    save(fullfile(outputDir, 'n64_collective_rank_construction.mat'), 'results');
    writetable(summary, ...
        fullfile(outputDir, 'n64_collective_rank_construction.csv'));
    writetable(familySummary, ...
        fullfile(outputDir, 'n64_collective_rank_construction_families.csv'));
    if ~isempty(findingTable)
        writetable(findingTable, ...
            fullfile(outputDir, 'n64_collective_rank_construction_findings.csv'));
    end
    plot_results(summary, familySummary, outputDir);

    disp(familySummary);
    if isempty(findingTable)
        fprintf('No N=64 collective-only event found in the structured set.\n');
    else
        fprintf('Found %d N=64 collective-only events.\n', height(findingTable));
        disp(findingTable(1:min(20, height(findingTable)), :));
    end
end

function [events, amplitudes, offsets] = build_four_block_events(N)
    blockLength = N / 4;
    amplitudePatterns = enumerate_ternary_vectors(4);
    numEvents = size(amplitudePatterns, 1) * blockLength;
    events = zeros(N, numEvents);
    amplitudes = zeros(numEvents, 4);
    offsets = zeros(numEvents, 1);
    cursor = 0;

    for offset = 1:blockLength
        positions = offset + (0:3) * blockLength;
        for patternIdx = 1:size(amplitudePatterns, 1)
            cursor = cursor + 1;
            amplitudes(cursor, :) = amplitudePatterns(patternIdx, :);
            offsets(cursor) = offset;
            events(positions, cursor) = 1i * amplitudePatterns(patternIdx, :).';
        end
    end
    events = events ./ vecnorm(events, 2, 1);
end

function patterns = enumerate_ternary_vectors(lengthValue)
    patterns = zeros(3^lengthValue - 1, lengthValue);
    alphabet = [-1, 0, 1];
    cursor = 0;
    for value = 0:3^lengthValue-1
        digits = zeros(1, lengthValue);
        state = value;
        for idx = 1:lengthValue
            digits(idx) = mod(state, 3) + 1;
            state = floor(state / 3);
        end
        vector = alphabet(digits);
        if all(vector == 0)
            continue;
        end
        cursor = cursor + 1;
        patterns(cursor, :) = vector;
    end
    patterns = patterns(1:cursor, :);
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

function responses = build_responses(N, c1, c2, support, eventSet)
    responses = cell(4, 1);
    for pathIdx = 1:4
        pathMatrix = afdm.analysis.build_path_matrix( ...
            N, c1, c2, support(pathIdx, 1), support(pathIdx, 2));
        responses{pathIdx} = pathMatrix * eventSet;
    end
end

function audit = audit_four_path_events(responses, tolerance)
    numEvents = size(responses{1}, 2);
    pairwiseDependent = false(numEvents, 1);
    maxPairCorrelation = zeros(numEvents, 1);
    gram = repmat(eye(4), 1, 1, numEvents);

    for firstIdx = 1:3
        for secondIdx = firstIdx+1:4
            correlation = sum( ...
                conj(responses{firstIdx}) .* responses{secondIdx}, 1).';
            gram(firstIdx, secondIdx, :) = correlation;
            gram(secondIdx, firstIdx, :) = conj(correlation);
            magnitude = abs(correlation);
            maxPairCorrelation = max(maxPairCorrelation, magnitude);
            pairwiseDependent = pairwiseDependent | magnitude >= 1 - tolerance;
        end
    end

    minLambda = zeros(numEvents, 1);
    for eventIdx = 1:numEvents
        minLambda(eventIdx) = max(min(real(eig(gram(:, :, eventIdx)))), 0);
    end
    audit.pairwise_dependent = pairwiseDependent;
    audit.max_pair_correlation = maxPairCorrelation;
    audit.min_lambda = minLambda;
    audit.collective_only = minLambda < tolerance & ~pairwiseDependent;
end

function summary = summarize_families(patternSummary)
    schemes = unique(patternSummary.scheme, 'stable');
    rows = repmat(struct( ...
        'scheme', "", 'num_patterns', NaN, ...
        'patterns_with_collective_only_loss', NaN, ...
        'collective_only_events', NaN, ...
        'best_min_pairwise_safe_lambda', NaN, ...
        'worst_min_pairwise_safe_lambda', NaN), numel(schemes), 1);

    for schemeIdx = 1:numel(schemes)
        selected = patternSummary.scheme == schemes(schemeIdx);
        rows(schemeIdx).scheme = schemes(schemeIdx);
        rows(schemeIdx).num_patterns = sum(selected);
        rows(schemeIdx).patterns_with_collective_only_loss = sum( ...
            selected & patternSummary.collective_only_events > 0);
        rows(schemeIdx).collective_only_events = sum( ...
            patternSummary.collective_only_events(selected));
        rows(schemeIdx).best_min_pairwise_safe_lambda = max( ...
            patternSummary.min_pairwise_safe_lambda(selected));
        rows(schemeIdx).worst_min_pairwise_safe_lambda = min( ...
            patternSummary.min_pairwise_safe_lambda(selected));
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

function plot_results(summary, familySummary, outputDir)
    fig = figure('Name', 'N64 collective rank construction', 'Color', 'w');
    tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    counts = familySummary.patterns_with_collective_only_loss;
    bar(categorical(familySummary.scheme), counts);
    grid on;
    ylim([0, max(1, max(counts) + 1)]);
    if all(counts == 0)
        text(2.5, 0.5, 'No exact collective-only event found', ...
            'HorizontalAlignment', 'center');
    end
    ylabel('patterns with collective-only loss');
    title('N=64 structured counterexamples');

    nexttile;
    values = max(summary.min_pairwise_safe_lambda, eps);
    boxchart(categorical(summary.scheme), values);
    set(gca, 'YScale', 'log');
    grid on;
    ylabel('minimum P=4 \lambda_{min}');
    title('Pairwise-independent event margin');

    exportgraphics(fig, ...
        fullfile(outputDir, 'n64_collective_rank_construction.png'), ...
        'Resolution', 200);
    close(fig);
end

function row = empty_row()
    row = struct( ...
        'scheme', "", ...
        'pattern', "", ...
        'num_events', NaN, ...
        'pairwise_rank_loss_events', NaN, ...
        'collective_only_events', NaN, ...
        'min_exact_lambda', NaN, ...
        'min_pairwise_safe_lambda', NaN, ...
        'max_pair_correlation_at_collective', NaN, ...
        'worst_pairwise_safe_event_index', NaN, ...
        'worst_block_offset', NaN, ...
        'worst_block_amplitudes', "", ...
        'max_pair_correlation_at_worst', NaN);
end

function finding = empty_finding()
    finding = struct( ...
        'scheme', "", ...
        'pattern', "", ...
        'event_index', NaN, ...
        'block_offset', NaN, ...
        'block_amplitudes', "", ...
        'event_weight', NaN, ...
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
