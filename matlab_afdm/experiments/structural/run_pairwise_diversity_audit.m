function results = run_pairwise_diversity_audit(options)
%RUN_PAIRWISE_DIVERSITY_AUDIT Audit pairwise diversity over pre-chirp modes.
%   The default N=8 BPSK run exhausts every nonzero constellation
%   difference vector and every V=4 GPS/proposed group pattern.

    if nargin < 1
        options = struct();
    end

    rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    setup_paths(rootDir);

    N = get_option(options, 'N', 8);
    V = get_option(options, 'V', 4);
    MMod = get_option(options, 'M_mod', 2);
    modType = get_option(options, 'mod_type', 'psk');
    maxWeight = get_option(options, 'max_weight', default_max_weight(N, MMod));
    alphaMax = get_option(options, 'alpha_max', 2);
    c1 = get_option(options, 'c1', (2 * alphaMax + 1) / (2 * N));
    c2Base = get_option(options, 'c2_base', sqrt(2) / (10 * N));
    proposedDelta = get_option(options, 'proposed_delta', c2Base / 16);
    delays = get_option(options, 'delays', [0, 2]);
    dopplers = get_option(options, 'dopplers', [0, 2]);

    if numel(delays) ~= numel(dopplers)
        error('delays and dopplers must have the same length.');
    end
    if mod(N, V) ~= 0
        error('N must be divisible by V.');
    end

    deltaSet = afdm.delta.generate_sparse_difference_set( ...
        N, MMod, modType, maxWeight);
    deltaNorm = sqrt(deltaSet.norm2);
    deltaNormalized = deltaSet.delta ./ deltaNorm;

    gpsPatterns = get_option(options, 'gps_patterns', enumerate_patterns(2, V));
    proposedPatterns = get_option(options, ...
        'proposed_patterns', enumerate_patterns(3, V));
    if size(gpsPatterns, 2) ~= V || size(proposedPatterns, 2) ~= V
        error('Every supplied pattern must have V entries.');
    end

    rows = repmat(empty_row(), ...
        1 + size(gpsPatterns, 1) + size(proposedPatterns, 1), 1);
    cursor = 0;

    cursor = cursor + 1;
    rows(cursor) = audit_pattern('baseline', zeros(1, V), c2Base, ...
        delays, dopplers, N, c1, deltaNormalized, deltaSet);

    for patternIdx = 1:size(gpsPatterns, 1)
        pattern = gpsPatterns(patternIdx, :);
        c2 = afdm.chirp.build_gps_pattern(N, V, pattern);
        cursor = cursor + 1;
        rows(cursor) = audit_pattern('GPS', pattern, c2, ...
            delays, dopplers, N, c1, deltaNormalized, deltaSet);
    end

    for patternIdx = 1:size(proposedPatterns, 1)
        pattern = proposedPatterns(patternIdx, :);
        c2 = afdm.chirp.build_proposed_pattern( ...
            N, V, pattern, c2Base, proposedDelta);
        cursor = cursor + 1;
        rows(cursor) = audit_pattern('proposed', pattern, c2, ...
            delays, dopplers, N, c1, deltaNormalized, deltaSet);
    end

    auditTable = struct2table(rows);
    results.config = struct( ...
        'N', N, 'V', V, 'M_mod', MMod, 'mod_type', modType, ...
        'max_weight', maxWeight, 'alpha_max', alphaMax, 'c1', c1, ...
        'c2_base', c2Base, 'proposed_delta', proposedDelta, ...
        'delays', delays, 'dopplers', dopplers);
    results.delta_set = deltaSet;
    results.audit_table = auditTable;

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    outputTag = sprintf('pairwise_diversity_audit_N%d_M%d_w%d', ...
        N, MMod, maxWeight);
    save(fullfile(outputDir, [outputTag '.mat']), 'results');
    writetable(auditTable, fullfile(outputDir, [outputTag '.csv']));

    print_summary(auditTable, size(deltaSet.delta, 2), numel(delays));
end

function row = audit_pattern(scheme, pattern, c2, delays, dopplers, ...
        N, c1, deltaNormalized, deltaSet)
    numPaths = numel(delays);
    pathMatrices = cell(1, numPaths);
    for pathIdx = 1:numPaths
        pathMatrices{pathIdx} = afdm.analysis.build_path_matrix( ...
            N, c1, c2, delays(pathIdx), dopplers(pathIdx));
    end

    minRank = numPaths;
    minSigma = Inf;
    minGramDet = Inf;
    rankLossCount = 0;
    worstEventIndex = 1;

    for eventIdx = 1:size(deltaNormalized, 2)
        psi = afdm.analysis.build_pairwise_path_response( ...
            pathMatrices, deltaNormalized(:, eventIdx));
        metrics = afdm.analysis.pairwise_diversity_metrics(psi);
        if metrics.rank < numPaths
            rankLossCount = rankLossCount + 1;
        end
        if metrics.min_singular_value < minSigma
            minSigma = metrics.min_singular_value;
            worstEventIndex = eventIdx;
        end
        minRank = min(minRank, metrics.rank);
        minGramDet = min(minGramDet, metrics.gram_determinant);
    end

    row = empty_row();
    row.scheme = string(scheme);
    row.pattern = pattern_label(pattern, scheme);
    row.effective_phase_count = count_effective_phases(c2, N);
    row.min_rank = minRank;
    row.rank_loss_events = rankLossCount;
    row.rank_loss_fraction = rankLossCount / size(deltaNormalized, 2);
    row.min_sigma = minSigma;
    row.min_gram_det = minGramDet;
    row.worst_event_index = worstEventIndex;
    row.worst_event_weight = deltaSet.weight(worstEventIndex);
end

function patterns = enumerate_patterns(numCandidates, numGroups)
    numPatterns = numCandidates^numGroups;
    patterns = ones(numPatterns, numGroups);
    for patternIdx = 0:numPatterns-1
        value = patternIdx;
        for groupIdx = 1:numGroups
            patterns(patternIdx + 1, groupIdx) = mod(value, numCandidates) + 1;
            value = floor(value / numCandidates);
        end
    end
end

function count = count_effective_phases(c2, N)
    if isscalar(c2)
        c2 = repmat(c2, N, 1);
    end
    m = (0:N-1).';
    phase = exp(1i * 2 * pi * c2(:) .* (m.^2));
    rounded = round(real(phase), 10) + 1i * round(imag(phase), 10);
    count = numel(unique(rounded));
end

function label = pattern_label(pattern, scheme)
    if strcmpi(scheme, 'baseline')
        label = "nominal";
    else
        label = join(string(pattern), '-');
    end
end

function print_summary(T, numEvents, numPaths)
    fprintf('Pairwise diversity audit: %d error events, %d paths\n', ...
        numEvents, numPaths);
    schemes = unique(T.scheme, 'stable');
    for idx = 1:numel(schemes)
        subset = T(T.scheme == schemes(idx), :);
        fprintf('%-9s patterns=%3d min-rank=%d rank-loss-patterns=%d min-sigma=%.3e\n', ...
            char(schemes(idx)), height(subset), min(subset.min_rank), ...
            sum(subset.rank_loss_events > 0), min(subset.min_sigma));
    end

    [~, order] = sort(T.min_sigma, 'ascend');
    disp(T(order(1:min(12, height(T))), ...
        {'scheme', 'pattern', 'effective_phase_count', 'min_rank', ...
        'rank_loss_events', 'min_sigma', 'min_gram_det', ...
        'worst_event_weight'}));
end

function row = empty_row()
    row = struct( ...
        'scheme', "", ...
        'pattern', "", ...
        'effective_phase_count', NaN, ...
        'min_rank', NaN, ...
        'rank_loss_events', NaN, ...
        'rank_loss_fraction', NaN, ...
        'min_sigma', NaN, ...
        'min_gram_det', NaN, ...
        'worst_event_index', NaN, ...
        'worst_event_weight', NaN);
end

function value = default_max_weight(N, MMod)
    if MMod == 2 && N <= 10
        value = N;
    else
        value = min(N, 2);
    end
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = defaultValue;
    end
end
