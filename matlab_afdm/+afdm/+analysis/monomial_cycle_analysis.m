function analysis = monomial_cycle_analysis(relativeOperator, differenceAlphabet, tolerance)
%MONOMIAL_CYCLE_ANALYSIS Find constellation-compatible monomial eigenmodes.
%   relativeOperator maps one unit-gain path response into another. A
%   compatible eigenvector delta satisfies H2*delta=lambda*H1*delta and
%   therefore makes the two columns of the pairwise path-response matrix
%   linearly dependent.

    if nargin < 3 || isempty(tolerance)
        tolerance = 1e-9;
    end

    [N, M] = size(relativeOperator);
    if N ~= M
        error('relativeOperator must be square.');
    end
    differenceAlphabet = differenceAlphabet(:);
    differenceAlphabet = differenceAlphabet(abs(differenceAlphabet) > tolerance);
    if isempty(differenceAlphabet)
        error('differenceAlphabet must contain a nonzero value.');
    end

    permutation = zeros(N, 1);
    dominantValues = zeros(N, 1);
    leakagePower = zeros(N, 1);
    for columnIdx = 1:N
        [~, rowIdx] = max(abs(relativeOperator(:, columnIdx)));
        permutation(columnIdx) = rowIdx;
        dominantValues(columnIdx) = relativeOperator(rowIdx, columnIdx);
        residual = relativeOperator(:, columnIdx);
        residual(rowIdx) = 0;
        leakagePower(columnIdx) = sum(abs(residual).^2);
    end
    if numel(unique(permutation)) ~= N
        error('Dominant entries do not define a permutation.');
    end

    cycles = extract_cycles(permutation);
    compatibleDeltas = {};
    compatibleEigenvalues = [];
    compatibleCycleIds = [];
    minAlphabetMismatch = Inf;
    closestDelta = [];
    closestEigenvalue = NaN;
    closestCycleId = NaN;
    cycleRows = repmat(empty_cycle_row(), numel(cycles), 1);

    for cycleIdx = 1:numel(cycles)
        cycle = cycles{cycleIdx};
        subOperator = relativeOperator(cycle, cycle);
        [vectors, values] = eig(subOperator, 'vector');
        compatibleCount = 0;

        for eigenIdx = 1:size(vectors, 2)
            [isCompatible, scaledVector, mismatch, quantizedVector] = ...
                match_difference_alphabet( ...
                vectors(:, eigenIdx), differenceAlphabet, tolerance);
            if mismatch < minAlphabetMismatch
                minAlphabetMismatch = mismatch;
                closestDelta = zeros(N, 1, 'like', relativeOperator);
                closestDelta(cycle) = quantizedVector;
                closestEigenvalue = values(eigenIdx);
                closestCycleId = cycleIdx;
            end
            if isCompatible
                delta = zeros(N, 1, 'like', relativeOperator);
                delta(cycle) = scaledVector;
                compatibleDeltas{end + 1, 1} = delta; %#ok<AGROW>
                compatibleEigenvalues(end + 1, 1) = values(eigenIdx); %#ok<AGROW>
                compatibleCycleIds(end + 1, 1) = cycleIdx; %#ok<AGROW>
                compatibleCount = compatibleCount + 1;
            end
        end

        row = empty_cycle_row();
        row.cycle_id = cycleIdx;
        row.cycle_length = numel(cycle);
        row.indices_zero_based = join(string(cycle - 1), '-');
        row.compatible_eigenvectors = compatibleCount;
        cycleRows(cycleIdx) = row;
    end

    shiftValues = mod(permutation - (1:N).', N);
    uniqueShifts = unique(shiftValues);
    if numel(uniqueShifts) == 1
        constantShift = uniqueShifts;
    else
        constantShift = NaN;
    end

    cycleTable = struct2table(cycleRows);
    analysis.permutation = permutation;
    analysis.dominant_values = dominantValues;
    analysis.max_monomial_leakage_power = max(leakagePower);
    analysis.constant_shift = constantShift;
    analysis.shift_values = shiftValues;
    analysis.cycles = cycles;
    analysis.cycle_table = cycleTable;
    analysis.compatible_deltas = compatibleDeltas;
    analysis.compatible_eigenvalues = compatibleEigenvalues;
    analysis.compatible_cycle_ids = compatibleCycleIds;
    analysis.num_compatible_eigenvectors = numel(compatibleDeltas);
    analysis.min_difference_alphabet_mismatch = minAlphabetMismatch;
    analysis.closest_delta = closestDelta;
    analysis.closest_eigenvalue = closestEigenvalue;
    analysis.closest_cycle_id = closestCycleId;
    if isempty(compatibleCycleIds)
        analysis.min_compatible_weight = NaN;
    else
        lengths = cycleTable.cycle_length(compatibleCycleIds);
        analysis.min_compatible_weight = min(lengths);
    end
    analysis.difference_alphabet = differenceAlphabet;
    analysis.tolerance = tolerance;
end

function cycles = extract_cycles(permutation)
    N = numel(permutation);
    visited = false(N, 1);
    cycles = {};
    for startIdx = 1:N
        if visited(startIdx)
            continue;
        end
        cycle = [];
        current = startIdx;
        while ~visited(current)
            cycle(end + 1) = current; %#ok<AGROW>
            visited(current) = true;
            current = permutation(current);
        end
        cycles{end + 1, 1} = cycle; %#ok<AGROW>
    end
end

function [compatible, scaledVector, bestMismatch, bestQuantizedVector] = ...
        match_difference_alphabet( ...
        vector, differenceAlphabet, tolerance)
    compatible = false;
    scaledVector = [];
    bestMismatch = Inf;
    bestQuantizedVector = [];
    [~, referenceIdx] = max(abs(vector));
    referenceValue = vector(referenceIdx);

    for alphabetIdx = 1:numel(differenceAlphabet)
        scale = differenceAlphabet(alphabetIdx) / referenceValue;
        candidate = scale * vector;
        distance = abs(candidate - differenceAlphabet.');
        [minDistance, nearestIdx] = min(distance, [], 2);
        mismatch = max(minDistance) / max(abs(differenceAlphabet));
        if mismatch < bestMismatch
            bestMismatch = mismatch;
            bestQuantizedVector = differenceAlphabet(nearestIdx);
        end
        threshold = tolerance * max(1, max(abs(candidate)));
        if all(minDistance <= threshold)
            compatible = true;
            scaledVector = candidate;
            return;
        end
    end
end

function row = empty_cycle_row()
    row = struct( ...
        'cycle_id', NaN, ...
        'cycle_length', NaN, ...
        'indices_zero_based', "", ...
        'compatible_eigenvectors', NaN);
end
