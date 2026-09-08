function result = monomial_rotational_cycle_condition( ...
        relativeOperator, rotationalOrder, tolerance)
%MONOMIAL_ROTATIONAL_CYCLE_CONDITION Test exact phase-cycle closure.
%   A cycle is compatible with a Q-fold rotational alphabet when there is
%   an eigenvalue phase phi such that every a_r/exp(j*phi) is a Q-th root
%   of unity and the product of these ratios around the cycle equals one.

    if nargin < 3 || isempty(tolerance)
        tolerance = 1e-10;
    end
    validateattributes(rotationalOrder, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});

    [N, M] = size(relativeOperator);
    if N ~= M
        error('relativeOperator must be square.');
    end

    permutation = zeros(N, 1);
    coefficients = zeros(N, 1);
    for columnIdx = 1:N
        [~, rowIdx] = max(abs(relativeOperator(:, columnIdx)));
        permutation(columnIdx) = rowIdx;
        coefficients(columnIdx) = relativeOperator(rowIdx, columnIdx);
    end
    cycles = extract_cycles(permutation);

    compatible = false(numel(cycles), 1);
    compatiblePhi = NaN(numel(cycles), 1);
    localResidual = Inf(numel(cycles), 1);
    closureResidual = Inf(numel(cycles), 1);

    for cycleIdx = 1:numel(cycles)
        cycle = cycles{cycleIdx};
        beta = angle(coefficients(cycle));
        cycleLength = numel(cycle);
        for rootIdx = 0:rotationalOrder-1
            phi = beta(1) - 2 * pi * rootIdx / rotationalOrder;
            local = max(abs(exp(1i * rotationalOrder * (beta - phi)) - 1));
            closure = abs(exp(1i * (cycleLength * phi - sum(beta))) - 1);
            if local < localResidual(cycleIdx) || ...
                    (abs(local - localResidual(cycleIdx)) <= tolerance && ...
                    closure < closureResidual(cycleIdx))
                localResidual(cycleIdx) = local;
                closureResidual(cycleIdx) = closure;
            end
            if local <= tolerance && closure <= tolerance
                compatible(cycleIdx) = true;
                compatiblePhi(cycleIdx) = phi;
                localResidual(cycleIdx) = local;
                closureResidual(cycleIdx) = closure;
                break;
            end
        end
    end

    result.cycles = cycles;
    result.compatible_cycles = compatible;
    result.compatible_phi = compatiblePhi;
    result.local_residual = localResidual;
    result.closure_residual = closureResidual;
    result.num_compatible_cycles = sum(compatible);
    result.predicts_compatible_cycle = any(compatible);
    result.rotational_order = rotationalOrder;
    result.tolerance = tolerance;
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
