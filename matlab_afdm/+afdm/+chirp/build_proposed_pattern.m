function [c2m, d] = build_proposed_pattern(N, V, pattern, c2Base, delta)
%BUILD_PROPOSED_PATTERN Build proposed grouped small-perturbation c2_m.
%   pattern(v)=1 selects c2Base; 2 selects c2Base-delta; 3 selects
%   c2Base+delta.

    validateattributes(N, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(V, {'numeric'}, {'scalar', 'integer', 'positive'});
    if numel(pattern) ~= V
        error('pattern length must equal V.');
    end

    groupLength = N / V;
    if abs(groupLength - round(groupLength)) > 0
        error('N must be divisible by V.');
    end

    candidates = [c2Base, c2Base - delta, c2Base + delta];
    c2m = zeros(N, 1);
    for groupId = 1:V
        candidateId = pattern(groupId);
        if candidateId < 1 || candidateId > numel(candidates)
            error('Proposed pattern entries must be 1, 2, or 3.');
        end
        c2m((groupId - 1) * groupLength + 1:groupId * groupLength) = candidates(candidateId);
    end

    m = (0:N-1).';
    d = exp(-1i * 2 * pi * c2m .* (m.^2));
end
