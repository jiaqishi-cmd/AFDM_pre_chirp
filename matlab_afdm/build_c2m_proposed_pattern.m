function [c2m, d_m] = build_c2m_proposed_pattern(N, V, pattern, c2_base, delta)
%BUILD_C2M_PROPOSED_PATTERN Build proposed grouped small-perturbation c2_m.
%   pattern(v)=1 selects c2_base; 2 selects c2_base-delta; 3 selects
%   c2_base+delta.

    validateattributes(N, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(V, {'numeric'}, {'scalar', 'integer', 'positive'});
    if numel(pattern) ~= V
        error('pattern length must equal V.');
    end

    M = N / V;
    if abs(M - round(M)) > 0
        error('N must be divisible by V.');
    end

    candidates = [c2_base, c2_base - delta, c2_base + delta];
    c2m = zeros(N, 1);
    for group_id = 1:V
        candidate_id = pattern(group_id);
        if candidate_id < 1 || candidate_id > numel(candidates)
            error('Proposed pattern entries must be 1, 2, or 3.');
        end
        c2m((group_id - 1) * M + 1:group_id * M) = candidates(candidate_id);
    end

    m = (0:N-1).';
    d_m = exp(-1i * 2 * pi * c2m .* (m.^2));
end
