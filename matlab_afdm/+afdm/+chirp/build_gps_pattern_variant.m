function [c2m, d, metadata] = build_gps_pattern_variant( ...
        N, V, pattern, variant, options)
%BUILD_GPS_PATTERN_VARIANT Build one grouped GPS realization.

    if nargin < 5
        options = struct();
    end
    if numel(pattern) ~= V
        error('pattern length must equal V.');
    end

    [candidateSet, metadata] = afdm.chirp.gps_candidate_set_variant( ...
        N, variant, options);
    groupIndex = afdm.chirp.group_index(N, V, 'contiguous');
    c2m = zeros(N, 1);
    for groupIdx = 1:V
        candidateIdx = pattern(groupIdx);
        if candidateIdx < 1 || candidateIdx > 2
            error('GPS pattern entries must be 1 or 2.');
        end
        selected = groupIndex == groupIdx;
        c2m(selected) = candidateSet(selected, candidateIdx);
    end

    m = (0:N-1).';
    d = exp(-1i * 2 * pi * c2m .* (m.^2));
    metadata.pattern = pattern;
end
