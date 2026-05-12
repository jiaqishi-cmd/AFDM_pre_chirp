function [c2m, d] = build_gps_pattern(N, V, pattern)
%BUILD_GPS_PATTERN Build Yuan GPS W=2 c2_m and d_m for fixed groups.
%   Mathematical index m is 0:N-1. pattern(v)=1 selects +1/(4*m^2);
%   pattern(v)=2 selects -1/(4*m^2). At m=0, c2_m is set to 0.

    validateattributes(N, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(V, {'numeric'}, {'scalar', 'integer', 'positive'});
    if numel(pattern) ~= V
        error('pattern length must equal V.');
    end

    groupLength = N / V;
    if abs(groupLength - round(groupLength)) > 0
        error('N must be divisible by V.');
    end

    c2m = zeros(N, 1);
    m = (0:N-1).';

    for groupId = 1:V
        firstIdx = (groupId - 1) * groupLength + 1;
        lastIdx = groupId * groupLength;
        groupM = m(firstIdx:lastIdx);

        if pattern(groupId) == 1
            signVal = 1;
        elseif pattern(groupId) == 2
            signVal = -1;
        else
            error('GPS pattern entries must be 1 or 2.');
        end

        nonzero = groupM ~= 0;
        c2Group = zeros(numel(groupM), 1);
        c2Group(nonzero) = signVal ./ (4 * groupM(nonzero).^2);
        c2m(firstIdx:lastIdx) = c2Group;
    end

    d = exp(-1i * 2 * pi * c2m .* (m.^2));
end
