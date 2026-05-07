function [c2m, d_m] = build_c2m_gps_pattern(N, V, pattern)
%BUILD_C2M_GPS_PATTERN Build Yuan GPS W=2 c2_m and d_m for fixed groups.
%   Mathematical index m is 0:N-1. pattern(v)=1 selects +1/(4*m^2);
%   pattern(v)=2 selects -1/(4*m^2). At m=0, c2_m is set to 0.

    validateattributes(N, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(V, {'numeric'}, {'scalar', 'integer', 'positive'});
    if numel(pattern) ~= V
        error('pattern length must equal V.');
    end

    M = N / V;
    if abs(M - round(M)) > 0
        error('N must be divisible by V.');
    end

    c2m = zeros(N, 1);
    m = (0:N-1).';

    for group_id = 1:V
        first_idx = (group_id - 1) * M + 1;
        last_idx = group_id * M;
        group_m = m(first_idx:last_idx);

        sign_val = 1;
        if pattern(group_id) == 2
            sign_val = -1;
        elseif pattern(group_id) ~= 1
            error('GPS pattern entries must be 1 or 2.');
        end

        nonzero = group_m ~= 0;
        c2_group = zeros(numel(group_m), 1);
        c2_group(nonzero) = sign_val ./ (4 * group_m(nonzero).^2);
        c2m(first_idx:last_idx) = c2_group;
    end

    d_m = exp(-1i * 2 * pi * c2m .* (m.^2));
end
