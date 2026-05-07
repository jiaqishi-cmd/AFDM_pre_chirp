function H_path = build_H_path_general_c2m(N, c1, c2m, l, alpha)
%BUILD_H_PATH_GENERAL_C2M Build unit-gain single-path AFDM matrix.
%   Supports scalar uniform c2 or per-subcarrier c2_m vector.

    if isscalar(c2m)
        c2_arg = c2m;
    else
        if numel(c2m) ~= N
            error('c2m must be scalar or length N.');
        end
        c2_arg = c2m(:);
    end

    H_path = estimate_effective_channel(N, c1, c2_arg, 1, l, alpha / N);
end
