function hPath = build_path_matrix(N, c1, c2m, delay, doppler)
%BUILD_PATH_MATRIX Build a unit-gain single-path AFDM channel matrix.
%   Supports scalar uniform c2 or a per-subcarrier c2_m vector.

    if isscalar(c2m)
        c2Arg = c2m;
    else
        if numel(c2m) ~= N
            error('c2m must be scalar or length N.');
        end
        c2Arg = c2m(:);
    end

    hPath = afdm.rx.estimate_effective_channel(N, c1, c2Arg, 1, delay, doppler / N);
end
