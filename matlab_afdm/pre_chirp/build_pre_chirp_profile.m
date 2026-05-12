function profile = build_pre_chirp_profile(scheme, numSubcarriers, params)
%BUILD_PRE_CHIRP_PROFILE Compatibility wrapper.
%   Prefer afdm.chirp.build_profile(scheme, numSubcarriers, params).

    if nargin < 3
        params = struct();
    end
    profile = afdm.chirp.build_profile(scheme, numSubcarriers, params);
end
