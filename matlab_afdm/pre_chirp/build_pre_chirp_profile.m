function profile = build_pre_chirp_profile(scheme, numSubcarriers, params)
%BUILD_PRE_CHIRP_PROFILE Build a static pre-chirp profile definition.
%   The returned profile contains the default c2, grouping, and candidates.
%   Frame-dependent greedy selection is handled separately.

    if nargin < 1 || isempty(scheme)
        scheme = 'baseline';
    end
    if nargin < 3
        params = struct();
    end

    switch lower(scheme)
        case 'baseline'
            profile = baseline_profile(numSubcarriers, params);
        case 'paper_grouping'
            profile = gps_profile(numSubcarriers, params);
        case 'proposed_grouping'
            profile = proposed_profile(numSubcarriers, params);
        otherwise
            error('Unsupported pre-chirp scheme: %s', scheme);
    end
end
