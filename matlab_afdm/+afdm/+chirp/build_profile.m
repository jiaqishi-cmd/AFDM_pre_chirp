function profile = build_profile(scheme, numSubcarriers, params)
%BUILD_PROFILE Build a static pre-chirp profile definition.
%   Frame-dependent greedy selection is handled separately.

    if nargin < 1 || isempty(scheme)
        scheme = 'baseline';
    end
    if nargin < 3
        params = struct();
    end

    switch lower(scheme)
        case 'baseline'
            profile = afdm.chirp.baseline_profile(numSubcarriers, params);
        case 'paper_grouping'
            profile = afdm.chirp.gps_profile(numSubcarriers, params);
        case 'proposed_grouping'
            profile = afdm.chirp.proposed_profile(numSubcarriers, params);
        otherwise
            error('Unsupported pre-chirp scheme: %s', scheme);
    end
end
