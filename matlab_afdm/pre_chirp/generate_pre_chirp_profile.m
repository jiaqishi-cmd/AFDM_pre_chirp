function profile = generate_pre_chirp_profile(scheme, numSubcarriers, params)
%GENERATE_PRE_CHIRP_PROFILE Build a pre-chirp coefficient profile.
%   profile = generate_pre_chirp_profile(scheme, numSubcarriers, params)

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
            profile = paper_grouping_profile(numSubcarriers, params);
        case 'proposed_grouping'
            profile = proposed_grouping_profile(numSubcarriers, params);
        otherwise
            error('Unsupported pre-chirp scheme: %s', scheme);
    end
end
