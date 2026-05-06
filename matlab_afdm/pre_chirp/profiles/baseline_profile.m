function profile = baseline_profile(numSubcarriers, params)
%BASELINE_PROFILE Use one scalar pre-chirp coefficient for all subcarriers.

    base_c2 = get_param(params, 'base_c2', sqrt(2) / (10 * numSubcarriers));

    profile.scheme = 'baseline';
    profile.definition.type = 'static';
    profile.definition.grouping = 'none';
    profile.c2 = base_c2;
    profile.group_index = ones(numSubcarriers, 1);
    profile.description = 'Scalar pre-chirp baseline.';
end
