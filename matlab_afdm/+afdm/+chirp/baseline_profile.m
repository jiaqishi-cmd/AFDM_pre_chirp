function profile = baseline_profile(numSubcarriers, params)
%BASELINE_PROFILE Use one scalar pre-chirp coefficient for all subcarriers.

    baseC2 = afdm.chirp.get_param(params, 'base_c2', sqrt(2) / (10 * numSubcarriers));

    profile.scheme = 'baseline';
    profile.definition.type = 'static';
    profile.definition.grouping = 'none';
    profile.c2 = baseC2;
    profile.group_index = ones(numSubcarriers, 1);
    profile.description = 'Scalar pre-chirp baseline.';
end
