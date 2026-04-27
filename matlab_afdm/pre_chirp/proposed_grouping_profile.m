function profile = proposed_grouping_profile(numSubcarriers, params)
%PROPOSED_GROUPING_PROFILE Build small-perturbation c2 candidates.
%   c2_m = c2^(0) + Delta^(w_v), Delta in {-delta, 0, +delta}, m in G_v.

    base_c2 = get_param(params, 'base_c2', sqrt(2) / (10 * numSubcarriers));
    num_groups = get_param(params, 'num_groups', 4);
    delta = get_param(params, 'delta', base_c2 / 16);

    group_index = contiguous_groups(numSubcarriers, num_groups);
    candidate_offsets = [0, -delta, delta];
    candidate_set = base_c2 + repmat(candidate_offsets, numSubcarriers, 1);

    profile.scheme = 'proposed_grouping';
    profile.c2 = candidate_set(:, 2);
    profile.group_index = group_index(:);
    profile.candidate_set = candidate_set;
    profile.delta = delta;
    profile.description = 'Small-perturbation grouped pre-chirp profile.';
end
