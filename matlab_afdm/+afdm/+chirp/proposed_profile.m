function profile = proposed_profile(numSubcarriers, params)
%PROPOSED_PROFILE Build small-perturbation c2 candidates.
%   c2_m = c2^(0) + Delta^(w_v), Delta in {0, -delta, +delta}, m in G_v.

    baseC2 = afdm.chirp.get_param(params, 'base_c2', sqrt(2) / (10 * numSubcarriers));
    numGroups = afdm.chirp.get_param(params, 'num_groups', 4);
    delta = afdm.chirp.get_param(params, 'delta', baseC2 / 16);
    grouping = afdm.chirp.get_param(params, 'grouping', 'contiguous');

    groupIndex = afdm.chirp.group_index(numSubcarriers, numGroups, grouping);
    candidateOffsets = [0, -delta, delta];
    candidateSet = baseC2 + repmat(candidateOffsets, numSubcarriers, 1);

    profile.scheme = 'proposed_grouping';
    profile.definition.type = 'candidate_set';
    profile.definition.grouping = grouping;
    profile.c2 = candidateSet(:, 1);
    profile.group_index = groupIndex(:);
    profile.candidate_set = candidateSet;
    profile.delta = delta;
    profile.description = 'Small-perturbation grouped pre-chirp profile.';
end
