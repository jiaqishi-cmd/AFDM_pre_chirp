function profile = paper_grouping_profile(numSubcarriers, params)
%PAPER_GROUPING_PROFILE Build GPS candidate set and contiguous groups.
%   The final c2 vector is selected per frame by greedy_gps_profile.

    num_groups = get_param(params, 'num_groups', 4);
    num_candidates = get_param(params, 'num_candidates', 2);

    candidate_set = gps_candidate_set(numSubcarriers, num_candidates);
    group_index = contiguous_groups(numSubcarriers, num_groups);

    profile.scheme = 'paper_grouping';
    profile.c2 = candidate_set(:, 1);
    profile.group_index = group_index(:);
    profile.candidate_set = candidate_set;
    profile.description = 'GPS candidate set and contiguous groups for paper reproduction.';
end
