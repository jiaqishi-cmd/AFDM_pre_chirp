function profile = gps_profile(numSubcarriers, params)
%GPS_PROFILE Build the GPS candidate set and grouping definition.
%   The final c2 vector is selected per frame by select_greedy_profile.

    num_groups = get_param(params, 'num_groups', 4);
    num_candidates = get_param(params, 'num_candidates', 2);
    grouping = get_param(params, 'grouping', 'contiguous');

    candidate_set = gps_candidate_set(numSubcarriers, num_candidates);
    group_index = build_group_index(numSubcarriers, num_groups, grouping);

    profile.scheme = 'paper_grouping';
    profile.definition.type = 'candidate_set';
    profile.definition.grouping = grouping;
    profile.c2 = candidate_set(:, 1);
    profile.group_index = group_index(:);
    profile.candidate_set = candidate_set;
    profile.description = 'GPS candidate set and grouping for paper reproduction.';
end
