function profile = gps_profile(numSubcarriers, params)
%GPS_PROFILE Build the GPS candidate set and grouping definition.
%   The final c2 vector is selected per frame by select_greedy_profile.

    numGroups = afdm.chirp.get_param(params, 'num_groups', 4);
    numCandidates = afdm.chirp.get_param(params, 'num_candidates', 2);
    grouping = afdm.chirp.get_param(params, 'grouping', 'contiguous');

    candidateSet = afdm.chirp.gps_candidate_set(numSubcarriers, numCandidates);
    groupIndex = afdm.chirp.group_index(numSubcarriers, numGroups, grouping);

    profile.scheme = 'paper_grouping';
    profile.definition.type = 'candidate_set';
    profile.definition.grouping = grouping;
    profile.c2 = candidateSet(:, 1);
    profile.group_index = groupIndex(:);
    profile.candidate_set = candidateSet;
    profile.description = 'GPS candidate set and grouping for paper reproduction.';
end
