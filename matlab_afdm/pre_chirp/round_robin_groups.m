function group_index = round_robin_groups(numSubcarriers, numGroups)
%ROUND_ROBIN_GROUPS Assign subcarriers to groups in an interleaved pattern.

    validate_group_count(numSubcarriers, numGroups);
    group_index = mod((0:numSubcarriers-1).', numGroups) + 1;
end
