function group_index = contiguous_groups(numSubcarriers, numGroups)
%CONTIGUOUS_GROUPS Assign neighboring subcarriers to the same group.

    validate_group_count(numSubcarriers, numGroups);

    edges = round(linspace(0, numSubcarriers, numGroups + 1));
    group_index = zeros(numSubcarriers, 1);

    for group_id = 1:numGroups
        idx_start = edges(group_id) + 1;
        idx_end = edges(group_id + 1);
        group_index(idx_start:idx_end) = group_id;
    end
end
