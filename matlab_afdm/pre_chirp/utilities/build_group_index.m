function group_index = build_group_index(numSubcarriers, numGroups, grouping)
%BUILD_GROUP_INDEX Assign subcarriers to pre-chirp groups.

    if nargin < 3 || isempty(grouping)
        grouping = 'contiguous';
    end

    validateattributes(numSubcarriers, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(numGroups, {'numeric'}, {'scalar', 'integer', 'positive'});

    if numGroups > numSubcarriers
        error('num_groups must not exceed the number of subcarriers.');
    end

    switch lower(grouping)
        case 'contiguous'
            edges = round(linspace(0, numSubcarriers, numGroups + 1));
            group_index = zeros(numSubcarriers, 1);

            for group_id = 1:numGroups
                idx_start = edges(group_id) + 1;
                idx_end = edges(group_id + 1);
                group_index(idx_start:idx_end) = group_id;
            end
        case 'round_robin'
            group_index = mod((0:numSubcarriers-1).', numGroups) + 1;
        otherwise
            error('Unsupported pre-chirp grouping: %s', grouping);
    end
end
