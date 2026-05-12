function groupIndex = group_index(numSubcarriers, numGroups, grouping)
%GROUP_INDEX Assign subcarriers to pre-chirp groups.

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
            groupIndex = zeros(numSubcarriers, 1);
            for groupId = 1:numGroups
                idxStart = edges(groupId) + 1;
                idxEnd = edges(groupId + 1);
                groupIndex(idxStart:idxEnd) = groupId;
            end
        case 'round_robin'
            groupIndex = mod((0:numSubcarriers-1).', numGroups) + 1;
        otherwise
            error('Unsupported pre-chirp grouping: %s', grouping);
    end
end
