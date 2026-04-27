function validate_group_count(numSubcarriers, numGroups)
%VALIDATE_GROUP_COUNT Validate pre-chirp grouping dimensions.

    validateattributes(numSubcarriers, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(numGroups, {'numeric'}, {'scalar', 'integer', 'positive'});

    if numGroups > numSubcarriers
        error('num_groups must not exceed the number of subcarriers.');
    end
end
