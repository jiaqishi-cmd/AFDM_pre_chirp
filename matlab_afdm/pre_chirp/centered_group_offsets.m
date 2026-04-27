function offsets = centered_group_offsets(numGroups, spacing)
%CENTERED_GROUP_OFFSETS Return zero-centered coefficient offsets.

    validateattributes(numGroups, {'numeric'}, {'scalar', 'integer', 'positive'});
    center = (numGroups + 1) / 2;
    offsets = ((1:numGroups).' - center) * spacing;
end
