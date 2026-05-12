function s = combine_partial_waveform(sPart, pattern)
%COMBINE_PARTIAL_WAVEFORM Combine precomputed partial waveforms.

    s = zeros(size(sPart{1, 1}));
    for groupId = 1:numel(pattern)
        s = s + sPart{groupId, pattern(groupId)};
    end
end
