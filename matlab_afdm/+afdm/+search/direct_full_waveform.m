function s = direct_full_waveform(symbols, baseC2, offsets, groupIndex, pattern, os)
%DIRECT_FULL_WAVEFORM Build waveform from a complete group pattern.

    M = numel(symbols);
    c2Vec = baseC2 * ones(M, 1);
    for groupId = 1:numel(pattern)
        c2Vec(groupIndex == groupId) = baseC2 + offsets(pattern(groupId));
    end
    s = afdm.search.full_waveform(symbols, c2Vec, os);
end
