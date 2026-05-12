function sPart = precompute_partial_waveforms(symbols, baseC2, offsets, groupIndex, os)
%PRECOMPUTE_PARTIAL_WAVEFORMS Precompute group/candidate partial waveforms.

    V = max(groupIndex);
    W = numel(offsets);
    M = numel(symbols);
    sPart = cell(V, W);

    for groupId = 1:V
        groupMask = groupIndex == groupId;
        for candId = 1:W
            c2Vec = baseC2 * ones(M, 1);
            c2Vec(groupMask) = baseC2 + offsets(candId);
            sPart{groupId, candId} = partial_waveform_for_group(symbols, c2Vec, groupMask, os);
        end
    end
end

function s = partial_waveform_for_group(symbols, c2Vec, groupMask, os)
    M = numel(symbols);
    m = (0:M-1).';
    xPre = zeros(M, 1);
    xPre(groupMask) = symbols(groupMask) .* exp(1i * 2 * pi .* c2Vec(groupMask) .* (m(groupMask).^2));
    s = afdm.search.ifft_oversampled(xPre, os);
end
