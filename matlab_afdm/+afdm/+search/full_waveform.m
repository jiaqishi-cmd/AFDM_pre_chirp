function s = full_waveform(symbols, c2Vec, os)
%FULL_WAVEFORM Build the full pre-chirped waveform for PAPR search.

    M = numel(symbols);
    m = (0:M-1).';
    xPre = symbols(:) .* exp(1i * 2 * pi .* c2Vec(:) .* (m.^2));
    s = afdm.search.ifft_oversampled(xPre, os);
end
