function decoded = decode_mode_index(receivedSymbols, numModes, repetition)
%DECODE_MODE_INDEX Decode a repeated BPSK waveform-mode index.

    if nargin < 3 || isempty(repetition)
        repetition = 1;
    end
    validateattributes(numModes, {'numeric'}, {'scalar', 'integer', '>=', 2});
    validateattributes(repetition, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});

    numBits = ceil(log2(numModes));
    receivedSymbols = receivedSymbols(:);
    if numel(receivedSymbols) ~= numBits * repetition
        error('receivedSymbols must contain ceil(log2(numModes))*repetition samples.');
    end

    repeated = reshape(real(receivedSymbols), repetition, numBits);
    softValues = sum(repeated, 1).';
    bits = softValues < 0;
    weights = 2.^(numBits-1:-1:0).';
    rawIndex = 1 + sum(double(bits) .* weights);
    valid = rawIndex <= numModes;

    decoded.raw_mode_index = rawIndex;
    decoded.mode_index = min(rawIndex, numModes);
    decoded.valid = valid;
    decoded.bits = bits;
    decoded.soft_values = softValues;
end
