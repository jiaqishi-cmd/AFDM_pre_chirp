function encoded = encode_mode_index(modeIndex, numModes, repetition)
%ENCODE_MODE_INDEX Encode an explicit waveform-mode index with repetition.

    if nargin < 3 || isempty(repetition)
        repetition = 1;
    end
    validateattributes(numModes, {'numeric'}, {'scalar', 'integer', '>=', 2});
    validateattributes(modeIndex, {'numeric'}, ...
        {'scalar', 'integer', '>=', 1, '<=', numModes});
    validateattributes(repetition, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});

    numBits = ceil(log2(numModes));
    value = modeIndex - 1;
    bits = bitget(value, numBits:-1:1).';
    codedBits = repelem(bits, repetition);

    encoded.mode_index = modeIndex;
    encoded.num_modes = numModes;
    encoded.num_information_bits = numBits;
    encoded.repetition = repetition;
    encoded.bits = bits;
    encoded.coded_bits = codedBits;
    encoded.symbols = 1 - 2 * codedBits;
end
