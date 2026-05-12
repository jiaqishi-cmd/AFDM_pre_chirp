function [symbols, bits] = random_data(numSymbols, MMod, modType, bits)
%RANDOM_DATA Generate random modulated symbols or map supplied bits.

    if nargin < 3 || isempty(modType)
        modType = 'qam';
    end
    if nargin < 2 || isempty(MMod)
        MMod = 4;
    end
    if nargin < 1 || isempty(numSymbols)
        numSymbols = 32;
    end

    validateattributes(numSymbols, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(MMod, {'numeric'}, {'scalar', 'integer', '>', 1});

    k = log2(MMod);
    if k ~= floor(k)
        error('M_mod must be a power of two.');
    end

    totalBits = numSymbols * k;
    if nargin < 4 || isempty(bits)
        bits = randi([0, 1], totalBits, 1);
    else
        bits = bits(:);
        if numel(bits) ~= totalBits
            error('Input bits length must equal numSymbols * log2(M_mod).');
        end
        if any(bits ~= 0 & bits ~= 1)
            error('Input bits must contain only 0 and 1.');
        end
    end

    switch lower(modType)
        case 'qam'
            symbols = qammod(bits, MMod, 'InputType', 'bit', 'UnitAveragePower', true);
        case 'psk'
            symbols = pskmod(bits, MMod, pi/MMod, 'InputType', 'bit');
        otherwise
            error('Unsupported modulation type: %s', modType);
    end
end
