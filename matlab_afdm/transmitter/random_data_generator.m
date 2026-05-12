function [symbols, bits] = random_data_generator(numSymbols, M_mod, modType, bits)
%RANDOM_DATA_GENERATOR Compatibility wrapper.
%   Prefer afdm.tx.random_data(numSymbols, M_mod, modType, bits).

    if nargin < 4
        bits = [];
    end
    if nargin < 3
        modType = [];
    end
    if nargin < 2
        M_mod = [];
    end
    if nargin < 1
        numSymbols = [];
    end
    [symbols, bits] = afdm.tx.random_data(numSymbols, M_mod, modType, bits);
end
