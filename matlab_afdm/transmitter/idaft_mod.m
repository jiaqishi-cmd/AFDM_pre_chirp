function signal = idaft_mod(symbols, numSubcarriers, c1, c2)
%IDAFT_MOD Compatibility wrapper.
%   Prefer afdm.tx.idaft_mod(symbols, numSubcarriers, c1, c2).

    if nargin < 4
        c2 = [];
    end
    if nargin < 3
        c1 = [];
    end
    signal = afdm.tx.idaft_mod(symbols, numSubcarriers, c1, c2);
end
