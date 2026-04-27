function [symbols, bits] = random_data_generator(numSymbols, M_mod, modType)
%RANDOM_DATA_GENERATOR 生成随机调制符号。
%   [symbols, bits] = random_data_generator(numSymbols, M_mod, modType)
%   numSymbols: 符号数量，默认 32
%   M_mod: 调制阶数（如 4 for QPSK, 16 for 16-QAM），默认 4
%   modType: 调制类型 ('qam' 或 'psk')，默认 'qam'
%
%   symbols: 调制后的复符号向量
%   bits: 随机生成的符号索引

    if nargin < 3 || isempty(modType)
        modType = 'qam';
    end
    if nargin < 2 || isempty(M_mod)
        M_mod = 4;
    end
    if nargin < 1 || isempty(numSymbols)
        numSymbols = 32;
    end

    k = log2(M_mod);
    total_bits = numSymbols * k;
    bits = randi([0, 1], total_bits, 1);
    
    switch lower(modType)
        case 'qam'
            symbols = qammod(bits, M_mod, 'InputType', 'bit', 'UnitAveragePower', true);
        case 'psk'
            symbols = pskmod(bits, M_mod, 'InputType', 'bit');
        otherwise
            error('Unsupported modulation type: %s', modType);
    end
end
