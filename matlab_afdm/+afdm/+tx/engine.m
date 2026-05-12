function [signalCpp, papr, bits, txState] = engine(config)
%ENGINE Run the AFDM transmitter chain.

    if nargin < 1 || isempty(config)
        config = afdm_config();
    end

    MMod = config.modulation.M_mod;
    modType = config.modulation.modType;
    c1 = config.waveform.c1;
    numSubcarriers = config.waveform.NumSubcarriers;
    cppLength = config.waveform.CPPLength;

    [symbols, bits] = build_tx_symbols(config, numSubcarriers, MMod, modType);
    txState.pre_chirp_profile = select_pre_chirp_for_symbols(symbols, config);
    txState.c2 = txState.pre_chirp_profile.c2;
    txState.scheme = txState.pre_chirp_profile.scheme;
    txState.symbols = symbols;

    baseband = afdm.tx.idaft_mod(symbols, numSubcarriers, c1, txState.c2);
    papr = afdm.tx.compute_papr(baseband);
    signalCpp = afdm.tx.add_cpp(baseband, cppLength, c1);
end

function [symbols, bits] = build_tx_symbols(config, numSubcarriers, MMod, modType)
    if isfield(config, 'tx') && isfield(config.tx, 'symbols') && ~isempty(config.tx.symbols)
        symbols = config.tx.symbols(:);
        if numel(symbols) ~= numSubcarriers
            error('config.tx.symbols length must equal the number of subcarriers.');
        end
        if ~isfield(config.tx, 'bits') || isempty(config.tx.bits)
            error('config.tx.bits must be provided when config.tx.symbols is provided.');
        end
        bits = validate_tx_bits(config.tx.bits, numSubcarriers, MMod);
        return;
    end

    if isfield(config, 'tx') && isfield(config.tx, 'bits') && ~isempty(config.tx.bits)
        [symbols, bits] = afdm.tx.random_data(numSubcarriers, MMod, modType, config.tx.bits);
    else
        [symbols, bits] = afdm.tx.random_data(numSubcarriers, MMod, modType);
    end
end

function bits = validate_tx_bits(bits, numSymbols, MMod)
    k = log2(MMod);
    expectedBits = numSymbols * k;
    bits = bits(:);

    if numel(bits) ~= expectedBits
        error('config.tx.bits length must equal numSymbols * log2(M_mod).');
    end
    if any(bits ~= 0 & bits ~= 1)
        error('config.tx.bits must contain only 0 and 1.');
    end
end
