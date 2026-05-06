function [signal_cpp, papr, bits, tx_state] = afdm_tx_engine(config)
%AFDM_TX_ENGINE Run the AFDM transmitter chain.
%   [signal_cpp, papr, bits, tx_state] = afdm_tx_engine(config)

    if nargin < 1 || isempty(config)
        config = afdm_config();
    end

    M_mod = config.modulation.M_mod;
    modType = config.modulation.modType;
    c1 = config.waveform.c1;
    N_subcarriers = config.waveform.NumSubcarriers;
    cpp_length = config.waveform.CPPLength;

    [symbols, bits] = build_tx_symbols(config, N_subcarriers, M_mod, modType);
    tx_state.pre_chirp_profile = select_pre_chirp_for_symbols(symbols, config);
    tx_state.c2 = tx_state.pre_chirp_profile.c2;
    tx_state.scheme = tx_state.pre_chirp_profile.scheme;
    tx_state.symbols = symbols;

    c2 = tx_state.c2;
    baseband = idaft_mod(symbols, N_subcarriers, c1, c2);
    papr = compute_papr(baseband);
    signal_cpp = add_cpp(baseband, cpp_length, c1);
end

function [symbols, bits] = build_tx_symbols(config, numSubcarriers, M_mod, modType)
    if isfield(config, 'tx') && isfield(config.tx, 'symbols') && ~isempty(config.tx.symbols)
        symbols = config.tx.symbols(:);
        if numel(symbols) ~= numSubcarriers
            error('config.tx.symbols length must equal the number of subcarriers.');
        end
        if ~isfield(config.tx, 'bits') || isempty(config.tx.bits)
            error('config.tx.bits must be provided when config.tx.symbols is provided.');
        end
        bits = validate_tx_bits(config.tx.bits, numSubcarriers, M_mod);
        return;
    end

    if isfield(config, 'tx') && isfield(config.tx, 'bits') && ~isempty(config.tx.bits)
        [symbols, bits] = random_data_generator(numSubcarriers, M_mod, modType, config.tx.bits);
    else
        [symbols, bits] = random_data_generator(numSubcarriers, M_mod, modType);
    end
end

function bits = validate_tx_bits(bits, numSymbols, M_mod)
    k = log2(M_mod);
    expectedBits = numSymbols * k;
    bits = bits(:);

    if numel(bits) ~= expectedBits
        error('config.tx.bits length must equal numSymbols * log2(M_mod).');
    end
    if any(bits ~= 0 & bits ~= 1)
        error('config.tx.bits must contain only 0 and 1.');
    end
end
