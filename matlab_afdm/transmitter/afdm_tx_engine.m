function [signal_cpp, papr, bits] = afdm_tx_engine(config)
%AFDM_TX_ENGINE Run the AFDM transmitter chain.
%   [signal_cpp, papr, bits] = afdm_tx_engine(config)

    if nargin < 1 || isempty(config)
        config = afdm_config();
    end

    scriptDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(scriptDir));

    M_mod = config.modulation.M_mod;
    modType = config.modulation.modType;
    c1 = config.waveform.c1;
    c2 = config.waveform.c2;
    N_subcarriers = config.waveform.NumSubcarriers;
    cpp_length = config.waveform.CPPLength;

    [symbols, bits] = random_data_generator(N_subcarriers, M_mod, modType);
    baseband = idaft_mod(symbols, N_subcarriers, c1, c2);
    papr = compute_papr(baseband);
    signal_cpp = add_cpp(baseband, cpp_length, c1);
end
