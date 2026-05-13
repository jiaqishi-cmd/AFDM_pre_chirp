function [xDec, errBits, totalBits] = engine(rSignal, config, xRef, rxState)
%ENGINE Run the AFDM receiver chain and compute bit errors.

    if nargin >= 4 && isstruct(rxState) && isfield(rxState, 'c2')
        c2 = rxState.c2;
    else
        c2 = config.waveform.c2;
    end

    rData = afdm.rx.remove_cpp(rSignal, config.waveform.CPPLength);

    yDaft = afdm.rx.daft_demod( ...
        rData, ...
        config.waveform.NumSubcarriers, ...
        config.waveform.c1, ...
        c2);

    if isfield(config.channel, 'multipath') && ~config.channel.multipath
        hEff = eye(config.waveform.NumSubcarriers);
    else
        hEff = afdm.rx.estimate_effective_channel( ...
            config.waveform.NumSubcarriers, ...
            config.waveform.c1, ...
            c2, ...
            config.channel.chan_coef, ...
            config.channel.delay_taps, ...
            config.channel.doppler_freq);
    end

    noiseVar = resolve_noise_var(config);
    xEst = afdm.rx.equalize_symbols(yDaft, hEff, noiseVar, config);
    xDec = afdm.rx.symbol_decision(xEst, config.modulation.M_mod, config.modulation.modType);
    [errBits, totalBits] = afdm.rx.compute_bit_errors(xDec, xRef);
end

function noiseVar = resolve_noise_var(config)
    if ~isfield(config.channel, 'snr_db') || isempty(config.channel.snr_db)
        error('config.channel.snr_db must be defined.');
    end

    if isfield(config.channel, 'noise_var') && ~isempty(config.channel.noise_var)
        noiseVar = config.channel.noise_var;
    elseif isfield(config.channel, 'add_noise') && ~config.channel.add_noise
        noiseVar = 0;
    else
        noiseVar = 1 / (10^(config.channel.snr_db / 10));
    end
end
