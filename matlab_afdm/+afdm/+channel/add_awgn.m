function [y, noiseVar] = add_awgn(signal, config)
%ADD_AWGN Add complex AWGN using config.channel.snr_db.

    if ~isfield(config, 'channel') || ~isfield(config.channel, 'snr_db')
        error('config.channel.snr_db not found.');
    end

    signal = complex(signal);
    snrDb = config.channel.snr_db;

    if ~config.channel.add_noise
        y = signal;
        noiseVar = 0;
        return;
    end

    power = mean(abs(signal).^2);
    snr = 10^(snrDb / 10);
    noiseVar = power / snr;

    noise = sqrt(noiseVar / 2) * (randn(size(signal)) + 1i * randn(size(signal)));
    y = signal + noise;
end
