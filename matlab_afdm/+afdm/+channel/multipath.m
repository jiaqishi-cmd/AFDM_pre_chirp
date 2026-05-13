function y = multipath(signal, config)
%MULTIPATH Apply the configured time-domain multipath Doppler channel.

    if ~isfield(config, 'channel') || ~isfield(config.channel, 'multipath')
        error('config.channel.multipath not found.');
    end

    if ~config.channel.multipath
        y = signal;
        return;
    end

    delayTaps = config.channel.delay_taps;
    dopplerFreq = config.channel.doppler_freq;
    chanCoef = config.channel.chan_coef;
    cppLength = config.waveform.CPPLength;

    taps = numel(delayTaps);
    if numel(dopplerFreq) ~= taps || numel(chanCoef) ~= taps
        error('delay_taps, doppler_freq, and chan_coef must have the same length.');
    end

    signal = signal(:);
    nTotal = numel(signal);
    y = zeros(nTotal, 1);

    for nAbs = 1:nTotal
        t = nAbs - (cppLength + 1);
        for pathIdx = 1:taps
            delay = delayTaps(pathIdx);
            doppler = dopplerFreq(pathIdx);
            if (nAbs - delay) > 0
                sample = signal(nAbs - delay);
                hVal = chanCoef(pathIdx) * exp(-1i * 2 * pi * doppler * t);
                y(nAbs) = y(nAbs) + hVal * sample;
            end
        end
    end
end
