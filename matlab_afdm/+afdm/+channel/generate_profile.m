function config = generate_profile(config, profile)
%GENERATE_PROFILE Populate channel taps for a named profile.

    if nargin < 2 || isempty(profile)
        profile = config.channel.profile;
    end

    N = config.waveform.NumSubcarriers;
    config.channel.profile = profile;

    switch lower(profile)
        case 'random_3path'
            taps = config.channel.num_taps;
            maxDelay = config.channel.max_delay;
            maxDoppler = config.channel.max_doppler;

            config.channel.delay_taps = randi([0, maxDelay], 1, taps);
            config.channel.doppler_taps = (2 * rand(1, taps) - 1) * maxDoppler;
            config.channel.doppler_freq = config.channel.doppler_taps / N;
            config.channel.chan_coef = (randn(1, taps) + 1i * randn(1, taps)) / sqrt(2);
            config.channel.chan_coef = config.channel.chan_coef / norm(config.channel.chan_coef);

        case 'bemani_21path'
            maxDelay = config.channel.max_integer_delay;
            maxDoppler = config.channel.max_integer_doppler;

            delays = 0:maxDelay;
            dopplers = -maxDoppler:maxDoppler;
            [delayGrid, dopplerGrid] = ndgrid(delays, dopplers);

            config.channel.delay_taps = delayGrid(:).';
            config.channel.doppler_taps = dopplerGrid(:).';
            config.channel.doppler_freq = config.channel.doppler_taps / N;

            paths = numel(config.channel.delay_taps);
            config.channel.chan_coef = (randn(1, paths) + 1i * randn(1, paths)) / sqrt(2 * paths);

        otherwise
            error('Unsupported channel profile: %s', profile);
    end
end
