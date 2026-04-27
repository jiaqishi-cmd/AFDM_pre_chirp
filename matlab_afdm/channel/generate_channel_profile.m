function config = generate_channel_profile(config, profile)
%GENERATE_CHANNEL_PROFILE Populate channel taps for a named profile.
%   random_3path keeps the original lightweight demo channel.
%   bemani_21path follows the Fig. 5 setup in Bemani et al.:
%   delays 0:2, integer Doppler shifts -3:3 for each delay tap.

    if nargin < 2 || isempty(profile)
        profile = config.channel.profile;
    end

    N = config.waveform.NumSubcarriers;
    config.channel.profile = profile;

    switch lower(profile)
        case 'random_3path'
            taps = config.channel.num_taps;
            max_delay = config.channel.max_delay;
            max_doppler = config.channel.max_doppler;

            config.channel.delay_taps = randi([0, max_delay], 1, taps);
            config.channel.doppler_taps = (2 * rand(1, taps) - 1) * max_doppler;
            config.channel.doppler_freq = config.channel.doppler_taps / N;
            config.channel.chan_coef = (randn(1, taps) + 1i * randn(1, taps)) / sqrt(2);
            config.channel.chan_coef = config.channel.chan_coef / norm(config.channel.chan_coef);

        case 'bemani_21path'
            max_delay = config.channel.max_integer_delay;
            max_doppler = config.channel.max_integer_doppler;

            delays = 0:max_delay;
            dopplers = -max_doppler:max_doppler;
            [delay_grid, doppler_grid] = ndgrid(delays, dopplers);

            config.channel.delay_taps = delay_grid(:).';
            config.channel.doppler_taps = doppler_grid(:).';
            config.channel.doppler_freq = config.channel.doppler_taps / N;

            paths = numel(config.channel.delay_taps);
            config.channel.chan_coef = (randn(1, paths) + 1i * randn(1, paths)) / sqrt(2 * paths);

        otherwise
            error('Unsupported channel profile: %s', profile);
    end
end
