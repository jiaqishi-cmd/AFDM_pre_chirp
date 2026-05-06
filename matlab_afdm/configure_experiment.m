function config = configure_experiment(config, options)
%CONFIGURE_EXPERIMENT Apply optional experiment overrides to a config.

    if nargin < 2
        options = struct();
    end

    if isfield(options, 'M_mod')
        config.modulation.M_mod = options.M_mod;
    end
    if isfield(options, 'modType')
        config.modulation.modType = options.modType;
    end
    if isfield(options, 'channel_profile')
        config = generate_channel_profile(config, options.channel_profile);
        config.waveform.c1 = ...
            (2 * (floor(max(abs(config.channel.doppler_taps))) + 1) + 1) ...
            / (2 * config.waveform.NumSubcarriers);
    end
    if isfield(options, 'refresh_channel_per_frame')
        config.simulation.refresh_channel_per_frame = options.refresh_channel_per_frame;
    end

    config = apply_pre_chirp_scheme(config, config.pre_chirp.scheme);
end
