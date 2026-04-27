function config = afdm_config()
%AFDM_CONFIG Build the default AFDM simulation configuration.

    config.simulation.random_seed = 20260427;
    rng(config.simulation.random_seed, 'twister');

    config.waveform.NumSubcarriers = 64;
    config.waveform.CPPLength = 3;

    config.modulation.M_mod = 16;
    config.modulation.modType = 'qam';

    config.channel.snr_db = 30;
    config.channel.add_noise = true;
    config.channel.multipath = true;

    taps = 3;
    max_delay = 3;
    k_max = 2.5;
    config.channel.delay_taps = randi([0, max_delay], 1, taps);

    doppler_taps = (2 * rand(1, taps) - 1) * k_max;
    config.channel.doppler_freq = doppler_taps / config.waveform.NumSubcarriers;

    config.channel.chan_coef = (randn(1, taps) + 1i * randn(1, taps)) / sqrt(2);
    config.channel.chan_coef = config.channel.chan_coef / norm(config.channel.chan_coef);

    k_v = 1;
    config.waveform.c1 = ...
        (2 * (floor(max(abs(doppler_taps))) + k_v) + 1) ...
        / (2 * config.waveform.NumSubcarriers);

    base_c2 = sqrt(2) / (10 * config.waveform.NumSubcarriers);
    config.pre_chirp.scheme = 'baseline';
    config.pre_chirp.base_c2 = base_c2;
    config.pre_chirp.num_groups = 4;
    config.pre_chirp.num_candidates = 2;
    config.pre_chirp.group_spacing = base_c2 / 4;
    config.pre_chirp.delta = base_c2 / 16;
    config = apply_pre_chirp_scheme(config, config.pre_chirp.scheme);

    config.simulation.enable_snr_loop = true;
    config.simulation.snr_range = 0:5:30;
    config.simulation.num_trials = 1;
    config.simulation.enable_montecarlo = false;
end
