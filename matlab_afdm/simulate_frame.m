function frame = simulate_frame(config, frameSeed, options)
%SIMULATE_FRAME Run one AFDM transmit-channel-receive frame.
%   frame = simulate_frame(config, frameSeed, options)

    if nargin < 1 || isempty(config)
        config = afdm_config();
    end
    if nargin < 2
        frameSeed = [];
    end
    if nargin < 3
        options = struct();
    end

    oldRngState = [];
    if ~isempty(frameSeed)
        oldRngState = rng;
        rng(frameSeed, 'twister');
    end

    if get_option(options, 'refresh_channel', false)
        config = generate_channel_profile(config, config.channel.profile);
    end

    [signal_cpp, papr, tx_bits, tx_state] = afdm_tx_engine(config);
    r_signal = multipath_channel(signal_cpp, config);
    r_signal = add_awgn(r_signal, config);
    [x_dec, err_bits, total_bits] = afdm_rx_engine(r_signal, config, tx_bits, tx_state);

    frame.papr = papr;
    frame.err_bits = err_bits;
    frame.total_bits = total_bits;
    frame.ber = err_bits / total_bits;
    frame.x_dec = x_dec;
    frame.tx_bits = tx_bits;
    frame.tx_state = tx_state;
    frame.config = config;

    if ~isempty(oldRngState)
        rng(oldRngState);
    end
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name)
        value = options.(name);
    else
        value = defaultValue;
    end
end
