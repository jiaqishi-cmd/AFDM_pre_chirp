% Main single-run AFDM simulation script.

addpath(genpath('channel'));
addpath(genpath('receive'));
addpath(genpath('transmitter'));

% 1. Load configuration.
config = afdm_config();

% 2. Transmitter.
[signal_cpp, papr, tx_bits] = afdm_tx_engine(config);
fprintf('Transmit signal PAPR: %.2f dB\n', papr);

% 3. Channel.
r_signal = multipath_channel(signal_cpp, config);
r_signal = add_awgn(r_signal, config);

% 4. Receiver.
[x_dec, err_bits, total_bits] = afdm_rx_engine(r_signal, config, tx_bits);

% 5. Results.
ber = err_bits / total_bits;
fprintf('Total bits: %d\n', total_bits);
fprintf('Error bits: %d\n', err_bits);
fprintf('BER: %.2e\n', ber);
fprintf('Simulation complete.\n');
