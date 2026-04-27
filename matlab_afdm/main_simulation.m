% Main single-run AFDM simulation script.

rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
setup_paths(rootDir);

config = afdm_config();

[signal_cpp, papr, tx_bits, tx_state] = afdm_tx_engine(config);
fprintf('Transmit signal PAPR: %.2f dB\n', papr);

r_signal = multipath_channel(signal_cpp, config);
r_signal = add_awgn(r_signal, config);

[~, err_bits, total_bits] = afdm_rx_engine(r_signal, config, tx_bits, tx_state);

ber = err_bits / total_bits;
fprintf('Total bits: %d\n', total_bits);
fprintf('Error bits: %d\n', err_bits);
fprintf('BER: %.2e\n', ber);
fprintf('Simulation complete.\n');
