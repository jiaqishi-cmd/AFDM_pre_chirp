% Main single-run AFDM simulation script.

rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
setup_paths(rootDir);

config = afdm_config();
frame = simulate_frame(config, config.simulation.random_seed);

fprintf('Transmit signal PAPR: %.2f dB\n', frame.papr);
fprintf('Total bits: %d\n', frame.total_bits);
fprintf('Error bits: %d\n', frame.err_bits);
fprintf('BER: %.2e\n', frame.ber);
fprintf('Simulation complete.\n');
