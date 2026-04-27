% AFDM simulation loop.
% Supports SNR sweep, Monte Carlo evaluation, and single-run modes.

rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
setup_paths(rootDir);

config = afdm_config();

if config.simulation.enable_snr_loop
    snr_values = config.simulation.snr_range;
    num_trials = config.simulation.num_trials;

    fprintf('========== SNR sweep simulation ==========\n');
    fprintf('SNR range: %s dB\n', mat2str(snr_values));
    fprintf('Trials per SNR point: %d\n\n', num_trials);

    ber_results = zeros(length(snr_values), 1);
    papr_results = zeros(length(snr_values), 1);

    for snr_idx = 1:length(snr_values)
        current_snr = snr_values(snr_idx);
        config.channel.snr_db = current_snr;

        total_err_bits = 0;
        total_bits = 0;
        total_papr = 0;

        for trial = 1:num_trials
            [signal_cpp, papr, tx_bits, tx_state] = afdm_tx_engine(config);
            total_papr = total_papr + papr;

            r_signal = multipath_channel(signal_cpp, config);
            r_signal = add_awgn(r_signal, config);

            [~, err_bits, num_bits] = afdm_rx_engine(r_signal, config, tx_bits, tx_state);

            total_err_bits = total_err_bits + err_bits;
            total_bits = total_bits + num_bits;
        end

        ber_results(snr_idx) = total_err_bits / total_bits;
        papr_results(snr_idx) = total_papr / num_trials;

        fprintf('SNR = %6.1f dB | BER = %.2e | Avg PAPR = %.2f dB\n', ...
            current_snr, ber_results(snr_idx), papr_results(snr_idx));
    end

    figure('Name', 'SNR sweep results');

    subplot(1, 2, 1);
    semilogy(snr_values, ber_results, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
    grid on;
    xlabel('SNR (dB)');
    ylabel('BER');
    title('AFDM BER performance');

    subplot(1, 2, 2);
    plot(snr_values, papr_results, 'rs-', 'LineWidth', 2, 'MarkerSize', 8);
    grid on;
    xlabel('SNR (dB)');
    ylabel('PAPR (dB)');
    title('PAPR vs SNR');

    fprintf('\n========== Simulation complete ==========\n');

elseif config.simulation.enable_montecarlo
    num_trials = config.simulation.num_trials;

    fprintf('========== Monte Carlo simulation (SNR = %.1f dB) ==========\n', ...
        config.channel.snr_db);
    fprintf('Trials: %d\n\n', num_trials);

    total_err_bits = 0;
    total_bits = 0;
    total_papr = 0;

    for trial = 1:num_trials
        [signal_cpp, papr, tx_bits, tx_state] = afdm_tx_engine(config);
        total_papr = total_papr + papr;

        r_signal = multipath_channel(signal_cpp, config);
        r_signal = add_awgn(r_signal, config);

        [~, err_bits, num_bits] = afdm_rx_engine(r_signal, config, tx_bits, tx_state);

        total_err_bits = total_err_bits + err_bits;
        total_bits = total_bits + num_bits;

        if mod(trial, 10) == 0 || trial == 1
            fprintf('Trial %4d/%d complete\n', trial, num_trials);
        end
    end

    ber = total_err_bits / total_bits;
    avg_papr = total_papr / num_trials;

    fprintf('\n========== Results ==========\n');
    fprintf('Total bits: %d\n', total_bits);
    fprintf('Error bits: %d\n', total_err_bits);
    fprintf('BER: %.2e\n', ber);
    fprintf('Avg PAPR: %.2f dB\n', avg_papr);
    fprintf('========== Simulation complete ==========\n');

else
    fprintf('========== Single-run simulation ==========\n');

    [signal_cpp, papr, tx_bits, tx_state] = afdm_tx_engine(config);
    fprintf('Transmit signal PAPR: %.2f dB\n', papr);

    r_signal = multipath_channel(signal_cpp, config);
    r_signal = add_awgn(r_signal, config);

    [~, err_bits, total_bits] = afdm_rx_engine(r_signal, config, tx_bits, tx_state);

    ber = err_bits / total_bits;
    fprintf('Total bits: %d\n', total_bits);
    fprintf('Error bits: %d\n', err_bits);
    fprintf('BER: %.2e\n', ber);
    fprintf('========== Simulation complete ==========\n');
end
