function results = run_ber_comparison(numFrames, snr_values, schemes)
%RUN_BER_COMPARISON Run BER curves for selected pre-chirp schemes.
%   results = run_ber_comparison(numFrames, snr_values, schemes)

    rootDir = fileparts(mfilename('fullpath'));
    addpath(rootDir);
    setup_paths(rootDir);

    if nargin < 1 || isempty(numFrames)
        numFrames = 200;
    end
    if nargin < 2 || isempty(snr_values)
        snr_values = 15:2.5:30;
    end
    if nargin < 3 || isempty(schemes)
        schemes = {'baseline', 'paper_grouping', 'proposed_grouping'};
    end

    scheme_labels = scheme_display_labels(schemes);
    num_schemes = numel(schemes);
    base_config = afdm_config();
    base_seed = base_config.simulation.random_seed;

    results.schemes = schemes;
    results.scheme_labels = scheme_labels;
    results.snr_values = snr_values;
    results.num_frames = numFrames;
    results.ber = zeros(numel(snr_values), num_schemes);
    results.error_bits = zeros(numel(snr_values), num_schemes);
    results.total_bits = zeros(numel(snr_values), num_schemes);

    fprintf('========== BER comparison ==========\n');
    fprintf('Frames per SNR per scheme: %d\n\n', numFrames);

    for snr_idx = 1:numel(snr_values)
        snr_db = snr_values(snr_idx);

        for scheme_idx = 1:num_schemes
            cfg = apply_pre_chirp_scheme(base_config, schemes{scheme_idx});
            cfg.channel.snr_db = snr_db;

            total_err_bits = 0;
            total_bits = 0;

            for frame_idx = 1:numFrames
                rng(base_seed + 1000000 * snr_idx + 1000 * scheme_idx + frame_idx, 'twister');

                [signal_cpp, ~, tx_bits, tx_state] = afdm_tx_engine(cfg);
                r_signal = multipath_channel(signal_cpp, cfg);
                r_signal = add_awgn(r_signal, cfg);
                [~, err_bits, num_bits] = afdm_rx_engine(r_signal, cfg, tx_bits, tx_state);

                total_err_bits = total_err_bits + err_bits;
                total_bits = total_bits + num_bits;
            end

            results.error_bits(snr_idx, scheme_idx) = total_err_bits;
            results.total_bits(snr_idx, scheme_idx) = total_bits;
            results.ber(snr_idx, scheme_idx) = total_err_bits / total_bits;

            fprintf('%-9s | SNR = %5.1f dB | errors = %5d / %-7d | BER = %.3e\n', ...
                scheme_labels{scheme_idx}, snr_db, total_err_bits, total_bits, ...
                results.ber(snr_idx, scheme_idx));
        end
    end

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    results.output_dir = outputDir;
    results.timestamp = timestamp;

    figure('Name', 'BER comparison', 'Color', 'w');
    semilogy(results.snr_values, results.ber, 'o-', 'LineWidth', 2);
    grid on;
    xlabel('SNR (dB)');
    ylabel('BER');
    title(sprintf('BER comparison (%d frames per point)', numFrames));
    legend(results.scheme_labels, 'Location', 'southwest');
    saveas(gcf, fullfile(outputDir, ['ber_comparison_' timestamp '.png']));

    matPath = fullfile(outputDir, ['ber_comparison_' timestamp '.mat']);
    save(matPath, 'results');
    fprintf('\nSaved results to %s\n', matPath);
end

function labels = scheme_display_labels(schemes)
    labels = schemes;

    for idx = 1:numel(schemes)
        switch lower(schemes{idx})
            case 'baseline'
                labels{idx} = 'Baseline';
            case 'paper_grouping'
                labels{idx} = 'GPS';
            case 'proposed_grouping'
                labels{idx} = 'Proposed';
            otherwise
                labels{idx} = schemes{idx};
        end
    end
end
