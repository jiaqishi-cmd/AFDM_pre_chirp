function results = run_scheme_comparison(numPaprFrames, numBerFrames, snr_values)
%RUN_SCHEME_COMPARISON Compare baseline, GPS, and proposed pre-chirp schemes.
%   results = run_scheme_comparison(numPaprFrames, numBerFrames, snr_values)

    rootDir = fileparts(mfilename('fullpath'));
    addpath(rootDir);
    setup_paths(rootDir);

    if nargin < 1 || isempty(numPaprFrames)
        numPaprFrames = 100;
    end
    if nargin < 2 || isempty(numBerFrames)
        numBerFrames = 3;
    end
    if nargin < 3 || isempty(snr_values)
        snr_values = 0:5:30;
    end

    schemes = {'baseline', 'paper_grouping', 'proposed_grouping'};
    scheme_labels = {'Baseline', 'GPS', 'Proposed'};
    num_schemes = numel(schemes);

    base_config = afdm_config();
    base_seed = base_config.simulation.random_seed;
    frame_options.refresh_channel = base_config.simulation.refresh_channel_per_frame;

    results.schemes = schemes;
    results.scheme_labels = scheme_labels;
    results.snr_values = snr_values;
    results.num_papr_frames = numPaprFrames;
    results.num_ber_frames = numBerFrames;
    results.papr_samples = zeros(numPaprFrames, num_schemes);
    results.ber = zeros(numel(snr_values), num_schemes);

    fprintf('========== PAPR sampling ==========\n');
    for frame_idx = 1:numPaprFrames
        for scheme_idx = 1:num_schemes
            cfg = apply_pre_chirp_scheme(base_config, schemes{scheme_idx});
            frame = simulate_frame(cfg, base_seed + frame_idx, frame_options);
            results.papr_samples(frame_idx, scheme_idx) = frame.papr;
        end

        if mod(frame_idx, 25) == 0 || frame_idx == numPaprFrames
            fprintf('PAPR frame %4d/%d complete\n', frame_idx, numPaprFrames);
        end
    end

    fprintf('\n========== BER sweep ==========\n');
    for snr_idx = 1:numel(snr_values)
        snr_db = snr_values(snr_idx);

        for scheme_idx = 1:num_schemes
            total_err_bits = 0;
            total_bits = 0;

            cfg = apply_pre_chirp_scheme(base_config, schemes{scheme_idx});
            cfg.channel.snr_db = snr_db;

            for frame_idx = 1:numBerFrames
                frame = simulate_frame(cfg, base_seed + 100000 * snr_idx + frame_idx, frame_options);
                total_err_bits = total_err_bits + frame.err_bits;
                total_bits = total_bits + frame.total_bits;
            end

            results.ber(snr_idx, scheme_idx) = total_err_bits / total_bits;
            fprintf('%-9s | SNR = %5.1f dB | BER = %.2e\n', ...
                scheme_labels{scheme_idx}, snr_db, results.ber(snr_idx, scheme_idx));
        end
    end

    papr_thresholds = 0:0.25:12;
    results.papr_thresholds = papr_thresholds;
    results.papr_ccdf = compute_papr_ccdf(results.papr_samples, papr_thresholds);

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    results.output_dir = outputDir;
    results.timestamp = timestamp;

    plot_scheme_results(results);

    matPath = fullfile(outputDir, ['scheme_comparison_' timestamp '.mat']);
    save(matPath, 'results');
    fprintf('\nSaved results to %s\n', matPath);
end

function ccdf = compute_papr_ccdf(papr_samples, thresholds)
    num_schemes = size(papr_samples, 2);
    ccdf = zeros(numel(thresholds), num_schemes);

    for scheme_idx = 1:num_schemes
        for threshold_idx = 1:numel(thresholds)
            ccdf(threshold_idx, scheme_idx) = ...
                mean(papr_samples(:, scheme_idx) > thresholds(threshold_idx));
        end
    end
end

function plot_scheme_results(results)
    outputDir = results.output_dir;
    timestamp = results.timestamp;

    figure('Name', 'PAPR CCDF', 'Color', 'w');
    semilogy(results.papr_thresholds, results.papr_ccdf, 'LineWidth', 2);
    grid on;
    xlabel('PAPR threshold (dB)');
    ylabel('Pr(PAPR > threshold)');
    title('PAPR CCDF');
    legend(results.scheme_labels, 'Location', 'southwest');
    saveas(gcf, fullfile(outputDir, ['papr_ccdf_' timestamp '.png']));

    figure('Name', 'BER vs SNR', 'Color', 'w');
    semilogy(results.snr_values, results.ber, 'o-', 'LineWidth', 2);
    grid on;
    xlabel('SNR (dB)');
    ylabel('BER');
    title('BER vs SNR');
    legend(results.scheme_labels, 'Location', 'southwest');
    saveas(gcf, fullfile(outputDir, ['ber_vs_snr_' timestamp '.png']));
end
