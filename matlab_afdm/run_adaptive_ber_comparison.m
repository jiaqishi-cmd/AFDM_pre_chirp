function results = run_adaptive_ber_comparison(snr_values, targetErrors, maxFrames, schemes)
%RUN_ADAPTIVE_BER_COMPARISON Run BER curves with target-error stopping.
%   results = run_adaptive_ber_comparison(snr_values, targetErrors, maxFrames, schemes)
%
%   Each SNR/scheme point stops when either targetErrors bit errors are
%   collected or maxFrames frames are simulated. This spends little time at
%   low SNR and reserves long runs for high-SNR tail points.

    rootDir = fileparts(mfilename('fullpath'));
    addpath(rootDir);
    setup_paths(rootDir);

    if nargin < 1 || isempty(snr_values)
        snr_values = [15, 17.5, 20, 22.5, 25, 27.5];
    end
    if nargin < 2 || isempty(targetErrors)
        targetErrors = default_target_errors(snr_values);
    end
    if nargin < 3 || isempty(maxFrames)
        maxFrames = default_max_frames(snr_values);
    end
    if nargin < 4 || isempty(schemes)
        schemes = {'baseline', 'paper_grouping', 'proposed_grouping'};
    end

    targetErrors = expand_point_parameter(targetErrors, snr_values, 'targetErrors');
    maxFrames = expand_point_parameter(maxFrames, snr_values, 'maxFrames');

    scheme_labels = scheme_display_labels(schemes);
    num_schemes = numel(schemes);
    base_config = afdm_config();
    base_seed = base_config.simulation.random_seed;

    results.schemes = schemes;
    results.scheme_labels = scheme_labels;
    results.snr_values = snr_values;
    results.target_errors = targetErrors;
    results.max_frames = maxFrames;
    results.ber = zeros(numel(snr_values), num_schemes);
    results.error_bits = zeros(numel(snr_values), num_schemes);
    results.total_bits = zeros(numel(snr_values), num_schemes);
    results.frames_run = zeros(numel(snr_values), num_schemes);
    results.hit_target = false(numel(snr_values), num_schemes);

    fprintf('========== Adaptive BER comparison ==========\n');
    for snr_idx = 1:numel(snr_values)
        fprintf('SNR %.1f dB | target errors %d | max frames %d\n', ...
            snr_values(snr_idx), targetErrors(snr_idx), maxFrames(snr_idx));
    end
    fprintf('\n');

    for snr_idx = 1:numel(snr_values)
        snr_db = snr_values(snr_idx);

        for scheme_idx = 1:num_schemes
            cfg = apply_pre_chirp_scheme(base_config, schemes{scheme_idx});
            cfg.channel.snr_db = snr_db;

            total_err_bits = 0;
            total_bits = 0;
            frame_idx = 0;

            while total_err_bits < targetErrors(snr_idx) && frame_idx < maxFrames(snr_idx)
                frame_idx = frame_idx + 1;
                rng(base_seed + 1000000 * snr_idx + 1000 * scheme_idx + frame_idx, 'twister');

                [signal_cpp, ~, tx_bits, tx_state] = afdm_tx_engine(cfg);
                r_signal = multipath_channel(signal_cpp, cfg);
                r_signal = add_awgn(r_signal, cfg);
                [~, err_bits, num_bits] = afdm_rx_engine(r_signal, cfg, tx_bits, tx_state);

                total_err_bits = total_err_bits + err_bits;
                total_bits = total_bits + num_bits;

                if mod(frame_idx, 1000) == 0
                    fprintf('%-9s | SNR = %5.1f dB | frame %7d | errors = %7d | BER = %.3e\n', ...
                        scheme_labels{scheme_idx}, snr_db, frame_idx, total_err_bits, ...
                        total_err_bits / total_bits);
                end
            end

            results.error_bits(snr_idx, scheme_idx) = total_err_bits;
            results.total_bits(snr_idx, scheme_idx) = total_bits;
            results.frames_run(snr_idx, scheme_idx) = frame_idx;
            results.hit_target(snr_idx, scheme_idx) = total_err_bits >= targetErrors(snr_idx);
            results.ber(snr_idx, scheme_idx) = total_err_bits / total_bits;

            fprintf('%-9s | SNR = %5.1f dB | frames = %7d | errors = %7d / %-9d | BER = %.3e | target=%d\n', ...
                scheme_labels{scheme_idx}, snr_db, frame_idx, total_err_bits, total_bits, ...
                results.ber(snr_idx, scheme_idx), results.hit_target(snr_idx, scheme_idx));
        end
    end

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    results.output_dir = outputDir;
    results.timestamp = timestamp;

    figure('Name', 'Adaptive BER comparison', 'Color', 'w');
    semilogy(results.snr_values, results.ber, 'o-', 'LineWidth', 2);
    grid on;
    xlabel('SNR (dB)');
    ylabel('BER');
    title('Adaptive BER comparison');
    legend(results.scheme_labels, 'Location', 'southwest');
    saveas(gcf, fullfile(outputDir, ['adaptive_ber_comparison_' timestamp '.png']));

    matPath = fullfile(outputDir, ['adaptive_ber_comparison_' timestamp '.mat']);
    save(matPath, 'results');
    fprintf('\nSaved results to %s\n', matPath);
end

function values = default_target_errors(snr_values)
    values = zeros(size(snr_values));

    for idx = 1:numel(snr_values)
        if snr_values(idx) <= 20
            values(idx) = 500;
        elseif snr_values(idx) <= 22.5
            values(idx) = 200;
        else
            values(idx) = 100;
        end
    end
end

function values = default_max_frames(snr_values)
    values = zeros(size(snr_values));

    for idx = 1:numel(snr_values)
        if snr_values(idx) <= 20
            values(idx) = 2000;
        elseif snr_values(idx) <= 22.5
            values(idx) = 10000;
        elseif snr_values(idx) <= 25
            values(idx) = 50000;
        else
            values(idx) = 100000;
        end
    end
end

function values = expand_point_parameter(values, snr_values, name)
    if isscalar(values)
        values = repmat(values, size(snr_values));
    end

    if numel(values) ~= numel(snr_values)
        error('%s must be scalar or match the number of SNR points.', name);
    end

    values = reshape(values, size(snr_values));
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
