function results = run_delta_sweep(deltaRatios, numPaprFrames, numBerFrames, snr_values)
%RUN_DELTA_SWEEP Sweep proposed pre-chirp perturbation magnitudes.
%   results = run_delta_sweep(deltaRatios, numPaprFrames, numBerFrames, snr_values)

    rootDir = fileparts(mfilename('fullpath'));
    addpath(rootDir);
    setup_paths(rootDir);

    if nargin < 1 || isempty(deltaRatios)
        deltaRatios = [128, 64, 32, 16, 8, 4];
    end
    if nargin < 2 || isempty(numPaprFrames)
        numPaprFrames = 500;
    end
    if nargin < 3 || isempty(numBerFrames)
        numBerFrames = 20;
    end
    if nargin < 4 || isempty(snr_values)
        snr_values = [15, 20, 25];
    end

    base_config = afdm_config();
    base_seed = base_config.simulation.random_seed;
    base_c2 = base_config.pre_chirp.base_c2;
    frame_options.refresh_channel = base_config.simulation.refresh_channel_per_frame;

    labels = build_delta_labels(deltaRatios);
    num_delta = numel(deltaRatios);

    results.delta_ratios = deltaRatios;
    results.delta_values = base_c2 ./ deltaRatios;
    results.labels = labels;
    results.snr_values = snr_values;
    results.num_papr_frames = numPaprFrames;
    results.num_ber_frames = numBerFrames;
    results.papr_samples = zeros(numPaprFrames, num_delta);
    results.ber = zeros(numel(snr_values), num_delta);
    results.reference = run_reference_schemes(base_config, base_seed, numPaprFrames, numBerFrames, snr_values);

    fprintf('========== Proposed delta PAPR sampling ==========\n');
    for frame_idx = 1:numPaprFrames
        for delta_idx = 1:num_delta
            cfg = proposed_config(base_config, results.delta_values(delta_idx));
            frame = simulate_frame(cfg, base_seed + frame_idx, frame_options);
            results.papr_samples(frame_idx, delta_idx) = frame.papr;
        end

        if mod(frame_idx, 50) == 0 || frame_idx == numPaprFrames
            fprintf('PAPR frame %4d/%d complete\n', frame_idx, numPaprFrames);
        end
    end

    fprintf('\n========== Proposed delta BER sweep ==========\n');
    for snr_idx = 1:numel(snr_values)
        snr_db = snr_values(snr_idx);

        for delta_idx = 1:num_delta
            cfg = proposed_config(base_config, results.delta_values(delta_idx));
            cfg.channel.snr_db = snr_db;

            [ber_value, ~, ~] = run_ber_frames(cfg, base_seed, snr_idx, numBerFrames, frame_options);
            results.ber(snr_idx, delta_idx) = ber_value;

            fprintf('%-12s | SNR = %5.1f dB | BER = %.2e\n', ...
                labels{delta_idx}, snr_db, ber_value);
        end
    end

    results.summary = summarize_papr(results.papr_samples);
    results.reference.summary = summarize_papr(results.reference.papr_samples);

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    results.output_dir = outputDir;
    results.timestamp = timestamp;

    plot_delta_sweep_results(results);

    matPath = fullfile(outputDir, ['delta_sweep_' timestamp '.mat']);
    save(matPath, 'results');
    fprintf('\nSaved results to %s\n', matPath);
end

function cfg = proposed_config(base_config, delta)
    cfg = base_config;
    cfg.pre_chirp.delta = delta;
    cfg = apply_pre_chirp_scheme(cfg, 'proposed_grouping');
end

function reference = run_reference_schemes(base_config, base_seed, numPaprFrames, numBerFrames, snr_values)
    schemes = {'baseline', 'paper_grouping'};
    labels = {'Baseline', 'GPS'};
    num_schemes = numel(schemes);

    reference.schemes = schemes;
    reference.labels = labels;
    reference.papr_samples = zeros(numPaprFrames, num_schemes);
    reference.ber = zeros(numel(snr_values), num_schemes);
    frame_options.refresh_channel = base_config.simulation.refresh_channel_per_frame;

    fprintf('========== Reference PAPR sampling ==========\n');
    for frame_idx = 1:numPaprFrames
        for scheme_idx = 1:num_schemes
            cfg = apply_pre_chirp_scheme(base_config, schemes{scheme_idx});
            frame = simulate_frame(cfg, base_seed + frame_idx, frame_options);
            reference.papr_samples(frame_idx, scheme_idx) = frame.papr;
        end
    end

    fprintf('========== Reference BER sweep ==========\n');
    for snr_idx = 1:numel(snr_values)
        snr_db = snr_values(snr_idx);

        for scheme_idx = 1:num_schemes
            cfg = apply_pre_chirp_scheme(base_config, schemes{scheme_idx});
            cfg.channel.snr_db = snr_db;

            [ber_value, ~, ~] = run_ber_frames(cfg, base_seed, snr_idx, numBerFrames, frame_options);
            reference.ber(snr_idx, scheme_idx) = ber_value;

            fprintf('%-9s | SNR = %5.1f dB | BER = %.2e\n', ...
                labels{scheme_idx}, snr_db, ber_value);
        end
    end
end

function [ber_value, total_err_bits, total_bits] = run_ber_frames(cfg, base_seed, snr_idx, numFrames, frame_options)
    total_err_bits = 0;
    total_bits = 0;

    for frame_idx = 1:numFrames
        frame = simulate_frame(cfg, base_seed + 100000 * snr_idx + frame_idx, frame_options);
        total_err_bits = total_err_bits + frame.err_bits;
        total_bits = total_bits + frame.total_bits;
    end

    ber_value = total_err_bits / total_bits;
end

function summary = summarize_papr(papr_samples)
    summary.mean = mean(papr_samples, 1);
    summary.p90 = percentile_by_sort(papr_samples, 0.90);
    summary.p99 = percentile_by_sort(papr_samples, 0.99);
    summary.p999 = percentile_by_sort(papr_samples, 0.999);
end

function values = percentile_by_sort(samples, probability)
    sorted_samples = sort(samples, 1);
    idx = max(1, ceil(probability * size(sorted_samples, 1)));
    values = sorted_samples(idx, :);
end

function labels = build_delta_labels(deltaRatios)
    labels = cell(1, numel(deltaRatios));

    for idx = 1:numel(deltaRatios)
        labels{idx} = sprintf('base/%g', deltaRatios(idx));
    end
end

function plot_delta_sweep_results(results)
    outputDir = results.output_dir;
    timestamp = results.timestamp;
    x = 1:numel(results.delta_ratios);

    figure('Name', 'Delta sweep PAPR summary', 'Color', 'w');
    plot(x, results.summary.mean, 'o-', 'LineWidth', 2);
    hold on;
    plot(x, results.summary.p99, 's-', 'LineWidth', 2);
    yline(results.reference.summary.mean(1), '--', 'Baseline mean');
    yline(results.reference.summary.mean(2), '--', 'GPS mean');
    yline(results.reference.summary.p99(1), ':', 'Baseline p99');
    yline(results.reference.summary.p99(2), ':', 'GPS p99');
    grid on;
    xticks(x);
    xticklabels(results.labels);
    xlabel('Proposed delta');
    ylabel('PAPR (dB)');
    title('Proposed delta sweep PAPR');
    legend({'Proposed mean', 'Proposed p99'}, 'Location', 'best');
    saveas(gcf, fullfile(outputDir, ['delta_sweep_papr_' timestamp '.png']));

    figure('Name', 'Delta sweep BER', 'Color', 'w');
    semilogy(x, results.ber.', 'o-', 'LineWidth', 2);
    grid on;
    xticks(x);
    xticklabels(results.labels);
    xlabel('Proposed delta');
    ylabel('BER');
    title('Proposed delta sweep BER');
    legend(compose('SNR %.1f dB', results.snr_values), 'Location', 'best');
    saveas(gcf, fullfile(outputDir, ['delta_sweep_ber_' timestamp '.png']));
end
