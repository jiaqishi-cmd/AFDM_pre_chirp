function run_caseA_theta_scan_resume()
%RUN_CASEA_THETA_SCAN_RESUME Resume Case A theta scan from saved progress.

    rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    setup_paths(rootDir);

    outputDir = fullfile(fileparts(rootDir), 'results');
    progressPath = fullfile(outputDir, 'results_caseA_theta_scan_progress.mat');
    if ~exist(progressPath, 'file')
        error('Progress file not found: %s', progressPath);
    end

    loaded = load(progressPath);
    theta_list = loaded.theta_list;
    SNR_dB = loaded.SNR_dB;
    schemes = loaded.schemes;
    ber = loaded.ber;
    ber_plot = loaded.ber_plot;
    error_count = loaded.error_count;
    total_bits = loaded.total_bits;
    frames_used = loaded.frames_used;
    minErrTarget = loaded.minErrTarget;
    maxBits = loaded.maxBits;
    startIdx = loaded.thetaIdx + 1;

    N = 64;
    V = 4;
    alpha_max = 3;
    c1 = (2 * alpha_max + 1) / (2 * N);
    c2_base = sqrt(2) / (10 * N);
    gps_pattern = [2 2 1 1];
    proposed_pattern = gps_pattern;
    proposed_delta = c2_base / 16;
    l1 = 0;
    alpha1 = 0;
    l2 = 2;
    alpha2 = 2;

    [c2_gps, ~] = afdm.chirp.build_gps_pattern(N, V, gps_pattern);
    [c2_prop, ~] = afdm.chirp.build_proposed_pattern(N, V, proposed_pattern, c2_base, proposed_delta);
    c2_values = {c2_base, c2_gps, c2_prop};
    numTheta = numel(theta_list);
    numSchemes = numel(schemes);
    maxFrames = ceil(maxBits / 64) + 10;

    fprintf('========== Resume Case A theta-scan BER ==========\n');
    fprintf('Resuming from theta index %d/%d, SNR=%g dB\n', startIdx, numTheta, SNR_dB);

    for idx = startIdx:numTheta
        theta = theta_list(idx);
        h = [1, exp(1i * theta)] / sqrt(2);

        for schemeIdx = 1:numSchemes
            cfg = build_caseA_config(N, c1, c2_values{schemeIdx}, SNR_dB, ...
                [l1, l2], [alpha1, alpha2], h);

            err = 0;
            bits = 0;
            frameIdx = 0;
            while err < minErrTarget && bits < maxBits && frameIdx < maxFrames
                frameIdx = frameIdx + 1;
                seed = 1 + 1000000 * idx + 10000 * schemeIdx + frameIdx;
                frame = simulate_frame(cfg, seed);
                err = err + frame.err_bits;
                bits = bits + frame.total_bits;
            end

            error_count(idx, schemeIdx) = err;
            total_bits(idx, schemeIdx) = bits;
            frames_used(idx, schemeIdx) = frameIdx;
            ber(idx, schemeIdx) = err / max(bits, 1);
            if err == 0
                ber_plot(idx, schemeIdx) = 0.5 / max(bits, 1);
            else
                ber_plot(idx, schemeIdx) = ber(idx, schemeIdx);
            end
        end

        thetaIdx = idx; %#ok<NASGU>
        save(progressPath, ...
            'theta_list', 'SNR_dB', 'schemes', 'ber', 'ber_plot', ...
            'error_count', 'total_bits', 'frames_used', ...
            'minErrTarget', 'maxBits', 'gps_pattern', 'thetaIdx');

        fprintf('theta/pi=%6.3f | base %.3e (%d/%d) | GPS %.3e (%d/%d) | prop %.3e (%d/%d)\n', ...
            theta / pi, ...
            ber(idx, 1), error_count(idx, 1), total_bits(idx, 1), ...
            ber(idx, 2), error_count(idx, 2), total_bits(idx, 2), ...
            ber(idx, 3), error_count(idx, 3), total_bits(idx, 3));
    end

    ber_base = ber(:, 1);
    ber_gps = ber(:, 2);
    ber_prop = ber(:, 3);
    ber_base_plot = ber_plot(:, 1);
    ber_gps_plot = ber_plot(:, 2);
    ber_prop_plot = ber_plot(:, 3);

    save(fullfile(outputDir, 'results_caseA_theta_scan.mat'), ...
        'theta_list', 'SNR_dB', 'schemes', 'ber', 'ber_plot', ...
        'ber_base', 'ber_gps', 'ber_prop', ...
        'ber_base_plot', 'ber_gps_plot', 'ber_prop_plot', ...
        'error_count', 'total_bits', 'frames_used', ...
        'minErrTarget', 'maxBits', 'gps_pattern');

    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    figure('Name', 'Case A theta scan BER', 'Color', 'w');
    semilogy(theta_list / pi, ber_base_plot, 'LineWidth', 2); hold on;
    semilogy(theta_list / pi, ber_gps_plot, 'LineWidth', 2);
    semilogy(theta_list / pi, ber_prop_plot, 'LineWidth', 2);
    xline(1, '--k', '\theta/\pi=1', 'LineWidth', 1.2);
    grid on;
    xlim([0, 2]);
    xlabel('\theta / \pi');
    ylabel('BER');
    title(sprintf('BER vs path phase theta on Case A fixed two-path channel, SNR=%g dB', SNR_dB));
    legend({'baseline', 'GPS', 'proposed'}, 'Location', 'best');
    saveas(gcf, fullfile(outputDir, ['caseA_theta_scan_BER_' timestamp '.png']));

    [worstGpsBer, worstIdx] = max(ber_gps);
    [~, piIdx] = min(abs(theta_list - pi));
    fprintf('\nWorst GPS theta/pi = %.3f, BER = %.3e\n', theta_list(worstIdx) / pi, worstGpsBer);
    fprintf('At theta=pi: baseline=%.3e, GPS=%.3e, proposed=%.3e\n', ...
        ber_base(piIdx), ber_gps(piIdx), ber_prop(piIdx));
    fprintf('Saved theta scan to %s\n', fullfile(outputDir, 'results_caseA_theta_scan.mat'));

    run('experiments/caseA/plot_caseA_theta_scan_enhanced.m');
end

function cfg = build_caseA_config(N, c1, c2, snrDb, delayTaps, dopplerTaps, gains)
    cfg = afdm_config();
    cfg.waveform.NumSubcarriers = N;
    cfg.waveform.CPPLength = max(delayTaps);
    cfg.waveform.c1 = c1;
    cfg.waveform.c2 = c2;
    cfg.modulation.M_mod = 2;
    cfg.modulation.modType = 'psk';
    cfg.channel.multipath = true;
    cfg.channel.add_noise = true;
    cfg.channel.snr_db = snrDb;
    cfg.channel.delay_taps = delayTaps;
    cfg.channel.doppler_taps = dopplerTaps;
    cfg.channel.doppler_freq = dopplerTaps / N;
    cfg.channel.chan_coef = gains;
    cfg.pre_chirp.scheme = 'baseline';
    cfg.pre_chirp.profile.scheme = 'baseline';
    cfg.pre_chirp.profile.c2 = c2;
end
