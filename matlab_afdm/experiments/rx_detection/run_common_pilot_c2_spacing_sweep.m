function summary = run_common_pilot_c2_spacing_sweep(numFrames, deltaRatios, options)
%RUN_COMMON_PILOT_C2_SPACING_SWEEP Sweep local c2 spacing at one hard point.

    rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    setup_paths(rootDir);

    if nargin < 1 || isempty(numFrames)
        numFrames = 300;
    end
    if nargin < 2 || isempty(deltaRatios)
        deltaRatios = [1/64, 1/32, 1/16, 1/8, 1/4];
    end
    if nargin < 3
        options = struct();
    end

    N = get_option(options, 'N', 32);
    snrDb = get_option(options, 'snr_db', 20);
    fractionalDoppler = get_option(options, 'fractional_doppler', 0.4);
    seed = get_option(options, 'seed', 20260909);
    c2Base = sqrt(2) / (10 * N);
    patterns = [ ...
        1, 1, 1, 1; ...
        2, 3, 2, 3; ...
        3, 2, 3, 2; ...
        2, 2, 3, 3];

    numPoints = numel(deltaRatios);
    indexAccuracy = zeros(numPoints, 1);
    paprGainDb = zeros(numPoints, 1);
    berKnown = zeros(numPoints, 1);
    berDetected = zeros(numPoints, 1);
    crossDistance = zeros(numPoints, 1);
    receiverMargin = zeros(numPoints, 1);

    fprintf('========== common-pilot c2 spacing sweep ==========\n');
    fprintf('Frames/point=%d, SNR=%.1f dB, fractional=%.2f\n', ...
        numFrames, snrDb, fractionalDoppler);

    for pointIdx = 1:numPoints
        delta = c2Base * deltaRatios(pointIdx);
        c2Codebook = build_codebook(N, c2Base, delta, patterns);
        runOptions = options;
        runOptions.N = N;
        runOptions.seed = seed;
        runOptions.c2_codebook = c2Codebook;
        runOptions.save_outputs = false;

        result = run_common_pilot_c2_reconstruction_validation( ...
            numFrames, snrDb, fractionalDoppler, runOptions);
        indexAccuracy(pointIdx) = result.index_accuracy(1);
        paprGainDb(pointIdx) = result.mean_PAPR_gain_dB(1);
        berKnown(pointIdx) = result.BER_estimated_CSI_known_c2(1);
        berDetected(pointIdx) = result.BER_estimated_CSI_detected_c2(1);
        crossDistance(pointIdx) = result.cross_mode_distance(1);
        receiverMargin(pointIdx) = result.receiver_margin(1);
    end

    summary = table( ...
        deltaRatios(:), repmat(snrDb, numPoints, 1), ...
        repmat(fractionalDoppler, numPoints, 1), repmat(numFrames, numPoints, 1), ...
        indexAccuracy, paprGainDb, berKnown, berDetected, crossDistance, receiverMargin, ...
        'VariableNames', {'delta_over_c2', 'SNR_dB', 'fractional_doppler', ...
        'num_frames', 'index_accuracy', 'mean_PAPR_gain_dB', ...
        'BER_estimated_CSI_known_c2', 'BER_estimated_CSI_detected_c2', ...
        'cross_mode_distance', 'receiver_margin'});

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    csvPath = fullfile(outputDir, ['common_pilot_c2_spacing_' timestamp '.csv']);
    matPath = fullfile(outputDir, ['common_pilot_c2_spacing_' timestamp '.mat']);
    pngPath = fullfile(outputDir, ['common_pilot_c2_spacing_' timestamp '.png']);

    writetable(summary, csvPath);
    save(matPath, 'summary', 'patterns', 'options');
    plot_summary(summary, pngPath);

    fprintf('\nSaved spacing sweep to %s\n', csvPath);
    disp(summary);
end

function c2Codebook = build_codebook(N, c2Base, delta, patterns)
    c2Codebook = zeros(N, size(patterns, 1));
    for modeIdx = 1:size(patterns, 1)
        c2Codebook(:, modeIdx) = afdm.chirp.build_proposed_pattern( ...
            N, size(patterns, 2), patterns(modeIdx, :), c2Base, delta);
    end
end

function plot_summary(summary, pngPath)
    fig = figure('Name', 'c2 spacing tradeoff', 'Color', 'w');
    tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plot(summary.delta_over_c2, summary.index_accuracy, 'o-', 'LineWidth', 1.6);
    grid on;
    ylim([0, 1.02]);
    xlabel('\delta/c_2');
    ylabel('mode accuracy');
    title('Receiver distinguishability');

    nexttile;
    plot(summary.delta_over_c2, summary.mean_PAPR_gain_dB, 's-', 'LineWidth', 1.6);
    grid on;
    xlabel('\delta/c_2');
    ylabel('mean PAPR gain (dB)');
    title('Transmitter gain');

    nexttile;
    semilogy(summary.delta_over_c2, summary.BER_estimated_CSI_known_c2, ...
        'o-', 'LineWidth', 1.6);
    hold on;
    semilogy(summary.delta_over_c2, summary.BER_estimated_CSI_detected_c2, ...
        's-', 'LineWidth', 1.6);
    set(gca, 'YScale', 'log');
    grid on;
    xlabel('\delta/c_2');
    ylabel('BER');
    legend({'known mode', 'detected mode'}, 'Location', 'best');
    title('Mode-error cost');

    nexttile;
    yyaxis left;
    plot(summary.delta_over_c2, summary.cross_mode_distance, 'o-', 'LineWidth', 1.6);
    ylabel('cross-mode distance');
    yyaxis right;
    plot(summary.delta_over_c2, summary.receiver_margin, 's-', 'LineWidth', 1.6);
    ylabel('MMSE receiver margin');
    grid on;
    xlabel('\delta/c_2');
    title('Candidate geometry');

    exportgraphics(fig, pngPath, 'Resolution', 200);
    close(fig);
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = defaultValue;
    end
end
