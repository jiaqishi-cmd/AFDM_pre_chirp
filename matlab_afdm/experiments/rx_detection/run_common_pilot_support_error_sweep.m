function summary = run_common_pilot_support_error_sweep(numFrames, errorValues, options)
%RUN_COMMON_PILOT_SUPPORT_ERROR_SWEEP Test coarse Doppler-support mismatch.

    rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    setup_paths(rootDir);

    if nargin < 1 || isempty(numFrames)
        numFrames = 300;
    end
    if nargin < 2 || isempty(errorValues)
        errorValues = [0, 0.02, 0.05, 0.1, 0.2];
    end
    if nargin < 3
        options = struct();
    end

    N = get_option(options, 'N', 32);
    snrDb = get_option(options, 'snr_db', 20);
    fractionalDoppler = get_option(options, 'fractional_doppler', 0.4);
    seed = get_option(options, 'seed', 20260909);
    numPoints = numel(errorValues);

    gainNmse = zeros(numPoints, 1);
    effectiveNmse = zeros(numPoints, 1);
    indexAccuracy = zeros(numPoints, 1);
    berKnown = zeros(numPoints, 1);
    berDetected = zeros(numPoints, 1);

    fprintf('========== common-pilot support-error sweep ==========\n');
    fprintf('Frames/point=%d, SNR=%.1f dB, fractional=%.2f\n', ...
        numFrames, snrDb, fractionalDoppler);

    for pointIdx = 1:numPoints
        runOptions = options;
        runOptions.N = N;
        runOptions.seed = seed;
        runOptions.receiver_doppler_error = errorValues(pointIdx);
        runOptions.save_outputs = false;

        result = run_common_pilot_c2_reconstruction_validation( ...
            numFrames, snrDb, fractionalDoppler, runOptions);
        gainNmse(pointIdx) = result.path_gain_NMSE(1);
        effectiveNmse(pointIdx) = result.effective_channel_NMSE(1);
        indexAccuracy(pointIdx) = result.index_accuracy(1);
        berKnown(pointIdx) = result.BER_estimated_CSI_known_c2(1);
        berDetected(pointIdx) = result.BER_estimated_CSI_detected_c2(1);
    end

    summary = table( ...
        errorValues(:), repmat(snrDb, numPoints, 1), ...
        repmat(fractionalDoppler, numPoints, 1), repmat(numFrames, numPoints, 1), ...
        gainNmse, effectiveNmse, indexAccuracy, berKnown, berDetected, ...
        'VariableNames', {'receiver_doppler_error', 'SNR_dB', ...
        'fractional_doppler', 'num_frames', 'path_gain_NMSE', ...
        'effective_channel_NMSE', 'index_accuracy', ...
        'BER_estimated_CSI_known_c2', 'BER_estimated_CSI_detected_c2'});

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    csvPath = fullfile(outputDir, ['common_pilot_support_error_' timestamp '.csv']);
    matPath = fullfile(outputDir, ['common_pilot_support_error_' timestamp '.mat']);
    pngPath = fullfile(outputDir, ['common_pilot_support_error_' timestamp '.png']);

    writetable(summary, csvPath);
    save(matPath, 'summary', 'options');
    plot_summary(summary, pngPath);

    fprintf('\nSaved support-error sweep to %s\n', csvPath);
    disp(summary);
end

function plot_summary(summary, pngPath)
    fig = figure('Name', 'Doppler support error', 'Color', 'w');
    tiledlayout(fig, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    semilogy(summary.receiver_doppler_error, summary.path_gain_NMSE, ...
        'o-', 'LineWidth', 1.6);
    hold on;
    semilogy(summary.receiver_doppler_error, summary.effective_channel_NMSE, ...
        's-', 'LineWidth', 1.6);
    set(gca, 'YScale', 'log');
    grid on;
    xlabel('Doppler support error (bins)');
    ylabel('NMSE');
    legend({'path gains', 'effective channel'}, 'Location', 'best');
    title('Reconstruction mismatch');

    nexttile;
    plot(summary.receiver_doppler_error, summary.index_accuracy, ...
        'o-', 'LineWidth', 1.6);
    grid on;
    ylim([0, 1.02]);
    xlabel('Doppler support error (bins)');
    ylabel('mode accuracy');
    title('Mode detection');

    nexttile;
    semilogy(summary.receiver_doppler_error, summary.BER_estimated_CSI_known_c2, ...
        'o-', 'LineWidth', 1.6);
    hold on;
    semilogy(summary.receiver_doppler_error, summary.BER_estimated_CSI_detected_c2, ...
        's-', 'LineWidth', 1.6);
    set(gca, 'YScale', 'log');
    grid on;
    xlabel('Doppler support error (bins)');
    ylabel('BER');
    legend({'known mode', 'detected mode'}, 'Location', 'best');
    title('Receiver robustness');

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
