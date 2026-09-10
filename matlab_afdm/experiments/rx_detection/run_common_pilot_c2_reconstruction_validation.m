function results = run_common_pilot_c2_reconstruction_validation(numFrames, snrValues, fractionalValues, options)
%RUN_COMMON_PILOT_C2_RECONSTRUCTION_VALIDATION Test a practical c2 receive chain.
%
% A mode-independent training block uses a fixed reference c2. Given coarse
% delay/Doppler support, the receiver estimates physical path gains by LS,
% reconstructs one effective channel per candidate c2, detects the selected
% mode, and finally equalizes the data. The experiment separates losses due
% to channel estimation from losses due to c2-mode detection.

    rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    setup_paths(rootDir);

    if nargin < 1 || isempty(numFrames)
        numFrames = 100;
    end
    if nargin < 2 || isempty(snrValues)
        snrValues = [15, 20, 25];
    end
    if nargin < 3 || isempty(fractionalValues)
        fractionalValues = [0, 0.2, 0.4];
    end
    if nargin < 4
        options = struct();
    end

    cfg = configure_base_config(options);
    seed = get_option(options, 'seed', cfg.simulation.random_seed);
    c2Base = cfg.pre_chirp.base_c2;
    c2Codebook = get_option(options, 'c2_codebook', build_default_codebook(cfg));
    if size(c2Codebook, 1) ~= cfg.waveform.NumSubcarriers
        error('c2_codebook must have N rows and one column per mode.');
    end
    pilotC2 = get_option(options, 'pilot_c2', c2Base);
    pilotSnrOffsetDb = get_option(options, 'pilot_snr_offset_db', 0);
    receiverDopplerError = get_option(options, 'receiver_doppler_error', 0);
    saveOutputs = get_option(options, 'save_outputs', true);

    numSnr = numel(snrValues);
    numFractional = numel(fractionalValues);
    numRows = numSnr * numFractional;
    row = 0;

    resultFractional = zeros(numRows, 1);
    resultSnr = zeros(numRows, 1);
    resultReceiverDopplerError = repmat(receiverDopplerError, numRows, 1);
    resultFrames = repmat(numFrames, numRows, 1);
    pilotCondition = zeros(numRows, 1);
    gainNmse = zeros(numRows, 1);
    effectiveNmse = zeros(numRows, 1);
    minModeSeparation = zeros(numRows, 1);
    crossModeDistance = zeros(numRows, 1);
    receiverMargin = zeros(numRows, 1);
    receiverMarginCorrect = NaN(numRows, 1);
    receiverMarginError = NaN(numRows, 1);
    indexAccuracy = zeros(numRows, 1);
    berPerfect = zeros(numRows, 1);
    berEstimatedKnown = zeros(numRows, 1);
    berEstimatedDetected = zeros(numRows, 1);
    paprGainDb = zeros(numRows, 1);
    numModes = size(c2Codebook, 2);
    confusion = zeros(numFractional, numSnr, numModes, numModes);

    fprintf('========== common-pilot c2 reconstruction validation ==========\n');
    fprintf('Frames/point=%d, N=%d, M=%d\n', ...
        numFrames, cfg.waveform.NumSubcarriers, cfg.modulation.M_mod);
    fprintf('c2 codebook: %d grouped-vector modes\n', numModes);
    fprintf('fractional Doppler magnitudes: %s\n', mat2str(fractionalValues, 4));

    oldRngState = rng;
    cleanup = onCleanup(@() rng(oldRngState));

    pilotSymbols = build_common_pilot(cfg.waveform.NumSubcarriers);

    for fracIdx = 1:numFractional
        cfgScenario = apply_fractional_profile(cfg, fractionalValues(fracIdx));
        cfgReceiver = apply_receiver_support_error(cfgScenario, receiverDopplerError, options);
        precomputed = precompute_models( ...
            cfgScenario, cfgReceiver, pilotSymbols, pilotC2, c2Codebook);

        for snrIdx = 1:numSnr
            row = row + 1;
            dataSnrDb = snrValues(snrIdx);
            pilotSnrDb = dataSnrDb + pilotSnrOffsetDb;
            accum = initialize_accumulator();

            fprintf('\nfractional=%.3f, data/pilot SNR=%.1f/%.1f dB\n', ...
                fractionalValues(fracIdx), dataSnrDb, pilotSnrDb);

            for frameIdx = 1:numFrames
                rng(seed + 1000000 * fracIdx + 10000 * snrIdx + frameIdx, 'twister');
                frame = simulate_frame( ...
                    cfgScenario, precomputed, c2Codebook, dataSnrDb, pilotSnrDb);

                accum.gain_nmse = accum.gain_nmse + frame.gain_nmse;
                accum.effective_nmse = accum.effective_nmse + frame.effective_nmse;
                accum.min_mode_separation = accum.min_mode_separation + frame.min_mode_separation;
                accum.cross_mode_distance = accum.cross_mode_distance + frame.cross_mode_distance;
                accum.receiver_margin = accum.receiver_margin + frame.receiver_margin;
                accum.index_correct = accum.index_correct + double(frame.detected_index == frame.tx_index);
                if frame.detected_index == frame.tx_index
                    accum.margin_correct = accum.margin_correct + frame.receiver_margin;
                    accum.num_correct = accum.num_correct + 1;
                else
                    accum.margin_error = accum.margin_error + frame.receiver_margin;
                    accum.num_error = accum.num_error + 1;
                end
                accum.err_perfect = accum.err_perfect + frame.perfect.err_bits;
                accum.err_estimated_known = accum.err_estimated_known + frame.estimated_known.err_bits;
                accum.err_estimated_detected = accum.err_estimated_detected + frame.estimated_detected.err_bits;
                accum.total_bits = accum.total_bits + frame.perfect.total_bits;
                accum.papr_gain_db = accum.papr_gain_db + frame.papr_gain_db;
                confusion(fracIdx, snrIdx, frame.tx_index, frame.detected_index) = ...
                    confusion(fracIdx, snrIdx, frame.tx_index, frame.detected_index) + 1;
            end

            resultFractional(row) = fractionalValues(fracIdx);
            resultSnr(row) = dataSnrDb;
            pilotCondition(row) = cond(precomputed.pilot_basis_rx);
            gainNmse(row) = accum.gain_nmse / numFrames;
            effectiveNmse(row) = accum.effective_nmse / numFrames;
            minModeSeparation(row) = accum.min_mode_separation / numFrames;
            crossModeDistance(row) = accum.cross_mode_distance / numFrames;
            receiverMargin(row) = accum.receiver_margin / numFrames;
            if accum.num_correct > 0
                receiverMarginCorrect(row) = accum.margin_correct / accum.num_correct;
            end
            if accum.num_error > 0
                receiverMarginError(row) = accum.margin_error / accum.num_error;
            end
            indexAccuracy(row) = accum.index_correct / numFrames;
            berPerfect(row) = accum.err_perfect / accum.total_bits;
            berEstimatedKnown(row) = accum.err_estimated_known / accum.total_bits;
            berEstimatedDetected(row) = accum.err_estimated_detected / accum.total_bits;
            paprGainDb(row) = accum.papr_gain_db / numFrames;

            fprintf('  gain/effective NMSE: %.3e / %.3e\n', gainNmse(row), effectiveNmse(row));
            fprintf('  min mode separation: %.4f, index accuracy: %.4f\n', ...
                minModeSeparation(row), indexAccuracy(row));
            fprintf('  cross-mode distance / receiver margin: %.4f / %.4e\n', ...
                crossModeDistance(row), receiverMargin(row));
            fprintf('  margin on correct/error frames: %.4e / %.4e\n', ...
                receiverMarginCorrect(row), receiverMarginError(row));
            fprintf('  BER perfect / estimated-known / estimated-detected: %.3e / %.3e / %.3e\n', ...
                berPerfect(row), berEstimatedKnown(row), berEstimatedDetected(row));
        end
    end

    results = table( ...
        resultFractional, resultSnr, resultReceiverDopplerError, resultFrames, ...
        pilotCondition, gainNmse, ...
        effectiveNmse, minModeSeparation, crossModeDistance, receiverMargin, ...
        receiverMarginCorrect, receiverMarginError, indexAccuracy, berPerfect, ...
        berEstimatedKnown, berEstimatedDetected, paprGainDb, ...
        'VariableNames', {'fractional_doppler', 'SNR_dB', 'receiver_doppler_error', ...
        'num_frames', ...
        'pilot_basis_condition', 'path_gain_NMSE', 'effective_channel_NMSE', ...
        'min_mode_separation', 'cross_mode_distance', 'receiver_margin', ...
        'receiver_margin_correct', 'receiver_margin_error', 'index_accuracy', ...
        'BER_perfect_CSI', ...
        'BER_estimated_CSI_known_c2', 'BER_estimated_CSI_detected_c2', ...
        'mean_PAPR_gain_dB'});

    if saveOutputs
        outputDir = fullfile(fileparts(rootDir), 'results');
        if ~exist(outputDir, 'dir')
            mkdir(outputDir);
        end
        timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
        csvPath = fullfile(outputDir, ['common_pilot_c2_validation_' timestamp '.csv']);
        matPath = fullfile(outputDir, ['common_pilot_c2_validation_' timestamp '.mat']);
        pngPath = fullfile(outputDir, ['common_pilot_c2_validation_' timestamp '.png']);
        logPath = fullfile(outputDir, ['common_pilot_c2_validation_' timestamp '.log']);

        writetable(results, csvPath);
        save(matPath, 'results', 'confusion', 'c2Codebook', 'cfg', ...
            'snrValues', 'fractionalValues', 'pilotC2', 'pilotSnrOffsetDb', ...
            'receiverDopplerError');
        plot_results(results, fractionalValues, pngPath);
        write_log(logPath, c2Codebook, pilotC2, pilotSnrOffsetDb);
        fprintf('\nSaved results to %s\n', csvPath);
    end
    disp(results);
end

function c2Codebook = build_default_codebook(cfg)
    N = cfg.waveform.NumSubcarriers;
    numGroups = 4;
    delta = cfg.pre_chirp.base_c2 / 16;
    patterns = [ ...
        1, 1, 1, 1; ...
        2, 3, 2, 3; ...
        3, 2, 3, 2; ...
        2, 2, 3, 3];
    c2Codebook = zeros(N, size(patterns, 1));
    for modeIdx = 1:size(patterns, 1)
        c2Codebook(:, modeIdx) = afdm.chirp.build_proposed_pattern( ...
            N, numGroups, patterns(modeIdx, :), cfg.pre_chirp.base_c2, delta);
    end
end

function cfg = configure_base_config(options)
    cfg = afdm_config();
    N = get_option(options, 'N', 32);

    cfg.waveform.NumSubcarriers = N;
    cfg.waveform.CPPLength = get_option(options, 'cpp_length', 3);
    cfg.waveform.c1 = 7 / (2 * N);
    cfg.modulation.M_mod = get_option(options, 'M_mod', 16);
    cfg.modulation.modType = get_option(options, 'modType', 'qam');
    cfg.channel.add_noise = true;
    cfg.channel.multipath = true;
    cfg.channel.delay_taps = get_option(options, 'delay_taps', [0, 1, 2]);
    cfg.channel.integer_doppler_taps = get_option(options, 'integer_doppler_taps', [-2, 0, 2]);
    cfg.channel.fractional_pattern = get_option(options, 'fractional_pattern', [1, -0.5, 0.75]);
    cfg.receiver.detector = 'full_mmse';

    c2Base = sqrt(2) / (10 * N);
    cfg.pre_chirp.base_c2 = c2Base;
    cfg.waveform.c2 = c2Base;
end

function cfg = apply_fractional_profile(cfg, fractionalMagnitude)
    dopplerTaps = cfg.channel.integer_doppler_taps + ...
        fractionalMagnitude * cfg.channel.fractional_pattern;
    cfg.channel.doppler_taps = dopplerTaps;
    cfg.channel.doppler_freq = dopplerTaps / cfg.waveform.NumSubcarriers;
end

function cfgReceiver = apply_receiver_support_error(cfg, errorMagnitude, options)
    cfgReceiver = cfg;
    defaultPattern = cfg.channel.fractional_pattern;
    errorPattern = get_option(options, 'receiver_doppler_error_pattern', defaultPattern);
    if numel(errorPattern) ~= numel(cfg.channel.doppler_taps)
        error('receiver_doppler_error_pattern must match the number of paths.');
    end
    cfgReceiver.channel.doppler_taps = cfg.channel.doppler_taps + ...
        errorMagnitude * errorPattern;
    cfgReceiver.channel.doppler_freq = ...
        cfgReceiver.channel.doppler_taps / cfg.waveform.NumSubcarriers;
end

function pilotSymbols = build_common_pilot(N)
    n = (0:N-1).';
    pilotSymbols = exp(1i * (pi / 2) * mod(n.^2 + 3 * n, 4));
end

function precomputed = precompute_models(cfgTrue, cfgRx, pilotSymbols, pilotC2, c2Codebook)
    numPaths = numel(cfgTrue.channel.delay_taps);
    numModes = size(c2Codebook, 2);
    N = cfgTrue.waveform.NumSubcarriers;

    pilotSignal = afdm.tx.idaft_mod(pilotSymbols, N, cfgTrue.waveform.c1, pilotC2);
    pilotSignalCpp = afdm.tx.add_cpp( ...
        pilotSignal, cfgTrue.waveform.CPPLength, cfgTrue.waveform.c1);
    precomputed.pilot_basis_true = zeros(numel(pilotSignalCpp), numPaths);
    precomputed.pilot_basis_rx = zeros(numel(pilotSignalCpp), numPaths);
    precomputed.effective_basis_true = cell(numModes, numPaths);
    precomputed.effective_basis_rx = cell(numModes, numPaths);

    for pathIdx = 1:numPaths
        unitGains = zeros(1, numPaths);
        unitGains(pathIdx) = 1;
        cfgTruePath = cfgTrue;
        cfgTruePath.channel.chan_coef = unitGains;
        cfgRxPath = cfgRx;
        cfgRxPath.channel.chan_coef = unitGains;
        precomputed.pilot_basis_true(:, pathIdx) = ...
            afdm.channel.multipath(pilotSignalCpp, cfgTruePath);
        precomputed.pilot_basis_rx(:, pathIdx) = ...
            afdm.channel.multipath(pilotSignalCpp, cfgRxPath);

        for modeIdx = 1:numModes
            precomputed.effective_basis_true{modeIdx, pathIdx} = ...
                afdm.rx.estimate_effective_channel( ...
                N, cfgTrue.waveform.c1, c2Codebook(:, modeIdx), unitGains, ...
                cfgTrue.channel.delay_taps, cfgTrue.channel.doppler_freq);
            precomputed.effective_basis_rx{modeIdx, pathIdx} = ...
                afdm.rx.estimate_effective_channel( ...
                N, cfgRx.waveform.c1, c2Codebook(:, modeIdx), unitGains, ...
                cfgRx.channel.delay_taps, cfgRx.channel.doppler_freq);
        end
    end
end

function frame = simulate_frame(cfg, precomputed, c2Codebook, dataSnrDb, pilotSnrDb)
    N = cfg.waveform.NumSubcarriers;
    numPaths = numel(cfg.channel.delay_taps);
    numModes = size(c2Codebook, 2);

    gains = (randn(numPaths, 1) + 1i * randn(numPaths, 1)) / sqrt(2);
    gains = gains / norm(gains);
    cfgFrame = cfg;
    cfgFrame.channel.chan_coef = gains.';

    pilotClean = precomputed.pilot_basis_true * gains;
    [pilotReceived, ~] = add_noise_at_snr(pilotClean, pilotSnrDb);
    gainsEstimated = precomputed.pilot_basis_rx \ pilotReceived;

    [symbols, txBits] = afdm.tx.random_data( ...
        N, cfg.modulation.M_mod, cfg.modulation.modType);
    [dataSignalCpp, paprValues, txIdx] = select_low_papr_signal( ...
        symbols, cfg, c2Codebook);
    dataClean = afdm.channel.multipath(dataSignalCpp, cfgFrame);
    [dataReceived, noiseVar] = add_noise_at_snr(dataClean, dataSnrDb);
    dataNoCpp = afdm.rx.remove_cpp(dataReceived, cfg.waveform.CPPLength);

    hTrue = cell(numModes, 1);
    hEstimated = cell(numModes, 1);
    for modeIdx = 1:numModes
        hTrue{modeIdx} = combine_effective_basis( ...
            precomputed.effective_basis_true(modeIdx, :), gains);
        hEstimated{modeIdx} = combine_effective_basis( ...
            precomputed.effective_basis_rx(modeIdx, :), gainsEstimated);
    end

    xPerfect = detect_symbols(dataNoCpp, hTrue{txIdx}, c2Codebook(:, txIdx), noiseVar, cfg);
    xEstimatedKnown = detect_symbols( ...
        dataNoCpp, hEstimated{txIdx}, c2Codebook(:, txIdx), noiseVar, cfg);

    metrics = zeros(numModes, 1);
    candidateSymbols = cell(numModes, 1);
    for modeIdx = 1:numModes
        [candidateSymbols{modeIdx}, metrics(modeIdx)] = detect_symbols( ...
            dataNoCpp, hEstimated{modeIdx}, c2Codebook(:, modeIdx), noiseVar, cfg);
    end
    [~, detectedIdx] = min(metrics);

    effectiveError = 0;
    for modeIdx = 1:numModes
        effectiveError = effectiveError + ...
            norm(hEstimated{modeIdx} - hTrue{modeIdx}, 'fro')^2 / ...
            max(norm(hTrue{modeIdx}, 'fro')^2, eps);
    end

    frame.tx_index = txIdx;
    frame.detected_index = detectedIdx;
    frame.gain_nmse = norm(gainsEstimated - gains)^2 / max(norm(gains)^2, eps);
    frame.effective_nmse = effectiveError / numModes;
    frame.min_mode_separation = minimum_mode_separation(hTrue);
    [frame.cross_mode_distance, frame.receiver_margin] = receiver_mode_metrics( ...
        hTrue, c2Codebook, txIdx, noiseVar);
    frame.papr_gain_db = paprValues(find_baseline_index(c2Codebook, cfg.pre_chirp.base_c2)) ...
        - paprValues(txIdx);
    frame.perfect = count_errors(xPerfect, txBits, cfg);
    frame.estimated_known = count_errors(xEstimatedKnown, txBits, cfg);
    frame.estimated_detected = count_errors(candidateSymbols{detectedIdx}, txBits, cfg);
end

function [minCrossDistance, receiverMargin] = receiver_mode_metrics( ...
        channels, c2Codebook, txIdx, noiseVar)
    numModes = numel(channels);
    N = size(channels{1}, 1);
    n = (0:N-1).';
    minCrossDistance = inf;
    expectedMse = zeros(numModes, 1);

    for txMode = 1:numModes
        for rxMode = 1:numModes
            if txMode == rxMode
                continue;
            end
            phase = exp(-1i * 2 * pi * ...
                (c2Codebook(:, rxMode) - c2Codebook(:, txMode)) .* (n.^2));
            crossChannel = phase .* channels{txMode};
            mismatch = crossChannel - channels{rxMode};
            distance = norm(mismatch, 'fro') / max(norm(crossChannel, 'fro'), eps);
            minCrossDistance = min(minCrossDistance, distance);
        end
    end

    for rxMode = 1:numModes
        phase = exp(-1i * 2 * pi * ...
            (c2Codebook(:, rxMode) - c2Codebook(:, txIdx)) .* (n.^2));
        actualChannel = phase .* channels{txIdx};
        modelChannel = channels{rxMode};
        equalizer = (modelChannel' * modelChannel + noiseVar * eye(N)) \ modelChannel';
        response = equalizer * actualChannel;
        expectedMse(rxMode) = norm(response - eye(N), 'fro')^2 / N + ...
            noiseVar * norm(equalizer, 'fro')^2 / N;
    end

    wrongModes = setdiff(1:numModes, txIdx);
    receiverMargin = min(expectedMse(wrongModes)) - expectedMse(txIdx);
end

function [xHat, metric] = detect_symbols(dataNoCpp, hEff, c2Rx, noiseVar, cfg)
    yDaft = afdm.rx.daft_demod( ...
        dataNoCpp, cfg.waveform.NumSubcarriers, cfg.waveform.c1, c2Rx);
    xHat = afdm.rx.mmse_equalize(yDaft, hEff, noiseVar);
    xSliced = slice_symbols(xHat, cfg.modulation.M_mod, cfg.modulation.modType);

    observationResidual = norm(yDaft - hEff * xSliced)^2 / max(norm(yDaft)^2, eps);
    constellationResidual = mean(abs(xHat - xSliced).^2);
    metric = observationResidual + constellationResidual;
end

function hEff = combine_effective_basis(pathBasis, gains)
    hEff = zeros(size(pathBasis{1}));
    for pathIdx = 1:numel(gains)
        hEff = hEff + gains(pathIdx) * pathBasis{pathIdx};
    end
end

function separation = minimum_mode_separation(channels)
    separation = inf;
    for firstIdx = 1:numel(channels)-1
        for secondIdx = firstIdx+1:numel(channels)
            denominator = sqrt( ...
                norm(channels{firstIdx}, 'fro') * norm(channels{secondIdx}, 'fro'));
            pairDistance = norm(channels{firstIdx} - channels{secondIdx}, 'fro') / ...
                max(denominator, eps);
            separation = min(separation, pairDistance);
        end
    end
end

function [noisySignal, noiseVar] = add_noise_at_snr(cleanSignal, snrDb)
    noiseVar = mean(abs(cleanSignal).^2) / (10^(snrDb / 10));
    noise = sqrt(noiseVar / 2) * ...
        (randn(size(cleanSignal)) + 1i * randn(size(cleanSignal)));
    noisySignal = cleanSignal + noise;
end

function xSliced = slice_symbols(xHat, MMod, modType)
    switch lower(modType)
        case 'qam'
            idx = qamdemod(xHat, MMod, 'OutputType', 'integer', 'UnitAveragePower', true);
            xSliced = qammod(idx, MMod, 'UnitAveragePower', true);
        case 'psk'
            idx = pskdemod(xHat, MMod, pi/MMod, 'OutputType', 'integer');
            xSliced = pskmod(idx, MMod, pi/MMod);
        otherwise
            error('Unsupported modulation type: %s.', modType);
    end
end

function out = count_errors(xHat, txBits, cfg)
    rxBits = afdm.rx.symbol_decision( ...
        xHat, cfg.modulation.M_mod, cfg.modulation.modType);
    [out.err_bits, out.total_bits] = afdm.rx.compute_bit_errors(rxBits, txBits);
end

function [signalCpp, paprValues, selectedIdx] = select_low_papr_signal(symbols, cfg, c2Codebook)
    numModes = size(c2Codebook, 2);
    signals = cell(numModes, 1);
    paprValues = zeros(numModes, 1);

    for modeIdx = 1:numModes
        signal = afdm.tx.idaft_mod( ...
            symbols, cfg.waveform.NumSubcarriers, cfg.waveform.c1, c2Codebook(:, modeIdx));
        signals{modeIdx} = afdm.tx.add_cpp( ...
            signal, cfg.waveform.CPPLength, cfg.waveform.c1);
        paprValues(modeIdx) = afdm.tx.compute_papr(signal);
    end

    [~, selectedIdx] = min(paprValues);
    signalCpp = signals{selectedIdx};
end

function baselineIdx = find_baseline_index(c2Codebook, c2Base)
    [~, baselineIdx] = min(vecnorm(c2Codebook - c2Base, 2, 1));
end

function accum = initialize_accumulator()
    accum = struct( ...
        'gain_nmse', 0, ...
        'effective_nmse', 0, ...
        'min_mode_separation', 0, ...
        'cross_mode_distance', 0, ...
        'receiver_margin', 0, ...
        'margin_correct', 0, ...
        'margin_error', 0, ...
        'num_correct', 0, ...
        'num_error', 0, ...
        'index_correct', 0, ...
        'err_perfect', 0, ...
        'err_estimated_known', 0, ...
        'err_estimated_detected', 0, ...
        'total_bits', 0, ...
        'papr_gain_db', 0);
end

function plot_results(results, fractionalValues, pngPath)
    fig = figure('Name', 'common-pilot c2 validation', 'Color', 'w');
    tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    colors = lines(numel(fractionalValues));

    nexttile;
    hold on;
    for idx = 1:numel(fractionalValues)
        rows = abs(results.fractional_doppler - fractionalValues(idx)) < 1e-12;
        plot(results.SNR_dB(rows), results.index_accuracy(rows), 'o-', ...
            'Color', colors(idx, :), 'LineWidth', 1.5);
    end
    grid on;
    ylim([0, 1.02]);
    xlabel('SNR (dB)');
    ylabel('mode accuracy');
    title('c2 mode detection');

    nexttile;
    hold on;
    for idx = 1:numel(fractionalValues)
        rows = abs(results.fractional_doppler - fractionalValues(idx)) < 1e-12;
        semilogy(results.SNR_dB(rows), max(results.path_gain_NMSE(rows), eps), 'o-', ...
            'Color', colors(idx, :), 'LineWidth', 1.5);
    end
    set(gca, 'YScale', 'log');
    ytickformat('%.1e');
    grid on;
    xlabel('SNR (dB)');
    ylabel('path-gain NMSE');
    title('common training LS');

    nexttile;
    hold on;
    for idx = 1:numel(fractionalValues)
        rows = abs(results.fractional_doppler - fractionalValues(idx)) < 1e-12;
        semilogy(results.SNR_dB(rows), max(results.BER_estimated_CSI_detected_c2(rows), eps), 'o-', ...
            'Color', colors(idx, :), 'LineWidth', 1.5);
    end
    set(gca, 'YScale', 'log');
    grid on;
    xlabel('SNR (dB)');
    ylabel('BER');
    title('estimated CSI and detected c2');

    nexttile;
    hold on;
    for idx = 1:numel(fractionalValues)
        rows = abs(results.fractional_doppler - fractionalValues(idx)) < 1e-12;
        semilogy(results.SNR_dB(rows), max(results.BER_perfect_CSI(rows), eps), '-', ...
            'Color', colors(idx, :), 'LineWidth', 1.5);
        semilogy(results.SNR_dB(rows), max(results.BER_estimated_CSI_known_c2(rows), eps), '--', ...
            'Color', colors(idx, :), 'LineWidth', 1.5);
    end
    set(gca, 'YScale', 'log');
    grid on;
    xlabel('SNR (dB)');
    ylabel('BER');
    title('perfect vs estimated CSI');

    legendLabels = compose('fractional=%.2f', fractionalValues);
    legend(nexttile(1), legendLabels, 'Location', 'southeast');
    exportgraphics(fig, pngPath, 'Resolution', 200);
    close(fig);
end

function write_log(logPath, c2Codebook, pilotC2, pilotSnrOffsetDb)
    fid = fopen(logPath, 'w');
    if fid < 0
        error('Unable to open log file: %s', logPath);
    end
    cleanup = onCleanup(@() fclose(fid));

    fprintf(fid, 'Common-pilot c2 reconstruction validation\n');
    fprintf(fid, 'c2 codebook size: %d x %d\n', size(c2Codebook, 1), size(c2Codebook, 2));
    fprintf(fid, 'pilot c2: %.12g\n', pilotC2);
    fprintf(fid, 'pilot SNR offset: %.2f dB\n', pilotSnrOffsetDb);
    fprintf(fid, 'Numerical summary is stored in the CSV file.\n');
    fprintf(fid, 'The mode-confusion tensor and c2 codebook are stored in the MAT file.\n');
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = defaultValue;
    end
end
