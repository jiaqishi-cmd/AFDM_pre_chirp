function results = run_semiblind_c2_candidate_detection(numFrames, snrValues, options)
%RUN_SEMIBLIND_C2_CANDIDATE_DETECTION Semi-blind c2 candidate detection.
%
% Purpose:
% The receiver does not know the c2 index selected by the proposed pre-chirp
% transmitter, but it knows the candidate set. For each candidate c2, the
% receiver runs DAFT demodulation and equalization, then selects the c2 index
% that minimizes the constellation residual
%   metric(k) = sum(abs(x_hat - slicer(x_hat)).^2).
% This experiment supports a practical receiver-side candidate detection
% scheme when explicit pattern signaling is not used.

    rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    setup_paths(rootDir);

    if nargin < 1 || isempty(numFrames)
        numFrames = 500;
    end
    if nargin < 2 || isempty(snrValues)
        snrValues = [15, 20, 25, 30];
    end
    if nargin < 3
        options = struct();
    end

    cfg = configure_base_config(options);
    seed = get_option(options, 'seed', cfg.simulation.random_seed);
    c2Base = cfg.pre_chirp.base_c2;
    c2Candidates = get_option(options, 'c2_candidates', c2Base * [0.8, 1.0, 1.2]);
    c2Candidates = c2Candidates(:).';
    numCandidates = numel(c2Candidates);
    numSnr = numel(snrValues);

    indexCorrect = zeros(numSnr, 1);
    errOracle = zeros(numSnr, 1);
    errDetected = zeros(numSnr, 1);
    errBaseline = zeros(numSnr, 1);
    bitCount = zeros(numSnr, 1);
    selectedCounts = zeros(numSnr, numCandidates);
    detectedCounts = zeros(numSnr, numCandidates);

    fprintf('========== semi-blind c2 candidate detection ==========\n');
    fprintf('Frames/SNR=%d, N=%d, M=%d\n', ...
        numFrames, cfg.waveform.NumSubcarriers, cfg.modulation.M_mod);
    fprintf('c2 candidates: %s\n', mat2str(c2Candidates, 6));

    oldRngState = rng;
    cleanup = onCleanup(@() rng(oldRngState)); %#ok<NASGU>

    for snrIdx = 1:numSnr
        cfg.channel.snr_db = snrValues(snrIdx);
        fprintf('\nSNR = %.1f dB\n', cfg.channel.snr_db);

        for frameIdx = 1:numFrames
            rng(seed + 100000 * snrIdx + frameIdx, 'twister');
            frame = simulate_detection_frame(cfg, c2Candidates);

            indexCorrect(snrIdx) = indexCorrect(snrIdx) + double(frame.detected_index == frame.tx_index);
            selectedCounts(snrIdx, frame.tx_index) = selectedCounts(snrIdx, frame.tx_index) + 1;
            detectedCounts(snrIdx, frame.detected_index) = detectedCounts(snrIdx, frame.detected_index) + 1;

            errOracle(snrIdx) = errOracle(snrIdx) + frame.oracle.err_bits;
            errDetected(snrIdx) = errDetected(snrIdx) + frame.detected.err_bits;
            errBaseline(snrIdx) = errBaseline(snrIdx) + frame.baseline.err_bits;
            bitCount(snrIdx) = bitCount(snrIdx) + frame.oracle.total_bits;

            if mod(frameIdx, 100) == 0 || frameIdx == numFrames
                fprintf('  Frame %5d/%d complete\n', frameIdx, numFrames);
            end
        end

        fprintf('  selected index counts: %s\n', mat2str(selectedCounts(snrIdx, :)));
        fprintf('  detected index counts: %s\n', mat2str(detectedCounts(snrIdx, :)));
        fprintf('  index accuracy: %.4f\n', indexCorrect(snrIdx) / numFrames);
    end

    results = table( ...
        snrValues(:), ...
        repmat(numFrames, numSnr, 1), ...
        indexCorrect ./ numFrames, ...
        errOracle ./ bitCount, ...
        errDetected ./ bitCount, ...
        errBaseline ./ bitCount, ...
        'VariableNames', {'SNR_dB', 'num_frames', 'index_accuracy', ...
        'BER_oracle', 'BER_detected', 'BER_baseline'});

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    csvPath = fullfile(outputDir, ['semiblind_c2_detection_' timestamp '.csv']);
    pngPath = fullfile(outputDir, ['semiblind_c2_detection_' timestamp '.png']);
    matPath = fullfile(outputDir, ['semiblind_c2_detection_' timestamp '.mat']);
    logPath = fullfile(outputDir, ['semiblind_c2_detection_' timestamp '.log']);

    writetable(results, csvPath);
    save(matPath, 'results', 'selectedCounts', 'detectedCounts', 'c2Candidates', 'cfg');
    write_log(logPath, results, selectedCounts, detectedCounts, c2Candidates);
    plot_results(results, pngPath);

    fprintf('\nSaved semi-blind c2 detection results to %s\n', csvPath);
    fprintf('Saved selected-index log to %s\n', logPath);
    disp(results);
end

function cfg = configure_base_config(options)
    cfg = afdm_config();
    N = get_option(options, 'N', 64);
    MMod = get_option(options, 'M_mod', 16);

    cfg.waveform.NumSubcarriers = N;
    cfg.waveform.CPPLength = get_option(options, 'cpp_length', 3);
    cfg.waveform.c1 = 7 / (2 * N);
    cfg.modulation.M_mod = MMod;
    cfg.modulation.modType = get_option(options, 'modType', 'qam');
    cfg.channel.add_noise = true;
    cfg.channel.multipath = get_option(options, 'multipath', true);
    cfg.receiver.detector = get_option(options, 'detector', 'full_mmse');

    c2Base = sqrt(2) / (10 * N);
    cfg.pre_chirp.scheme = 'baseline';
    cfg.pre_chirp.base_c2 = c2Base;
    cfg.pre_chirp.num_groups = 1;
    cfg.pre_chirp.num_candidates = 3;
    cfg.pre_chirp.delta = c2Base / 5;
    cfg.waveform.c2 = c2Base;
end

function frame = simulate_detection_frame(cfg, c2Candidates)
    N = cfg.waveform.NumSubcarriers;
    MMod = cfg.modulation.M_mod;
    modType = cfg.modulation.modType;

    [symbols, txBits] = afdm.tx.random_data(N, MMod, modType);
    [signalCpp, paprValues, txIdx] = select_low_papr_signal(symbols, cfg, c2Candidates);

    rSignal = afdm.channel.multipath(signalCpp, cfg);
    [rSignal, noiseVar] = afdm.channel.add_awgn(rSignal, cfg);
    cfg.channel.noise_var = noiseVar;

    numCandidates = numel(c2Candidates);
    metrics = zeros(numCandidates, 1);
    candidateBits = cell(numCandidates, 1);

    for candidateIdx = 1:numCandidates
        rx = receive_with_c2(rSignal, cfg, c2Candidates(candidateIdx));
        metrics(candidateIdx) = rx.metric;
        candidateBits{candidateIdx} = rx.bits;
    end

    [~, detectedIdx] = min(metrics);
    baselineIdx = find_baseline_index(c2Candidates, cfg.pre_chirp.base_c2);

    frame.tx_index = txIdx;
    frame.detected_index = detectedIdx;
    frame.baseline_index = baselineIdx;
    frame.papr_values = paprValues;
    frame.metrics = metrics;

    frame.oracle = count_errors(candidateBits{txIdx}, txBits);
    frame.detected = count_errors(candidateBits{detectedIdx}, txBits);
    frame.baseline = count_errors(candidateBits{baselineIdx}, txBits);
end

function [signalCpp, paprValues, selectedIdx] = select_low_papr_signal(symbols, cfg, c2Candidates)
    numCandidates = numel(c2Candidates);
    signals = cell(numCandidates, 1);
    paprValues = zeros(numCandidates, 1);

    for candidateIdx = 1:numCandidates
        signal = afdm.tx.idaft_mod( ...
            symbols, ...
            cfg.waveform.NumSubcarriers, ...
            cfg.waveform.c1, ...
            c2Candidates(candidateIdx));
        signals{candidateIdx} = afdm.tx.add_cpp(signal, cfg.waveform.CPPLength, cfg.waveform.c1);
        paprValues(candidateIdx) = afdm.tx.compute_papr(signal);
    end

    [~, selectedIdx] = min(paprValues);
    signalCpp = signals{selectedIdx};
end

function rx = receive_with_c2(rSignal, cfg, c2Rx)
    rData = afdm.rx.remove_cpp(rSignal, cfg.waveform.CPPLength);
    yDaft = afdm.rx.daft_demod( ...
        rData, ...
        cfg.waveform.NumSubcarriers, ...
        cfg.waveform.c1, ...
        c2Rx);

    if isfield(cfg.channel, 'multipath') && ~cfg.channel.multipath
        hEff = eye(cfg.waveform.NumSubcarriers);
    else
        hEff = afdm.rx.estimate_effective_channel( ...
            cfg.waveform.NumSubcarriers, ...
            cfg.waveform.c1, ...
            c2Rx, ...
            cfg.channel.chan_coef, ...
            cfg.channel.delay_taps, ...
            cfg.channel.doppler_freq);
    end

    xHat = afdm.rx.equalize_symbols(yDaft, hEff, cfg.channel.noise_var, cfg);
    xSliced = slice_symbols(xHat, cfg.modulation.M_mod, cfg.modulation.modType);
    rx.metric = sum(abs(xHat(:) - xSliced(:)).^2);
    rx.bits = afdm.rx.symbol_decision(xHat, cfg.modulation.M_mod, cfg.modulation.modType);
    rx.x_hat = xHat;
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

function out = count_errors(rxBits, txBits)
    [out.err_bits, out.total_bits] = afdm.rx.compute_bit_errors(rxBits, txBits);
end

function baselineIdx = find_baseline_index(c2Candidates, c2Base)
    [~, baselineIdx] = min(abs(c2Candidates - c2Base));
end

function plot_results(results, pngPath)
    fig = figure('Name', 'semi-blind c2 candidate detection', 'Color', 'w');
    tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plot(results.SNR_dB, results.index_accuracy, 'o-', 'LineWidth', 2);
    grid on;
    ylim([0, 1.02]);
    xlabel('SNR (dB)');
    ylabel('index accuracy');
    title('c2 index detection');

    nexttile;
    semilogy(results.SNR_dB, max(results.BER_oracle, eps), 'o-', 'LineWidth', 2);
    hold on;
    semilogy(results.SNR_dB, max(results.BER_detected, eps), 's-', 'LineWidth', 2);
    semilogy(results.SNR_dB, max(results.BER_baseline, eps), '^-', 'LineWidth', 2);
    grid on;
    xlabel('SNR (dB)');
    ylabel('BER');
    legend({'oracle c2', 'detected c2', 'baseline c2'}, 'Location', 'southwest');
    title('BER comparison');

    exportgraphics(fig, pngPath, 'Resolution', 200);
    close(fig);
end

function write_log(logPath, results, selectedCounts, detectedCounts, c2Candidates)
    fid = fopen(logPath, 'w');
    if fid < 0
        error('Unable to open log file: %s', logPath);
    end
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, 'Semi-blind c2 candidate detection\n');
    fprintf(fid, 'c2 candidates: %s\n\n', mat2str(c2Candidates, 12));
    for idx = 1:height(results)
        fprintf(fid, 'SNR %.1f dB\n', results.SNR_dB(idx));
        fprintf(fid, '  num_frames: %d\n', results.num_frames(idx));
        fprintf(fid, '  index_accuracy: %.6f\n', results.index_accuracy(idx));
        fprintf(fid, '  selected index distribution: %s\n', mat2str(selectedCounts(idx, :)));
        fprintf(fid, '  detected index distribution: %s\n', mat2str(detectedCounts(idx, :)));
        fprintf(fid, '  BER oracle/detected/baseline: %.6g / %.6g / %.6g\n\n', ...
            results.BER_oracle(idx), results.BER_detected(idx), results.BER_baseline(idx));
    end
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = defaultValue;
    end
end
