function results = run_explicit_c2_side_info_baseline(numFrames, snrValues, options)
%RUN_EXPLICIT_C2_SIDE_INFO_BASELINE Explicit mode-index transmission baseline.
%
% The side-information channel represents a fixed, already equalized control
% field. It uses BPSK and configurable repetition. The AFDM data channel uses
% exact physical CSI so the experiment isolates side-information failures.

    rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    setup_paths(rootDir);

    if nargin < 1 || isempty(numFrames)
        numFrames = 1000;
    end
    if nargin < 2 || isempty(snrValues)
        snrValues = [10, 15, 20];
    end
    if nargin < 3
        options = struct();
    end

    cfg = configure_base_config(options);
    seed = get_option(options, 'seed', cfg.simulation.random_seed);
    sideInfoSnrOffset = get_option(options, 'side_info_snr_offset_db', -10);
    repetitions = get_option(options, 'repetitions', [1, 3]);
    c2Codebook = get_option(options, 'c2_codebook', build_default_codebook(cfg));
    numModes = size(c2Codebook, 2);
    numModeBits = ceil(log2(numModes));
    methods = ["ideal", compose("BPSK-R%d", repetitions)];
    effectiveBasis = precompute_effective_basis(cfg, c2Codebook);

    numRows = numel(snrValues) * numel(methods);
    rows = repmat(empty_row(), numRows, 1);
    rowIdx = 0;

    fprintf('========== explicit c2 side-information baseline ==========\n');
    fprintf('Frames/SNR=%d, N=%d, M=%d, modes=%d (%d bits)\n', ...
        numFrames, cfg.waveform.NumSubcarriers, cfg.modulation.M_mod, ...
        numModes, numModeBits);

    oldRngState = rng;
    cleanup = onCleanup(@() rng(oldRngState));

    for snrIdx = 1:numel(snrValues)
        dataSnrDb = snrValues(snrIdx);
        sideInfoSnrDb = dataSnrDb + sideInfoSnrOffset;
        accumulators = repmat(initialize_accumulator(), numel(methods), 1);

        for frameIdx = 1:numFrames
            rng(seed + 100000 * snrIdx + frameIdx, 'twister');
            frame = simulate_frame( ...
                cfg, c2Codebook, effectiveBasis, ...
                dataSnrDb, sideInfoSnrDb, repetitions);

            for methodIdx = 1:numel(methods)
                accumulators(methodIdx).mode_errors = ...
                    accumulators(methodIdx).mode_errors + ...
                    double(frame(methodIdx).rx_mode ~= frame(methodIdx).tx_mode);
                accumulators(methodIdx).si_bit_errors = ...
                    accumulators(methodIdx).si_bit_errors + frame(methodIdx).si_bit_errors;
                accumulators(methodIdx).si_bits = ...
                    accumulators(methodIdx).si_bits + frame(methodIdx).si_bits;
                accumulators(methodIdx).data_bit_errors = ...
                    accumulators(methodIdx).data_bit_errors + frame(methodIdx).data_bit_errors;
                accumulators(methodIdx).data_bits = ...
                    accumulators(methodIdx).data_bits + frame(methodIdx).data_bits;
                accumulators(methodIdx).block_errors = ...
                    accumulators(methodIdx).block_errors + frame(methodIdx).block_error;
                accumulators(methodIdx).header_symbols = ...
                    accumulators(methodIdx).header_symbols + frame(methodIdx).header_symbols;
            end
        end

        for methodIdx = 1:numel(methods)
            rowIdx = rowIdx + 1;
            accum = accumulators(methodIdx);
            row = empty_row();
            row.method = methods(methodIdx);
            row.data_SNR_dB = dataSnrDb;
            row.side_info_SNR_dB = sideInfoSnrDb;
            row.num_frames = numFrames;
            row.mode_error_rate = accum.mode_errors / numFrames;
            if accum.si_bits > 0
                row.side_info_BER = accum.si_bit_errors / accum.si_bits;
            else
                row.side_info_BER = 0;
            end
            row.data_BER = accum.data_bit_errors / accum.data_bits;
            row.data_BLER = accum.block_errors / numFrames;
            row.header_symbols_per_block = accum.header_symbols / numFrames;
            row.header_overhead_fraction = row.header_symbols_per_block / ...
                (cfg.waveform.NumSubcarriers * log2(cfg.modulation.M_mod));
            rows(rowIdx) = row;
        end
    end

    results = struct2table(rows);
    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    save(fullfile(outputDir, 'explicit_c2_side_info_baseline.mat'), ...
        'results', 'cfg', 'c2Codebook', 'repetitions', 'sideInfoSnrOffset');
    writetable(results, ...
        fullfile(outputDir, 'explicit_c2_side_info_baseline.csv'));
    plot_results(results, methods, outputDir);

    disp(results);
end

function cfg = configure_base_config(options)
    cfg = afdm_config();
    N = get_option(options, 'N', 64);
    cfg.waveform.NumSubcarriers = N;
    cfg.waveform.CPPLength = get_option(options, 'cpp_length', 3);
    cfg.waveform.c1 = 7 / (2 * N);
    cfg.modulation.M_mod = get_option(options, 'M_mod', 16);
    cfg.modulation.modType = 'qam';
    cfg.channel.delay_taps = get_option(options, 'delay_taps', [0, 1, 2]);
    cfg.channel.doppler_taps = get_option(options, ...
        'doppler_taps', [-2.2, 0.35, 2.4]);
    cfg.channel.doppler_freq = cfg.channel.doppler_taps / N;
    cfg.channel.multipath = true;
    cfg.channel.add_noise = true;
    cfg.pre_chirp.base_c2 = sqrt(2) / (10 * N);
end

function c2Codebook = build_default_codebook(cfg)
    N = cfg.waveform.NumSubcarriers;
    c2Base = cfg.pre_chirp.base_c2;
    delta = c2Base / 16;
    patterns = [ ...
        1, 1, 1, 1; ...
        2, 3, 2, 3; ...
        3, 2, 3, 2; ...
        2, 2, 3, 3];
    c2Codebook = zeros(N, size(patterns, 1));
    for modeIdx = 1:size(patterns, 1)
        c2Codebook(:, modeIdx) = afdm.chirp.build_proposed_pattern( ...
            N, size(patterns, 2), patterns(modeIdx, :), c2Base, delta);
    end
end

function basis = precompute_effective_basis(cfg, c2Codebook)
    numModes = size(c2Codebook, 2);
    numPaths = numel(cfg.channel.delay_taps);
    basis = cell(numModes, numPaths);
    for modeIdx = 1:numModes
        for pathIdx = 1:numPaths
            gains = zeros(1, numPaths);
            gains(pathIdx) = 1;
            basis{modeIdx, pathIdx} = afdm.rx.estimate_effective_channel( ...
                cfg.waveform.NumSubcarriers, cfg.waveform.c1, ...
                c2Codebook(:, modeIdx), gains, cfg.channel.delay_taps, ...
                cfg.channel.doppler_freq);
        end
    end
end

function frame = simulate_frame( ...
        cfg, c2Codebook, effectiveBasis, dataSnrDb, sideInfoSnrDb, repetitions)
    N = cfg.waveform.NumSubcarriers;
    numModes = size(c2Codebook, 2);
    gains = (randn(1, numel(cfg.channel.delay_taps)) + ...
        1i * randn(1, numel(cfg.channel.delay_taps))) / sqrt(2);
    gains = gains / norm(gains);
    cfgFrame = cfg;
    cfgFrame.channel.chan_coef = gains;

    [symbols, txBits] = afdm.tx.random_data( ...
        N, cfg.modulation.M_mod, cfg.modulation.modType);
    [signalCpp, txMode] = select_low_papr_signal(symbols, cfgFrame, c2Codebook);
    cleanReceived = afdm.channel.multipath(signalCpp, cfgFrame);
    [dataReceived, noiseVariance] = add_noise_at_snr(cleanReceived, dataSnrDb);

    rxModes = zeros(1 + numel(repetitions), 1);
    siBitErrors = zeros(size(rxModes));
    siBits = zeros(size(rxModes));
    headerSymbols = zeros(size(rxModes));
    rxModes(1) = txMode;

    for repetitionIdx = 1:numel(repetitions)
        repetition = repetitions(repetitionIdx);
        encoded = afdm.tx.encode_mode_index(txMode, numModes, repetition);
        [receivedHeader, ~] = add_noise_at_snr(encoded.symbols, sideInfoSnrDb);
        decoded = afdm.rx.decode_mode_index(receivedHeader, numModes, repetition);
        methodIdx = repetitionIdx + 1;
        rxModes(methodIdx) = decoded.mode_index;
        siBitErrors(methodIdx) = sum(decoded.bits ~= encoded.bits);
        siBits(methodIdx) = encoded.num_information_bits;
        headerSymbols(methodIdx) = numel(encoded.symbols);
    end

    frame = repmat(struct( ...
        'tx_mode', txMode, 'rx_mode', NaN, 'si_bit_errors', 0, ...
        'si_bits', 0, 'data_bit_errors', 0, 'data_bits', numel(txBits), ...
        'block_error', 0, 'header_symbols', 0), numel(rxModes), 1);

    effectiveChannels = cell(numModes, 1);
    for modeIdx = 1:numModes
        effectiveChannels{modeIdx} = combine_effective_basis( ...
            effectiveBasis(modeIdx, :), gains);
    end

    for methodIdx = 1:numel(rxModes)
        rxMode = rxModes(methodIdx);
        rxBits = receive_data(dataReceived, noiseVariance, cfgFrame, ...
            c2Codebook(:, rxMode), effectiveChannels{rxMode});
        bitErrors = sum(rxBits(:) ~= txBits(:));
        frame(methodIdx).rx_mode = rxMode;
        frame(methodIdx).si_bit_errors = siBitErrors(methodIdx);
        frame(methodIdx).si_bits = siBits(methodIdx);
        frame(methodIdx).data_bit_errors = bitErrors;
        frame(methodIdx).block_error = double(bitErrors > 0);
        frame(methodIdx).header_symbols = headerSymbols(methodIdx);
    end
end

function hEffective = combine_effective_basis(pathBasis, gains)
    hEffective = zeros(size(pathBasis{1}));
    for pathIdx = 1:numel(gains)
        hEffective = hEffective + gains(pathIdx) * pathBasis{pathIdx};
    end
end

function [signalCpp, selectedMode] = select_low_papr_signal(symbols, cfg, c2Codebook)
    numModes = size(c2Codebook, 2);
    papr = zeros(numModes, 1);
    signals = cell(numModes, 1);
    for modeIdx = 1:numModes
        signal = afdm.tx.idaft_mod( ...
            symbols, cfg.waveform.NumSubcarriers, cfg.waveform.c1, ...
            c2Codebook(:, modeIdx));
        papr(modeIdx) = afdm.tx.compute_papr(signal);
        signals{modeIdx} = afdm.tx.add_cpp( ...
            signal, cfg.waveform.CPPLength, cfg.waveform.c1);
    end
    [~, selectedMode] = min(papr);
    signalCpp = signals{selectedMode};
end

function rxBits = receive_data(received, noiseVariance, cfg, c2, hEffective)
    data = afdm.rx.remove_cpp(received, cfg.waveform.CPPLength);
    yDaft = afdm.rx.daft_demod( ...
        data, cfg.waveform.NumSubcarriers, cfg.waveform.c1, c2);
    estimated = afdm.rx.mmse_equalize(yDaft, hEffective, noiseVariance);
    rxBits = afdm.rx.symbol_decision( ...
        estimated, cfg.modulation.M_mod, cfg.modulation.modType);
end

function [noisy, noiseVariance] = add_noise_at_snr(clean, snrDb)
    noiseVariance = mean(abs(clean).^2) / (10^(snrDb / 10));
    noise = sqrt(noiseVariance / 2) * ...
        (randn(size(clean)) + 1i * randn(size(clean)));
    noisy = clean + noise;
end

function accum = initialize_accumulator()
    accum = struct( ...
        'mode_errors', 0, ...
        'si_bit_errors', 0, ...
        'si_bits', 0, ...
        'data_bit_errors', 0, ...
        'data_bits', 0, ...
        'block_errors', 0, ...
        'header_symbols', 0);
end

function plot_results(results, methods, outputDir)
    fig = figure('Name', 'Explicit c2 side information', 'Color', 'w');
    tiledlayout(fig, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    colors = lines(numel(methods));

    nexttile;
    hold on;
    for methodIdx = 1:numel(methods)
        rows = results.method == methods(methodIdx);
        semilogy(results.data_SNR_dB(rows), ...
            max(results.mode_error_rate(rows), eps), 'o-', ...
            'Color', colors(methodIdx, :), 'LineWidth', 1.5);
    end
    set(gca, 'YScale', 'log');
    grid on;
    xlabel('data SNR (dB)');
    ylabel('mode error rate');
    title('Explicit mode recovery');

    nexttile;
    hold on;
    for methodIdx = 1:numel(methods)
        rows = results.method == methods(methodIdx);
        semilogy(results.data_SNR_dB(rows), max(results.data_BER(rows), eps), ...
            'o-', 'Color', colors(methodIdx, :), 'LineWidth', 1.5);
    end
    set(gca, 'YScale', 'log');
    grid on;
    xlabel('data SNR (dB)');
    ylabel('data BER');
    title('End-to-end BER');

    nexttile;
    hold on;
    for methodIdx = 1:numel(methods)
        rows = results.method == methods(methodIdx);
        plot(results.data_SNR_dB(rows), results.data_BLER(rows), ...
            'o-', 'Color', colors(methodIdx, :), 'LineWidth', 1.5);
    end
    grid on;
    xlabel('data SNR (dB)');
    ylabel('data BLER');
    title('Block-error propagation');

    legend(nexttile(1), methods, 'Location', 'best');
    exportgraphics(fig, ...
        fullfile(outputDir, 'explicit_c2_side_info_baseline.png'), ...
        'Resolution', 200);
    close(fig);
end

function row = empty_row()
    row = struct( ...
        'method', "", ...
        'data_SNR_dB', NaN, ...
        'side_info_SNR_dB', NaN, ...
        'num_frames', NaN, ...
        'mode_error_rate', NaN, ...
        'side_info_BER', NaN, ...
        'data_BER', NaN, ...
        'data_BLER', NaN, ...
        'header_symbols_per_block', NaN, ...
        'header_overhead_fraction', NaN);
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = defaultValue;
    end
end
