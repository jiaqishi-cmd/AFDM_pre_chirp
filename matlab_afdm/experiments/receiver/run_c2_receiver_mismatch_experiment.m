function results = run_c2_receiver_mismatch_experiment(numFrames, snrDb, options)
%RUN_C2_RECEIVER_MISMATCH_EXPERIMENT Compare BER under c2 receiver mismatch.
%
% Purpose:
% This experiment demonstrates that the selected pre-chirp c2 is necessary
% receiver-side demodulation information. The transmitter selects one c2 per
% frame by minimizing PAPR over a small proposed candidate set. The receiver
% is then evaluated with oracle c2, baseline c2, and an intentionally wrong
% candidate c2. A BER gap supports either explicit pattern signaling or
% receiver-side candidate detection for the proposed pre-chirp scheme.

    rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    setup_paths(rootDir);

    if nargin < 1 || isempty(numFrames)
        numFrames = 1000;
    end
    if nargin < 2 || isempty(snrDb)
        snrDb = 25;
    end
    if nargin < 3
        options = struct();
    end

    cfg = configure_base_config(snrDb, options);
    rngSeed = get_option(options, 'seed', cfg.simulation.random_seed);
    c2Base = cfg.pre_chirp.base_c2;
    c2Candidates = get_option(options, 'c2_candidates', c2Base * [0.8, 1.0, 1.2]);
    c2Candidates = c2Candidates(:).';

    modes = {'oracle', 'baseline', 'wrong_candidate'};
    errBits = zeros(numel(modes), 1);
    totalBits = zeros(numel(modes), 1);
    selectedCounts = zeros(numel(c2Candidates), 1);

    fprintf('========== c2 receiver mismatch experiment ==========\n');
    fprintf('Frames=%d, SNR=%.1f dB, N=%d, M=%d\n', ...
        numFrames, snrDb, cfg.waveform.NumSubcarriers, cfg.modulation.M_mod);
    fprintf('c2 candidates: %s\n', mat2str(c2Candidates, 6));

    oldRngState = rng;
    cleanup = onCleanup(@() rng(oldRngState)); %#ok<NASGU>

    for frameIdx = 1:numFrames
        rng(rngSeed + frameIdx, 'twister');
        frame = simulate_mismatch_frame(cfg, c2Candidates);
        selectedCounts(frame.tx_candidate_index) = selectedCounts(frame.tx_candidate_index) + 1;

        for modeIdx = 1:numel(modes)
            errBits(modeIdx) = errBits(modeIdx) + frame.(modes{modeIdx}).err_bits;
            totalBits(modeIdx) = totalBits(modeIdx) + frame.(modes{modeIdx}).total_bits;
        end

        if mod(frameIdx, 100) == 0 || frameIdx == numFrames
            fprintf('Frame %5d/%d complete\n', frameIdx, numFrames);
        end
    end

    ber = errBits ./ totalBits;
    results = table( ...
        string(modes(:)), ...
        repmat(numFrames, numel(modes), 1), ...
        repmat(snrDb, numel(modes), 1), ...
        ber, ...
        errBits, ...
        totalBits, ...
        'VariableNames', {'mode', 'num_frames', 'SNR_dB', 'BER', 'num_bit_errors', 'num_bits'});

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    csvPath = fullfile(outputDir, ['c2_receiver_mismatch_' timestamp '.csv']);
    pngPath = fullfile(outputDir, ['c2_receiver_mismatch_' timestamp '.png']);
    matPath = fullfile(outputDir, ['c2_receiver_mismatch_' timestamp '.mat']);

    writetable(results, csvPath);
    save(matPath, 'results', 'selectedCounts', 'c2Candidates', 'cfg');
    plot_results(results, selectedCounts, c2Candidates, pngPath);

    fprintf('Saved c2 mismatch results to %s\n', csvPath);
    disp(results);
end

function cfg = configure_base_config(snrDb, options)
    cfg = afdm_config();
    N = get_option(options, 'N', 64);
    MMod = get_option(options, 'M_mod', 16);

    cfg.waveform.NumSubcarriers = N;
    cfg.waveform.CPPLength = get_option(options, 'cpp_length', 3);
    cfg.waveform.c1 = 7 / (2 * N);
    cfg.modulation.M_mod = MMod;
    cfg.modulation.modType = get_option(options, 'modType', 'qam');
    cfg.channel.snr_db = snrDb;
    cfg.channel.add_noise = get_option(options, 'add_noise', true);
    cfg.channel.multipath = get_option(options, 'multipath', true);

    c2Base = sqrt(2) / (10 * N);
    cfg.pre_chirp.scheme = 'baseline';
    cfg.pre_chirp.base_c2 = c2Base;
    cfg.pre_chirp.num_groups = 1;
    cfg.pre_chirp.num_candidates = 3;
    cfg.pre_chirp.delta = c2Base / 5;
    cfg.waveform.c2 = c2Base;
end

function frame = simulate_mismatch_frame(cfg, c2Candidates)
    N = cfg.waveform.NumSubcarriers;
    MMod = cfg.modulation.M_mod;
    modType = cfg.modulation.modType;

    [symbols, txBits] = afdm.tx.random_data(N, MMod, modType);
    [signalCpp, paprValues, selectedIdx] = select_low_papr_signal(symbols, cfg, c2Candidates);
    c2Tx = c2Candidates(selectedIdx);

    rSignal = afdm.channel.multipath(signalCpp, cfg);
    [rSignal, noiseVar] = afdm.channel.add_awgn(rSignal, cfg);
    cfg.channel.noise_var = noiseVar;

    wrongIdx = select_wrong_index(selectedIdx, numel(c2Candidates));
    rxC2.oracle = c2Tx;
    rxC2.baseline = cfg.pre_chirp.base_c2;
    rxC2.wrong_candidate = c2Candidates(wrongIdx);

    modes = fieldnames(rxC2);
    for modeIdx = 1:numel(modes)
        mode = modes{modeIdx};
        rxState.c2 = rxC2.(mode);
        [~, errBits, totalBits] = afdm.rx.engine(rSignal, cfg, txBits, rxState);
        frame.(mode).err_bits = errBits;
        frame.(mode).total_bits = totalBits;
        frame.(mode).c2_rx = rxState.c2;
    end

    frame.tx_candidate_index = selectedIdx;
    frame.wrong_candidate_index = wrongIdx;
    frame.c2_tx = c2Tx;
    frame.papr_values = paprValues;
    frame.noise_var = noiseVar;
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

function wrongIdx = select_wrong_index(selectedIdx, numCandidates)
    wrongIdx = mod(selectedIdx, numCandidates) + 1;
end

function plot_results(results, selectedCounts, c2Candidates, pngPath)
    fig = figure('Name', 'c2 receiver mismatch BER', 'Color', 'w');
    tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    bar(categorical(results.mode), max(results.BER, 1 ./ max(results.num_bits, 1)));
    set(gca, 'YScale', 'log');
    grid on;
    ylabel('BER');
    title('Receiver c2 assumption');

    nexttile;
    bar(categorical(compose('%.3g', c2Candidates)), selectedCounts);
    grid on;
    xlabel('selected c2');
    ylabel('frames');
    title('Transmitter selection counts');

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
