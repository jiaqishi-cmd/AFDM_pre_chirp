function results = run_mimo_c2_papr_separability_gate(options)
%RUN_MIMO_C2_PAPR_SEPARABILITY_GATE Test whether uncoupled MIMO adds a problem.
%
% Independent antennas and a Cartesian c2 codebook make minimax PAPR
% selection separable. This experiment verifies the identity numerically and
% measures the cost of imposing one common c2 mode across all antennas.

    rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    setup_paths(rootDir);

    if nargin < 1
        options = struct();
    end

    N = get_option(options, 'N', 64);
    antennaCounts = get_option(options, 'antenna_counts', [2, 4, 8]);
    numFrames = get_option(options, 'num_frames', 5000);
    MMod = get_option(options, 'M_mod', 16);
    seed = get_option(options, 'seed', 20260910);
    c2Base = get_option(options, 'c2_base', sqrt(2) / (10 * N));
    delta = get_option(options, 'delta', c2Base / 16);
    c2Codebook = get_option(options, 'c2_codebook', ...
        build_default_codebook(N, c2Base, delta));

    maxAntennas = max(antennaCounts);
    numModes = size(c2Codebook, 2);
    independentPapr = zeros(numFrames, numel(antennaCounts));
    exhaustivePapr = zeros(numFrames, numel(antennaCounts));
    commonModePapr = zeros(numFrames, numel(antennaCounts));
    independentCommonAgreement = zeros(numFrames, numel(antennaCounts));

    n = (0:N-1).';
    phaseCodebook = exp(1i * 2 * pi * c2Codebook .* (n.^2));
    rng(seed, 'twister');

    for frameIdx = 1:numFrames
        [symbols, ~] = afdm.tx.random_data(N * maxAntennas, MMod, 'qam');
        symbols = reshape(symbols, N, maxAntennas);
        paprMatrix = zeros(maxAntennas, numModes);
        for modeIdx = 1:numModes
            waveforms = ifft(symbols .* phaseCodebook(:, modeIdx), [], 1) * sqrt(N);
            power = abs(waveforms).^2;
            paprMatrix(:, modeIdx) = 10 * log10( ...
                max(power, [], 1).' ./ mean(power, 1).');
        end

        for antennaIdx = 1:numel(antennaCounts)
            numAntennas = antennaCounts(antennaIdx);
            activePapr = paprMatrix(1:numAntennas, :);
            perAntennaMinimum = min(activePapr, [], 2);
            independentPapr(frameIdx, antennaIdx) = max(perAntennaMinimum);
            if numModes^numAntennas <= 256
                exhaustivePapr(frameIdx, antennaIdx) = ...
                    exhaustive_cartesian_minimax(activePapr);
            else
                exhaustivePapr(frameIdx, antennaIdx) = ...
                    independentPapr(frameIdx, antennaIdx);
            end
            commonModePapr(frameIdx, antennaIdx) = min(max(activePapr, [], 1));

            [~, independentModes] = min(activePapr, [], 2);
            independentCommonAgreement(frameIdx, antennaIdx) = ...
                isscalar(unique(independentModes));
        end
    end

    rows = repmat(empty_row(), numel(antennaCounts), 1);
    for antennaIdx = 1:numel(antennaCounts)
        difference = commonModePapr(:, antennaIdx) - independentPapr(:, antennaIdx);
        rows(antennaIdx).num_tx_antennas = antennaCounts(antennaIdx);
        rows(antennaIdx).num_modes = numModes;
        rows(antennaIdx).cartesian_combinations = numModes^antennaCounts(antennaIdx);
        rows(antennaIdx).independent_mode_bits = ...
            antennaCounts(antennaIdx) * log2(numModes);
        rows(antennaIdx).common_mode_bits = log2(numModes);
        rows(antennaIdx).mean_independent_minimax_papr = ...
            mean(independentPapr(:, antennaIdx));
        rows(antennaIdx).p99_independent_minimax_papr = ...
            percentile_by_sort(independentPapr(:, antennaIdx), 0.99);
        rows(antennaIdx).mean_common_mode_papr = ...
            mean(commonModePapr(:, antennaIdx));
        rows(antennaIdx).p99_common_mode_papr = ...
            percentile_by_sort(commonModePapr(:, antennaIdx), 0.99);
        rows(antennaIdx).mean_common_mode_penalty = mean(difference);
        rows(antennaIdx).p99_common_mode_penalty = ...
            percentile_by_sort(difference, 0.99);
        rows(antennaIdx).independent_common_agreement = ...
            mean(independentCommonAgreement(:, antennaIdx));
        rows(antennaIdx).exhaustive_identity_verified = ...
            numModes^antennaCounts(antennaIdx) <= 256;
        if rows(antennaIdx).exhaustive_identity_verified
            rows(antennaIdx).max_joint_identity_error = max(abs( ...
                independentPapr(:, antennaIdx) - exhaustivePapr(:, antennaIdx)));
        end
    end
    summary = struct2table(rows);

    results.config = struct( ...
        'N', N, 'antenna_counts', antennaCounts, 'num_frames', numFrames, ...
        'M_mod', MMod, 'seed', seed, 'c2_base', c2Base, 'delta', delta);
    results.c2_codebook = c2Codebook;
    results.summary = summary;

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    save(fullfile(outputDir, 'mimo_c2_papr_separability_gate.mat'), 'results');
    writetable(summary, ...
        fullfile(outputDir, 'mimo_c2_papr_separability_gate.csv'));
    plot_results(summary, outputDir);

    fprintf('MIMO c2 PAPR separability gate\n');
    disp(summary);
end

function c2Codebook = build_default_codebook(N, c2Base, delta)
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

function value = exhaustive_cartesian_minimax(paprMatrix)
    numAntennas = size(paprMatrix, 1);
    numModes = size(paprMatrix, 2);
    value = Inf;

    for combination = 0:numModes^numAntennas-1
        state = combination;
        candidateValue = -Inf;
        for antennaIdx = 1:numAntennas
            modeIdx = mod(state, numModes) + 1;
            state = floor(state / numModes);
            candidateValue = max(candidateValue, paprMatrix(antennaIdx, modeIdx));
            if candidateValue >= value
                break;
            end
        end
        value = min(value, candidateValue);
    end
end

function plot_results(summary, outputDir)
    fig = figure('Name', 'MIMO PAPR separability gate', 'Color', 'w');
    tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plot(summary.num_tx_antennas, summary.mean_independent_minimax_papr, ...
        'o-', 'LineWidth', 1.6);
    hold on;
    plot(summary.num_tx_antennas, summary.mean_common_mode_papr, ...
        's-', 'LineWidth', 1.6);
    grid on;
    xlabel('number of transmit antennas');
    ylabel('mean worst-antenna PAPR (dB)');
    legend({'independent modes', 'one common mode'}, 'Location', 'best');
    title('MIMO mode constraint');

    nexttile;
    yyaxis left;
    plot(summary.num_tx_antennas, summary.mean_common_mode_penalty, ...
        'o-', 'LineWidth', 1.6);
    ylabel('mean common-mode penalty (dB)');
    yyaxis right;
    plot(summary.num_tx_antennas, summary.independent_common_agreement, ...
        's-', 'LineWidth', 1.6);
    ylabel('independent selections all equal');
    ylim([0, 1.02]);
    grid on;
    xlabel('number of transmit antennas');
    title('Coupling cost');

    exportgraphics(fig, ...
        fullfile(outputDir, 'mimo_c2_papr_separability_gate.png'), ...
        'Resolution', 200);
    close(fig);
end

function value = percentile_by_sort(samples, probability)
    samples = sort(samples(:));
    index = max(1, min(numel(samples), ceil(probability * numel(samples))));
    value = samples(index);
end

function row = empty_row()
    row = struct( ...
        'num_tx_antennas', NaN, ...
        'num_modes', NaN, ...
        'cartesian_combinations', NaN, ...
        'independent_mode_bits', NaN, ...
        'common_mode_bits', NaN, ...
        'mean_independent_minimax_papr', NaN, ...
        'p99_independent_minimax_papr', NaN, ...
        'mean_common_mode_papr', NaN, ...
        'p99_common_mode_papr', NaN, ...
        'mean_common_mode_penalty', NaN, ...
        'p99_common_mode_penalty', NaN, ...
        'independent_common_agreement', NaN, ...
        'exhaustive_identity_verified', false, ...
        'max_joint_identity_error', NaN);
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = defaultValue;
    end
end
