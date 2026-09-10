function results = run_relative_fractional_doppler_closure(options)
%RUN_RELATIVE_FRACTIONAL_DOPPLER_CLOSURE Close the fractional-Doppler audit.

    rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    setup_paths(rootDir);

    if nargin < 1
        options = struct();
    end

    N = get_option(options, 'N', 64);
    V = get_option(options, 'V', 4);
    alphaMax = get_option(options, 'alpha_max', 3);
    c1 = get_option(options, 'c1', (2 * alphaMax + 1) / (2 * N));
    c2Base = get_option(options, 'c2_base', sqrt(2) / (10 * N));
    localDelta = get_option(options, 'proposed_delta', c2Base / 16);
    betaValues = get_option(options, 'beta_values', -0.5:0.025:0.5);
    snrValues = get_option(options, 'snr_values', [28, 40]);
    numFrames = get_option(options, 'num_frames', 3000);
    seed = get_option(options, 'seed', 20260910);
    delays = [0, 2];
    integerDopplers = [0, 2];

    schemes = build_schemes(N, V, c2Base, localDelta);
    dangerDelta = integer_danger_event( ...
        N, c1, schemes(2).c2, delays, integerDopplers);

    numBeta = numel(betaValues);
    firstPaths = cell(numBeta, 1);
    secondPaths = cell(numBeta, 1);
    for betaIdx = 1:numBeta
        firstPaths{betaIdx} = afdm.analysis.build_path_matrix( ...
            N, c1, schemes(2).c2, delays(1), ...
            integerDopplers(1) + betaValues(betaIdx));
        secondPaths{betaIdx} = afdm.analysis.build_path_matrix( ...
            N, c1, schemes(2).c2, delays(2), ...
            integerDopplers(2) + betaValues(betaIdx));
    end

    sigmaGrid = zeros(numBeta, numBeta);
    for firstIdx = 1:numBeta
        for secondIdx = 1:numBeta
            psi = afdm.analysis.build_pairwise_path_response( ...
                {firstPaths{firstIdx}, secondPaths{secondIdx}}, dangerDelta);
            metrics = afdm.analysis.pairwise_diversity_metrics(psi, 1e-10);
            sigmaGrid(secondIdx, firstIdx) = metrics.min_singular_value;
        end
    end

    relativeConsistency = relative_doppler_consistency(sigmaGrid);
    berCases = build_ber_cases();
    berTable = run_fixed_channel_ber( ...
        schemes, berCases, snrValues, numFrames, seed, N, c1, delays, ...
        integerDopplers);

    results.config = struct( ...
        'N', N, 'V', V, 'c1', c1, 'c2_base', c2Base, ...
        'local_delta', localDelta, 'beta_values', betaValues, ...
        'snr_values', snrValues, 'num_frames', numFrames, 'seed', seed);
    results.schemes = schemes;
    results.danger_delta = dangerDelta;
    results.sigma_grid = sigmaGrid;
    results.relative_consistency = relativeConsistency;
    results.ber_table = berTable;

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    save(fullfile(outputDir, 'relative_fractional_doppler_closure.mat'), 'results');
    writetable(berTable, ...
        fullfile(outputDir, 'relative_fractional_doppler_ber.csv'));
    plot_closure(betaValues, sigmaGrid, relativeConsistency, outputDir);
    plot_ber(berTable, schemes, snrValues, N, outputDir);

    fprintf('Relative fractional-Doppler closure audit\n');
    fprintf('Maximum sigma spread at equal delta-beta: %.3e\n', ...
        max(relativeConsistency.sigma_spread));
    disp(berTable);
end

function schemes = build_schemes(N, V, c2Base, localDelta)
    pattern = [2, 2, 1, 1];
    gpsRational = afdm.chirp.gps_candidate_set_variant(N, 'rational');
    gpsSafe = afdm.chirp.gps_candidate_set_variant( ...
        N, 'phase_offset', struct('phase_offset', 0.3));

    schemes = repmat(struct('name', "", 'c2', []), 4, 1);
    schemes(1).name = "baseline";
    schemes(1).c2 = c2Base * ones(N, 1);
    schemes(2).name = "GPS-rational";
    schemes(2).c2 = select_grouped_c2(gpsRational, N, V, pattern);
    schemes(3).name = "local-proposed";
    schemes(3).c2 = afdm.chirp.build_proposed_pattern( ...
        N, V, pattern, c2Base, localDelta);
    schemes(4).name = "integer-safe-GPS";
    schemes(4).c2 = select_grouped_c2(gpsSafe, N, V, pattern);
end

function c2 = select_grouped_c2(candidateSet, N, V, pattern)
    groupIndex = afdm.chirp.group_index(N, V, 'contiguous');
    c2 = zeros(N, 1);
    for groupIdx = 1:V
        selected = groupIndex == groupIdx;
        c2(selected) = candidateSet(selected, pattern(groupIdx));
    end
end

function delta = integer_danger_event(N, c1, c2, delays, dopplers)
    firstPath = afdm.analysis.build_path_matrix( ...
        N, c1, c2, delays(1), dopplers(1));
    secondPath = afdm.analysis.build_path_matrix( ...
        N, c1, c2, delays(2), dopplers(2));
    differenceAlphabet = afdm.delta.generate_sparse_difference_set( ...
        1, 2, 'psk', 1).difference_alphabet;
    cycle = afdm.analysis.monomial_cycle_analysis( ...
        firstPath' * secondPath, differenceAlphabet, 1e-12);
    delta = cycle.closest_delta;
    delta = delta / norm(delta);
end

function consistency = relative_doppler_consistency(sigmaGrid)
    gridSize = size(sigmaGrid, 1);
    offsets = -(gridSize - 1):(gridSize - 1);
    meanSigma = zeros(numel(offsets), 1);
    sigmaSpread = zeros(numel(offsets), 1);
    for offsetIdx = 1:numel(offsets)
        values = diag(sigmaGrid, offsets(offsetIdx));
        meanSigma(offsetIdx) = mean(values);
        sigmaSpread(offsetIdx) = max(values) - min(values);
    end
    consistency = table(offsets(:), meanSigma, sigmaSpread, ...
        'VariableNames', {'grid_offset', 'mean_sigma', 'sigma_spread'});
end

function cases = build_ber_cases()
    cases = struct( ...
        'name', {"aligned-0", "aligned-0.25", "aligned-0.50", ...
        "relative-0.05", "relative-0.25"}, ...
        'beta1', {0, 0.25, 0.5, 0, 0}, ...
        'beta2', {0, 0.25, 0.5, 0.05, 0.25});
end

function berTable = run_fixed_channel_ber( ...
        schemes, cases, snrValues, numFrames, seed, N, c1, delays, dopplers)
    rows = repmat(empty_ber_row(), ...
        numel(schemes) * numel(cases) * numel(snrValues), 1);
    rowIdx = 0;
    bitsPerRun = N * numFrames;

    for caseIdx = 1:numel(cases)
        pathDopplers = dopplers + [cases(caseIdx).beta1, cases(caseIdx).beta2];
        for snrIdx = 1:numel(snrValues)
            rng(seed + 10000 * caseIdx + snrIdx, 'twister');
            bits = randi([0, 1], N, numFrames);
            symbols = 1i * (1 - 2 * bits);
            noiseVariance = 10^(-snrValues(snrIdx) / 10);
            normalizedNoise = (randn(N, numFrames) + ...
                1i * randn(N, numFrames)) / sqrt(2);

            for schemeIdx = 1:numel(schemes)
                firstPath = afdm.analysis.build_path_matrix( ...
                    N, c1, schemes(schemeIdx).c2, delays(1), pathDopplers(1));
                secondPath = afdm.analysis.build_path_matrix( ...
                    N, c1, schemes(schemeIdx).c2, delays(2), pathDopplers(2));
                channel = (firstPath - secondPath) / sqrt(2);
                equalizer = (channel' * channel + noiseVariance * eye(N)) \ channel';
                received = channel * symbols + sqrt(noiseVariance) * normalizedNoise;
                estimated = equalizer * received;
                detectedBits = imag(estimated) < 0;

                rowIdx = rowIdx + 1;
                row = empty_ber_row();
                row.scheme = schemes(schemeIdx).name;
                row.case_name = cases(caseIdx).name;
                row.beta1 = cases(caseIdx).beta1;
                row.beta2 = cases(caseIdx).beta2;
                row.relative_beta = cases(caseIdx).beta2 - cases(caseIdx).beta1;
                row.SNR_dB = snrValues(snrIdx);
                row.num_frames = numFrames;
                row.BER = sum(detectedBits(:) ~= bits(:)) / bitsPerRun;
                rows(rowIdx) = row;
            end
        end
    end
    berTable = struct2table(rows);
end

function plot_closure(betaValues, sigmaGrid, consistency, outputDir)
    fig = figure('Name', 'Relative fractional Doppler closure', 'Color', 'w');
    tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    imagesc(betaValues, betaValues, log10(max(sigmaGrid, eps)));
    axis xy square;
    hold on;
    plot(betaValues, betaValues, 'w--', 'LineWidth', 1.2);
    colorbar;
    xlabel('\beta_1');
    ylabel('\beta_2');
    title('log_{10}(\sigma_{min})');

    nexttile;
    plot(consistency.grid_offset, consistency.mean_sigma, 'o-', 'LineWidth', 1.5);
    grid on;
    xlabel('grid offset proportional to \Delta\beta');
    ylabel('mean \sigma_{min}');
    title('Dependence on relative Doppler');

    exportgraphics(fig, ...
        fullfile(outputDir, 'relative_fractional_doppler_closure.png'), ...
        'Resolution', 200);
    close(fig);
end

function plot_ber(berTable, schemes, snrValues, N, outputDir)
    fig = figure('Name', 'Relative fractional Doppler BER', 'Color', 'w');
    tiledlayout(fig, 1, numel(snrValues), ...
        'Padding', 'compact', 'TileSpacing', 'compact');
    caseNames = unique(berTable.case_name, 'stable');
    colors = lines(numel(schemes));

    for snrIdx = 1:numel(snrValues)
        nexttile;
        hold on;
        for schemeIdx = 1:numel(schemes)
            rows = berTable.scheme == schemes(schemeIdx).name & ...
                berTable.SNR_dB == snrValues(snrIdx);
            values = zeros(numel(caseNames), 1);
            for caseIdx = 1:numel(caseNames)
                values(caseIdx) = berTable.BER(rows & ...
                    berTable.case_name == caseNames(caseIdx));
            end
            semilogy(1:numel(caseNames), max(values, 0.5 / ...
                (berTable.num_frames(find(rows, 1)) * N)), ...
                'o-', 'Color', colors(schemeIdx, :), ...
                'LineWidth', 1.4);
        end
        set(gca, 'YScale', 'log', 'XTick', 1:numel(caseNames), ...
            'XTickLabel', caseNames);
        xtickangle(25);
        grid on;
        ylabel('BER');
        title(sprintf('SNR = %.0f dB', snrValues(snrIdx)));
    end

    legend(nexttile(1), string({schemes.name}), 'Location', 'best');
    exportgraphics(fig, ...
        fullfile(outputDir, 'relative_fractional_doppler_ber.png'), ...
        'Resolution', 200);
    close(fig);
end

function row = empty_ber_row()
    row = struct( ...
        'scheme', "", ...
        'case_name', "", ...
        'beta1', NaN, ...
        'beta2', NaN, ...
        'relative_beta', NaN, ...
        'SNR_dB', NaN, ...
        'num_frames', NaN, ...
        'BER', NaN);
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = defaultValue;
    end
end
