function results = run_small_n_ml_diversity_validation(options)
%RUN_SMALL_N_ML_DIVERSITY_VALIDATION Validate diversity with exact ML.

    if nargin < 1
        options = struct();
    end

    rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    setup_paths(rootDir);

    N = get_option(options, 'N', 8);
    V = get_option(options, 'V', 4);
    alphaMax = get_option(options, 'alpha_max', 2);
    c1 = get_option(options, 'c1', (2 * alphaMax + 1) / (2 * N));
    c2Base = get_option(options, 'c2_base', sqrt(2) / (10 * N));
    localDelta = get_option(options, 'proposed_delta', c2Base / 16);
    pattern = get_option(options, 'pattern', [2 2 1 1]);
    delay2 = get_option(options, 'delay2', 2);
    doppler2 = get_option(options, 'doppler2', 2);
    snrDb = get_option(options, 'snr_db', 0:4:32);
    numFrames = get_option(options, 'num_frames', 10000);
    seed = get_option(options, 'seed', 20260908);

    [variants, c2Patterns] = build_variants( ...
        N, V, pattern, c2Base, localDelta);
    numVariants = numel(variants);
    numCodewords = 2^N;
    codewordBits = de2bi(0:numCodewords-1, N, 'left-msb').';
    codewordSymbols = 1i * (1 - 2 * codewordBits);

    h1Matrices = cell(numVariants, 1);
    h2Matrices = cell(numVariants, 1);
    for variantIdx = 1:numVariants
        h1Matrices{variantIdx} = afdm.analysis.build_path_matrix( ...
            N, c1, c2Patterns{variantIdx}, 0, 0);
        h2Matrices{variantIdx} = afdm.analysis.build_path_matrix( ...
            N, c1, c2Patterns{variantIdx}, delay2, doppler2);
    end

    errorBits = zeros(numel(snrDb), numVariants);
    totalBits = numFrames * N * ones(numel(snrDb), numVariants);
    oldRngState = rng;
    cleanup = onCleanup(@() rng(oldRngState)); %#ok<NASGU>

    for snrIdx = 1:numel(snrDb)
        rng(seed + 100000 * snrIdx, 'twister');
        noiseStd = sqrt(10^(-snrDb(snrIdx) / 10) / 2);
        txIndices = randi(numCodewords, numFrames, 1);
        gains = (randn(2, numFrames) + 1i * randn(2, numFrames)) / 2;
        unitNoise = randn(N, numFrames) + 1i * randn(N, numFrames);

        for frameIdx = 1:numFrames
            txIndex = txIndices(frameIdx);
            txBits = codewordBits(:, txIndex);
            txSymbols = codewordSymbols(:, txIndex);
            noise = noiseStd * unitNoise(:, frameIdx);

            for variantIdx = 1:numVariants
                hEff = gains(1, frameIdx) * h1Matrices{variantIdx} + ...
                    gains(2, frameIdx) * h2Matrices{variantIdx};
                received = hEff * txSymbols + noise;
                candidateRx = hEff * codewordSymbols;
                distances = sum(abs(candidateRx - received).^2, 1);
                [~, detectedIndex] = min(distances);
                errorBits(snrIdx, variantIdx) = ...
                    errorBits(snrIdx, variantIdx) + ...
                    sum(codewordBits(:, detectedIndex) ~= txBits);
            end
        end
        fprintf('ML SNR %g dB complete\n', snrDb(snrIdx));
    end

    ber = errorBits ./ totalBits;
    berPlot = max(ber, 0.5 ./ totalBits);
    diversityEstimate = estimate_slopes(snrDb, ber, errorBits);

    results.config = struct( ...
        'N', N, 'V', V, 'alpha_max', alphaMax, 'c1', c1, ...
        'pattern', pattern, 'delay2', delay2, 'doppler2', doppler2, ...
        'num_frames', numFrames, 'seed', seed);
    results.variants = variants;
    results.snr_db = snrDb;
    results.ber = ber;
    results.ber_plot = berPlot;
    results.error_bits = errorBits;
    results.total_bits = totalBits;
    results.diversity_estimate = diversityEstimate;

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    save(fullfile(outputDir, 'small_n_ml_diversity_validation.mat'), 'results');
    write_summary(results, outputDir);
    plot_results(results, outputDir);
end

function [variants, c2Patterns] = build_variants( ...
        N, V, pattern, c2Base, localDelta)
    variants = ["baseline", "GPS-rational", "GPS-irrational-k1", ...
        "GPS-irrational-k2", "proposed-2", "proposed-3"];
    c2Patterns = cell(numel(variants), 1);
    c2Patterns{1} = c2Base;
    c2Patterns{2} = afdm.chirp.build_gps_pattern_variant( ...
        N, V, pattern, 'rational');
    c2Patterns{3} = afdm.chirp.build_gps_pattern_variant( ...
        N, V, pattern, 'irrational', struct('precision_k', 1));
    c2Patterns{4} = afdm.chirp.build_gps_pattern_variant( ...
        N, V, pattern, 'irrational', struct('precision_k', 2));
    groupIndex = afdm.chirp.group_index(N, V, 'contiguous');
    candidate2 = c2Base + repmat([-localDelta, localDelta], N, 1);
    candidate3 = c2Base + repmat([0, -localDelta, localDelta], N, 1);
    c2Patterns{5} = select_pattern(candidate2, groupIndex, pattern);
    c2Patterns{6} = select_pattern(candidate3, groupIndex, pattern);
end

function c2 = select_pattern(candidateSet, groupIndex, pattern)
    c2 = zeros(size(groupIndex));
    for groupIdx = 1:numel(pattern)
        selected = groupIndex == groupIdx;
        c2(selected) = candidateSet(selected, pattern(groupIdx));
    end
end

function slopes = estimate_slopes(snrDb, ber, errorBits)
    snrLinear = 10.^(snrDb(:) / 10);
    slopes = NaN(1, size(ber, 2));
    for variantIdx = 1:size(ber, 2)
        valid = snrDb(:) >= 16 & ber(:, variantIdx) > 0 & ...
            errorBits(:, variantIdx) >= 20;
        if nnz(valid) >= 3
            fit = polyfit(log10(snrLinear(valid)), ...
                log10(ber(valid, variantIdx)), 1);
            slopes(variantIdx) = -fit(1);
        end
    end
end

function write_summary(results, outputDir)
    summary = table(results.variants(:), results.diversity_estimate(:), ...
        'VariableNames', {'variant', 'estimated_diversity'});
    for snrIdx = 1:numel(results.snr_db)
        name = sprintf('ber_%ddB', round(results.snr_db(snrIdx)));
        summary.(name) = results.ber(snrIdx, :).';
    end
    writetable(summary, ...
        fullfile(outputDir, 'small_n_ml_diversity_validation.csv'));
    disp(summary);
end

function plot_results(results, outputDir)
    figure('Name', 'Small-N exact ML diversity', 'Color', 'w');
    semilogy(results.snr_db, results.ber_plot, 'LineWidth', 1.8);
    grid on;
    xlabel('SNR (dB)');
    ylabel('BER');
    title(sprintf('Exact ML BER, N=%d, random two-path Rayleigh gains', ...
        results.config.N));
    legend(cellstr(results.variants), 'Location', 'southwest');
    saveas(gcf, fullfile(outputDir, 'small_n_ml_diversity_validation.png'));
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = defaultValue;
    end
end
