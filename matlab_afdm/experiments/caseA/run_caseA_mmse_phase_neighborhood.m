function results = run_caseA_mmse_phase_neighborhood(options)
%RUN_CASEA_MMSE_PHASE_NEIGHBORHOOD Scan finite-SNR BER around theta=pi.
%   All schemes use the same BPSK data and unit-variance noise samples.

    if nargin < 1
        options = struct();
    end

    rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    setup_paths(rootDir);

    N = get_option(options, 'N', 64);
    V = get_option(options, 'V', 4);
    alphaMax = get_option(options, 'alpha_max', 3);
    c1 = get_option(options, 'c1', (2 * alphaMax + 1) / (2 * N));
    c2Base = get_option(options, 'c2_base', sqrt(2) / (10 * N));
    delta = get_option(options, 'proposed_delta', c2Base / 16);
    pattern = get_option(options, 'pattern', [2 2 1 1]);
    thetaList = get_option(options, 'theta_list', linspace(0.9 * pi, 1.1 * pi, 41));
    snrDb = get_option(options, 'snr_db', 28);
    numFrames = get_option(options, 'num_frames', 2000);
    seed = get_option(options, 'seed', 20260908);
    delays = get_option(options, 'delays', [0, 2]);
    dopplers = get_option(options, 'dopplers', [0, 2]);

    c2Gps = afdm.chirp.build_gps_pattern(N, V, pattern);
    c2Prop = afdm.chirp.build_proposed_pattern( ...
        N, V, pattern, c2Base, delta);
    schemes = ["baseline", "GPS", "proposed"];
    c2Values = {c2Base, c2Gps, c2Prop};

    oldRngState = rng;
    cleanup = onCleanup(@() rng(oldRngState)); %#ok<NASGU>
    rng(seed, 'twister');
    bits = randi([0, 1], N, numFrames);
    symbols = 1i * (1 - 2 * bits);
    unitNoise = (randn(N, numFrames) + 1i * randn(N, numFrames)) / sqrt(2);

    pathMatrices = cell(numel(schemes), numel(delays));
    for schemeIdx = 1:numel(schemes)
        for pathIdx = 1:numel(delays)
            pathMatrices{schemeIdx, pathIdx} = afdm.analysis.build_path_matrix( ...
                N, c1, c2Values{schemeIdx}, ...
                delays(pathIdx), dopplers(pathIdx));
        end
    end

    numTheta = numel(thetaList);
    ber = zeros(numTheta, numel(schemes));
    berPlot = zeros(size(ber));
    errorBits = zeros(size(ber));
    p001Margin = zeros(size(ber));
    minSingular = zeros(size(ber));
    snrLinear = 10^(snrDb / 10);

    for thetaIdx = 1:numTheta
        gains = [1, exp(1i * thetaList(thetaIdx))] / sqrt(2);
        for schemeIdx = 1:numel(schemes)
            hEff = gains(1) * pathMatrices{schemeIdx, 1} + ...
                gains(2) * pathMatrices{schemeIdx, 2};
            cleanRx = hEff * symbols;
            signalPower = mean(abs(cleanRx).^2, 'all');
            noiseVar = signalPower / snrLinear;
            received = cleanRx + sqrt(noiseVar) * unitNoise;
            equalizer = (hEff' * hEff + noiseVar * eye(N)) \ hEff';
            estimate = equalizer * received;
            decisions = pskdemod(estimate(:), 2, pi / 2, 'OutputType', 'bit');
            errors = sum(decisions ~= bits(:));
            margins = real(conj(symbols(:)) .* estimate(:));
            singularValues = svd(hEff);

            errorBits(thetaIdx, schemeIdx) = errors;
            ber(thetaIdx, schemeIdx) = errors / numel(bits);
            berPlot(thetaIdx, schemeIdx) = max(ber(thetaIdx, schemeIdx), ...
                0.5 / numel(bits));
            p001Margin(thetaIdx, schemeIdx) = percentile_by_sort(margins, 0.001);
            minSingular(thetaIdx, schemeIdx) = singularValues(end);
        end
    end

    results.config = struct( ...
        'N', N, 'V', V, 'alpha_max', alphaMax, 'c1', c1, ...
        'c2_base', c2Base, 'delta', delta, 'pattern', pattern, ...
        'snr_db', snrDb, 'num_frames', numFrames, 'seed', seed, ...
        'delays', delays, 'dopplers', dopplers);
    results.theta_list = thetaList;
    results.schemes = schemes;
    results.ber = ber;
    results.ber_plot = berPlot;
    results.error_bits = errorBits;
    results.p001_margin = p001Margin;
    results.min_singular = minSingular;

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    save(fullfile(outputDir, 'caseA_mmse_phase_neighborhood.mat'), 'results');
    write_summary(results, outputDir);
    plot_results(results, outputDir);
end

function write_summary(results, outputDir)
    theta_over_pi = results.theta_list(:) / pi;
    ber_baseline = results.ber(:, 1);
    ber_GPS = results.ber(:, 2);
    ber_proposed = results.ber(:, 3);
    errors_baseline = results.error_bits(:, 1);
    errors_GPS = results.error_bits(:, 2);
    errors_proposed = results.error_bits(:, 3);
    summary = table(theta_over_pi, ber_baseline, ber_GPS, ber_proposed, ...
        errors_baseline, errors_GPS, errors_proposed);
    writetable(summary, fullfile(outputDir, 'caseA_mmse_phase_neighborhood.csv'));

    [peakGpsBer, peakIdx] = max(ber_GPS);
    fprintf('Case A MMSE phase neighborhood, SNR %.1f dB\n', ...
        results.config.snr_db);
    fprintf('Peak GPS BER %.4e at theta/pi %.5f\n', ...
        peakGpsBer, theta_over_pi(peakIdx));
    fprintf('At that point: baseline %.4e, proposed %.4e\n', ...
        ber_baseline(peakIdx), ber_proposed(peakIdx));
end

function plot_results(results, outputDir)
    figure('Name', 'Case A MMSE phase neighborhood', 'Color', 'w');
    semilogy(results.theta_list / pi, results.ber_plot, 'LineWidth', 2);
    grid on;
    xlabel('\theta / \pi');
    ylabel('BER');
    title(sprintf('Matched MMSE BER near the dangerous phase, SNR=%g dB', ...
        results.config.snr_db));
    legend(cellstr(results.schemes), 'Location', 'best');
    saveas(gcf, fullfile(outputDir, 'caseA_mmse_phase_neighborhood.png'));
end

function value = percentile_by_sort(samples, probability)
    samples = sort(samples(:));
    index = max(1, min(numel(samples), ceil(probability * numel(samples))));
    value = samples(index);
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = defaultValue;
    end
end
