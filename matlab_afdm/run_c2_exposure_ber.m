function results = run_c2_exposure_ber(diagnosticSource, scenarioIdx, caseRanks, snr_values, numFrames, receiverDetector)
%RUN_C2_EXPOSURE_BER Run BER curves for diagnostic top channel cases.
%   results = run_c2_exposure_ber()
%   results = run_c2_exposure_ber(diagnosticSource, scenarioIdx, caseRanks, snr_values, numFrames)
%
%   diagnosticSource can be a diagnostics struct or a path to a
%   c2_channel_diagnostics_*.mat file. If omitted, the newest diagnostics
%   file from ../results is loaded.

    rootDir = fileparts(mfilename('fullpath'));
    addpath(rootDir);
    setup_paths(rootDir);

    if nargin < 1 || isempty(diagnosticSource)
        diagnosticSource = newest_diagnostics_file(rootDir);
    end
    if nargin < 2 || isempty(scenarioIdx)
        scenarioIdx = 2;
    end
    if nargin < 3 || isempty(caseRanks)
        caseRanks = 1;
    end
    if nargin < 4 || isempty(snr_values)
        snr_values = [15, 20, 25, 30];
    end
    if nargin < 5 || isempty(numFrames)
        numFrames = 60;
    end
    if nargin < 6 || isempty(receiverDetector)
        receiverDetector = 'full_mmse';
    end

    diagnostics = load_diagnostics(diagnosticSource);
    scenario = diagnostics.scenarios{scenarioIdx};
    [~, order] = sort(scenario.gps_exposure_score, 'descend');

    schemes = {'baseline', 'paper_grouping', 'proposed_grouping'};
    labels = {'Baseline', 'GPS', 'Proposed'};
    base_config = afdm_config();
    base_config.receiver.detector = receiverDetector;
    base_seed = base_config.simulation.random_seed;
    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');

    results.diagnostic_source = diagnosticSource;
    results.scenario_idx = scenarioIdx;
    results.scenario_name = scenario.name;
    results.case_ranks = caseRanks;
    results.snr_values = snr_values;
    results.num_frames = numFrames;
    results.receiver_detector = receiverDetector;
    results.schemes = schemes;
    results.labels = labels;
    results.cases = cell(1, numel(caseRanks));
    results.ber = zeros(numel(snr_values), numel(schemes), numel(caseRanks));
    results.error_bits = zeros(numel(snr_values), numel(schemes), numel(caseRanks));
    results.total_bits = zeros(numel(snr_values), numel(schemes), numel(caseRanks));
    results.output_dir = outputDir;
    results.timestamp = timestamp;

    checkpointPath = fullfile(outputDir, ['c2_exposure_ber_checkpoint_' timestamp '.mat']);

    fprintf('========== c2 exposure BER ==========\n');
    fprintf('Scenario: %s\n', scenario.name);
    fprintf('Receiver detector: %s\n', receiverDetector);
    fprintf('Frames per SNR/scheme: %d\n\n', numFrames);

    for rank_pos = 1:numel(caseRanks)
        caseRank = caseRanks(rank_pos);
        caseIdx = order(caseRank);
        caseDef = scenario.cases{caseIdx};
        results.cases{rank_pos} = caseDef;

        fprintf('Case rank %d: %s | score %.3f\n', ...
            caseRank, caseDef.name, scenario.gps_exposure_score(caseIdx));
        fprintf('Diagnostic d2: baseline %.4g | GPS %.4g | proposed %.4g\n', ...
            scenario.two_symbol_distance(caseIdx, 1), ...
            scenario.two_symbol_distance(caseIdx, 2), ...
            scenario.two_symbol_distance(caseIdx, 3));

        for snr_idx = 1:numel(snr_values)
            for scheme_idx = 1:numel(schemes)
                cfg = apply_pre_chirp_scheme(base_config, schemes{scheme_idx});
                cfg = apply_case_to_config(cfg, caseDef);
                cfg.channel.multipath = true;
                cfg.channel.add_noise = true;
                cfg.channel.snr_db = snr_values(snr_idx);

                totalErr = 0;
                totalBits = 0;
                for frame_idx = 1:numFrames
                    frameSeed = base_seed + 1000000 * rank_pos + 100000 * snr_idx + frame_idx;
                    frame = simulate_frame(cfg, frameSeed);
                    totalErr = totalErr + frame.err_bits;
                    totalBits = totalBits + frame.total_bits;

                    if mod(frame_idx, 1000) == 0
                        fprintf('%-9s | case rank %d | SNR = %5.1f dB | frame %6d/%d | BER = %.3e\n', ...
                            labels{scheme_idx}, caseRank, snr_values(snr_idx), ...
                            frame_idx, numFrames, totalErr / totalBits);
                    end
                end

                results.error_bits(snr_idx, scheme_idx, rank_pos) = totalErr;
                results.total_bits(snr_idx, scheme_idx, rank_pos) = totalBits;
                results.ber(snr_idx, scheme_idx, rank_pos) = totalErr / totalBits;
                save(checkpointPath, 'results');

                fprintf('%-9s | SNR = %5.1f dB | errors = %5d / %-7d | BER = %.3e\n', ...
                    labels{scheme_idx}, snr_values(snr_idx), totalErr, totalBits, ...
                    results.ber(snr_idx, scheme_idx, rank_pos));
            end
        end
        fprintf('\n');
    end

    plot_exposure_ber(results);

    matPath = fullfile(outputDir, ['c2_exposure_ber_' timestamp '.mat']);
    save(matPath, 'results');
    fprintf('Saved exposure BER results to %s\n', matPath);
end

function filePath = newest_diagnostics_file(rootDir)
    outputDir = fullfile(fileparts(rootDir), 'results');
    files = dir(fullfile(outputDir, 'c2_channel_diagnostics_*.mat'));
    if isempty(files)
        error('No c2_channel_diagnostics_*.mat files found in %s.', outputDir);
    end
    [~, idx] = max([files.datenum]);
    filePath = fullfile(files(idx).folder, files(idx).name);
end

function diagnostics = load_diagnostics(source)
    if isstruct(source)
        diagnostics = source;
        return;
    end

    loaded = load(char(source), 'results');
    diagnostics = loaded.results;
end

function cfg = apply_case_to_config(cfg, caseDef)
    cfg.channel.profile = caseDef.name;
    cfg.channel.delay_taps = caseDef.delay_taps;
    cfg.channel.doppler_taps = caseDef.doppler_taps;
    cfg.channel.doppler_freq = caseDef.doppler_taps / cfg.waveform.NumSubcarriers;
    cfg.channel.chan_coef = caseDef.chan_coef / norm(caseDef.chan_coef);
end

function plot_exposure_ber(results)
    for case_pos = 1:numel(results.cases)
        safeName = regexprep(lower(results.cases{case_pos}.name), '[^a-z0-9]+', '_');
        safeName = regexprep(safeName, '^_|_$', '');

        figure('Name', ['C2 exposure BER - ' results.cases{case_pos}.name], 'Color', 'w');
        semilogy(results.snr_values, results.ber(:, :, case_pos), 'o-', 'LineWidth', 2);
        grid on;
        xlabel('SNR (dB)');
        ylabel('BER');
        title(sprintf('BER on diagnostic case: %s', results.cases{case_pos}.name), ...
            'Interpreter', 'none');
        legend(results.labels, 'Location', 'southwest');
        saveas(gcf, fullfile(results.output_dir, ...
            ['c2_exposure_ber_' safeName '_' results.timestamp '.png']));

        figure('Name', ['C2 exposure BER ratio - ' results.cases{case_pos}.name], 'Color', 'w');
        gpsBer = results.ber(:, 2, case_pos);
        proposedBer = results.ber(:, 3, case_pos);
        plot(results.snr_values, gpsBer ./ max(proposedBer, eps), 's-', 'LineWidth', 2);
        grid on;
        xlabel('SNR (dB)');
        ylabel('BER GPS / BER Proposed');
        title(sprintf('GPS vs proposed BER ratio: %s', results.cases{case_pos}.name), ...
            'Interpreter', 'none');
        yline(1, '--', 'Equal');
        saveas(gcf, fullfile(results.output_dir, ...
            ['c2_exposure_ber_ratio_' safeName '_' results.timestamp '.png']));
    end
end
