function results = run_papr_ccdf_comparison(numFrames, options)
%RUN_PAPR_CCDF_COMPARISON 鍙繍琛?PAPR CCDF 瀵规瘮锛屼笉鍋?BER sweep銆?%   姣旇緝 baseline / GPS(paper_grouping) / proposed_grouping銆?%
%   results = run_papr_ccdf_comparison(1000)

    rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    setup_paths(rootDir);

    if nargin < 1 || isempty(numFrames)
        numFrames = 1000;
    end
    if nargin < 2
        options = struct();
    end

    rngSeed = get_option(options, 'seed', 20260508);
    paprThresholds = get_option(options, 'papr_thresholds', 0:0.25:12);
    schemes = get_option(options, 'schemes', {'baseline', 'paper_grouping', 'proposed_grouping'});
    schemeLabels = get_option(options, 'scheme_labels', {'Baseline', 'GPS', 'Proposed'});

    baseConfig = afdm_config();
    baseConfig.channel.add_noise = false; %#ok<NASGU> PAPR 鍙湅鍙戝皠绔紝淇￠亾鏃犲叧銆?    numSchemes = numel(schemes);

    paprSamples = zeros(numFrames, numSchemes);
    selectedPatterns = cell(numFrames, numSchemes);

    fprintf('========== PAPR CCDF sampling ==========\n');
    fprintf('Frames=%d, schemes=%d\n', numFrames, numSchemes);

    for frameIdx = 1:numFrames
        rng(rngSeed + frameIdx, 'twister');
        txBits = randi([0, 1], baseConfig.waveform.NumSubcarriers * log2(baseConfig.modulation.M_mod), 1);

        for schemeIdx = 1:numSchemes
            cfg = apply_pre_chirp_scheme(baseConfig, schemes{schemeIdx});
            cfg.tx.bits = txBits;
            [~, papr, ~, txState] = afdm_tx_engine(cfg);
            paprSamples(frameIdx, schemeIdx) = papr;
            selectedPatterns{frameIdx, schemeIdx} = extract_pattern_label(txState.pre_chirp_profile);
        end

        if mod(frameIdx, 100) == 0 || frameIdx == numFrames
            fprintf('PAPR frame %5d/%d complete\n', frameIdx, numFrames);
        end
    end

    paprCcdf = compute_papr_ccdf(paprSamples, paprThresholds);
    summary = summarize_papr(paprSamples, schemeLabels);

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');

    results = struct();
    results.schemes = schemes;
    results.scheme_labels = schemeLabels;
    results.num_frames = numFrames;
    results.papr_thresholds = paprThresholds;
    results.papr_samples = paprSamples;
    results.papr_ccdf = paprCcdf;
    results.summary = summary;
    results.selected_patterns = selectedPatterns;
    results.timestamp = timestamp;

    plot_papr_ccdf_results(results, outputDir);

    matPath = fullfile(outputDir, ['papr_ccdf_comparison_' timestamp '.mat']);
    csvPath = fullfile(outputDir, ['papr_ccdf_summary_' timestamp '.csv']);
    save(matPath, 'results');
    writetable(summary, csvPath);
    fprintf('Saved PAPR CCDF results to %s\n', matPath);
    disp(summary);
end

function ccdf = compute_papr_ccdf(paprSamples, thresholds)
    ccdf = zeros(numel(thresholds), size(paprSamples, 2));
    for schemeIdx = 1:size(paprSamples, 2)
        for thresholdIdx = 1:numel(thresholds)
            ccdf(thresholdIdx, schemeIdx) = mean(paprSamples(:, schemeIdx) > thresholds(thresholdIdx));
        end
    end
end

function summary = summarize_papr(paprSamples, schemeLabels)
    method_name = string(schemeLabels(:));
    mean_papr = mean(paprSamples, 1).';
    median_papr = percentile_by_sort(paprSamples, 0.50).';
    p90_papr = percentile_by_sort(paprSamples, 0.90).';
    p99_papr = percentile_by_sort(paprSamples, 0.99).';
    max_papr = max(paprSamples, [], 1).';
    summary = table(method_name, mean_papr, median_papr, p90_papr, p99_papr, max_papr);
end

function values = percentile_by_sort(samples, q)
    sorted = sort(samples, 1);
    idx = max(1, min(size(sorted, 1), ceil(q * size(sorted, 1))));
    values = sorted(idx, :);
end

function plot_papr_ccdf_results(results, outputDir)
    figure('Name', 'PAPR CCDF comparison', 'Color', 'w');
    semilogy(results.papr_thresholds, results.papr_ccdf, 'LineWidth', 2);
    grid on;
    xlabel('PAPR threshold (dB)');
    ylabel('Pr(PAPR > threshold)');
    title(sprintf('PAPR CCDF comparison, %d frames', results.num_frames));
    legend(results.scheme_labels, 'Location', 'southwest');
    saveas(gcf, fullfile(outputDir, ['papr_ccdf_comparison_' results.timestamp '.png']));

    figure('Name', 'PAPR summary comparison', 'Color', 'w');
    bar([results.summary.mean_papr, results.summary.p90_papr, results.summary.p99_papr]);
    grid on;
    xticklabels(results.summary.method_name);
    ylabel('PAPR (dB)');
    title('PAPR summary');
    legend({'Mean', 'P90', 'P99'}, 'Location', 'northwest');
    saveas(gcf, fullfile(outputDir, ['papr_summary_comparison_' results.timestamp '.png']));
end

function label = extract_pattern_label(profile)
    if isfield(profile, 'selection') && isfield(profile.selection, 'selected_candidates')
        label = mat2str(profile.selection.selected_candidates(:).');
    elseif isfield(profile, 'scheme')
        label = profile.scheme;
    else
        label = '';
    end
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = defaultValue;
    end
end
