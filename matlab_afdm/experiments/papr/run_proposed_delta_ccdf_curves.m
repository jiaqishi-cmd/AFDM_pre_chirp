function results = run_proposed_delta_ccdf_curves(numFrames, options)
%RUN_PROPOSED_DELTA_CCDF_CURVES 姣旇緝 proposed 涓嶅悓 delta/c2 涓嬬殑 PAPR CCDF 鏇茬嚎銆?%   璇ュ疄楠屽洖褰掓甯搁殢鏈烘暟鎹抚锛屼笉浣跨敤 Case A fixed channel銆?
    rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
    addpath(rootDir);
    setup_paths(rootDir);

    if nargin < 1 || isempty(numFrames)
        numFrames = 1000;
    end
    if nargin < 2
        options = struct();
    end

    deltaRatioList = get_option(options, 'delta_ratio_list', [0.01 0.02 0.05 0.1 0.2 0.5 1.0]);
    paprThresholds = get_option(options, 'papr_thresholds', 0:0.05:12);
    seedBase = get_option(options, 'seed', 20260509);
    includeRefs = get_option(options, 'include_refs', true);

    baseConfig = afdm_config();
    N = baseConfig.waveform.NumSubcarriers;
    bitsPerFrame = N * log2(baseConfig.modulation.M_mod);
    c2Base = baseConfig.pre_chirp.base_c2;

    numDelta = numel(deltaRatioList);
    numRefs = 2 * includeRefs;
    numCurves = numRefs + numDelta;
    paprSamples = zeros(numFrames, numCurves);
    labels = strings(1, numCurves);

    fprintf('========== Proposed delta CCDF curves ==========\n');
    fprintf('Frames=%d, delta candidates=%d\n', numFrames, numDelta);

    for frameIdx = 1:numFrames
        rng(seedBase + frameIdx, 'twister');
        txBits = randi([0, 1], bitsPerFrame, 1);
        curveIdx = 0;

        if includeRefs
            curveIdx = curveIdx + 1;
            labels(curveIdx) = "Baseline";
            paprSamples(frameIdx, curveIdx) = frame_papr_for_scheme(baseConfig, 'baseline', txBits);

            curveIdx = curveIdx + 1;
            labels(curveIdx) = "GPS";
            paprSamples(frameIdx, curveIdx) = frame_papr_for_scheme(baseConfig, 'paper_grouping', txBits);
        end

        for deltaIdx = 1:numDelta
            cfg = baseConfig;
            cfg.pre_chirp.delta = deltaRatioList(deltaIdx) * c2Base;
            curveIdx = curveIdx + 1;
            labels(curveIdx) = sprintf("Proposed %.2g", deltaRatioList(deltaIdx));
            paprSamples(frameIdx, curveIdx) = frame_papr_for_scheme(cfg, 'proposed_grouping', txBits);
        end

        if mod(frameIdx, 100) == 0 || frameIdx == numFrames
            fprintf('PAPR frame %5d/%d complete\n', frameIdx, numFrames);
        end
    end

    paprCcdf = compute_papr_ccdf(paprSamples, paprThresholds);
    summary = summarize_papr(paprSamples, labels);

    outputDir = fullfile(fileparts(rootDir), 'results');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');

    results = struct();
    results.num_frames = numFrames;
    results.delta_ratio_list = deltaRatioList;
    results.labels = labels;
    results.papr_thresholds = paprThresholds;
    results.papr_samples = paprSamples;
    results.papr_ccdf = paprCcdf;
    results.summary = summary;
    results.timestamp = timestamp;

    plot_delta_ccdf(results, outputDir);
    save(fullfile(outputDir, ['proposed_delta_ccdf_curves_' timestamp '.mat']), 'results');
    writetable(summary, fullfile(outputDir, ['proposed_delta_ccdf_summary_' timestamp '.csv']));
    disp(summary);
end

function papr = frame_papr_for_scheme(baseConfig, scheme, txBits)
    cfg = apply_pre_chirp_scheme(baseConfig, scheme);
    cfg.tx.bits = txBits;
    [~, papr] = afdm_tx_engine(cfg);
end

function ccdf = compute_papr_ccdf(samples, thresholds)
    ccdf = zeros(numel(thresholds), size(samples, 2));
    for curveIdx = 1:size(samples, 2)
        for thresholdIdx = 1:numel(thresholds)
            ccdf(thresholdIdx, curveIdx) = mean(samples(:, curveIdx) > thresholds(thresholdIdx));
        end
    end
end

function summary = summarize_papr(samples, labels)
    method_name = labels(:);
    mean_papr = mean(samples, 1).';
    p90_papr = percentile_by_sort(samples, 0.90).';
    p99_papr = percentile_by_sort(samples, 0.99).';
    p999_papr = percentile_by_sort(samples, 0.999).';
    summary = table(method_name, mean_papr, p90_papr, p99_papr, p999_papr);
end

function values = percentile_by_sort(samples, q)
    sorted = sort(samples, 1);
    idx = max(1, min(size(sorted, 1), ceil(q * size(sorted, 1))));
    values = sorted(idx, :);
end

function plot_delta_ccdf(results, outputDir)
    figure('Name', 'Proposed delta PAPR CCDF curves', 'Color', 'w');
    semilogy(results.papr_thresholds, results.papr_ccdf, 'LineWidth', 1.8);
    grid on;
    xlabel('PAPR threshold (dB)');
    ylabel('Pr(PAPR > threshold)');
    title(sprintf('PAPR CCDF for proposed delta sweep, %d frames', results.num_frames));
    legend(cellstr(results.labels), 'Location', 'southwest');
    saveas(gcf, fullfile(outputDir, ['proposed_delta_ccdf_curves_' results.timestamp '.png']));
end

function value = get_option(options, name, defaultValue)
    if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
        value = options.(name);
    else
        value = defaultValue;
    end
end
