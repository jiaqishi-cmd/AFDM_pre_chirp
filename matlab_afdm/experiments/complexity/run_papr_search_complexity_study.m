% RUN_PAPR_SEARCH_COMPLEXITY_STUDY
% 鐮旂┒ proposed group-wise small c2 perturbation 鐨?PAPR 鎼滅储澶嶆潅搴︽姌涓€?% 涓嶄慨鏀圭幇鏈?PAPR/CCDF 涓昏剼鏈紱杩欓噷鐙珛瀹炵幇涓嶅悓 oversampling 鎼滅储绛栫暐銆?%
% 鎵€鏈夋渶缁堟姤鍛婄殑 PAPR 鍧囦娇鐢?OS=4 閲嶆柊璁＄畻锛屼繚璇佸叕骞炽€?% 鎼滅储闃舵鍙娇鐢?OS=1/2/4锛岀敤浜庣爺绌跺鏉傚害闄嶄綆銆?
rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
addpath(rootDir);
setup_paths(rootDir);

% ========================
% 鍙傛暟闆嗕腑璁剧疆
% ========================
rng(1, 'twister');
if ~exist('M', 'var'), M = 64; end
if ~exist('Lcpp', 'var'), Lcpp = 8; end %#ok<NASGU>
if ~exist('V', 'var'), V = 8; end
if ~exist('W', 'var'), W = 3; end
if ~exist('delta_ratio_list', 'var'), delta_ratio_list = [0.05 0.1 0.2]; end
if ~exist('numFrames', 'var'), numFrames = 3000; end
if ~exist('oversampling_final', 'var'), oversampling_final = 4; end
if ~exist('topK', 'var'), topK = 8; end

strategies = { ...
    struct('name', 'os4_full_greedy', 'search_os', 4, 'topK', 1); ...
    struct('name', 'os1_greedy_os4_refine', 'search_os', 1, 'topK', topK); ...
    struct('name', 'os2_greedy_os4_refine', 'search_os', 2, 'topK', topK) ...
    };

baseCfg = afdm_config();
base_c2 = baseCfg.pre_chirp.base_c2;
numStrategies = numel(strategies);
numDelta = numel(delta_ratio_list);

papr_all = zeros(numStrategies, numDelta, numFrames);
eval_count_all = zeros(numStrategies, numDelta, numFrames);
elapsed_time_total = zeros(numStrategies, numDelta);
selected_pattern_all = cell(numStrategies, numDelta, numFrames);

fprintf('========== PAPR search complexity study ==========\n');
fprintf('M=%d, V=%d, W=%d, frames=%d, topK=%d, final OS=%d\n', ...
    M, V, W, numFrames, topK, oversampling_final);

for deltaIdx = 1:numDelta
    deltaRatio = delta_ratio_list(deltaIdx);
    delta = deltaRatio * base_c2;
    [candidateSet, groupIndex] = build_proposed_candidates(M, V, base_c2, delta);

    fprintf('\nDelta/c2 = %.4g\n', deltaRatio);
    for strategyIdx = 1:numStrategies
        strategy = strategies{strategyIdx};
        fprintf('  Strategy %s\n', strategy.name);
        tic;
        for frameIdx = 1:numFrames
            rng(20260509 + 100000 * deltaIdx + frameIdx, 'twister');
            bits = randi([0, 1], M, 1);
            symbols = 1 - 2 * bits; % BPSK: 0->+1, 1->-1

            result = proposed_search_with_strategy( ...
                symbols, candidateSet, groupIndex, strategy.search_os, ...
                oversampling_final, strategy.topK);

            papr_all(strategyIdx, deltaIdx, frameIdx) = result.final_papr_os4;
            eval_count_all(strategyIdx, deltaIdx, frameIdx) = result.eval_count;
            selected_pattern_all{strategyIdx, deltaIdx, frameIdx} = result.selected_pattern;
        end
        elapsed_time_total(strategyIdx, deltaIdx) = toc;
    end
end

% ========================
% 鎸囨爣缁熻
% ========================
papr_at_1e2 = zeros(numStrategies, numDelta);
papr_at_1e3 = zeros(numStrategies, numDelta);
avg_eval_count = mean(eval_count_all, 3);
avg_runtime_per_frame = elapsed_time_total / numFrames;
loss_dB = zeros(numStrategies, numDelta);

for strategyIdx = 1:numStrategies
    for deltaIdx = 1:numDelta
        samples = squeeze(papr_all(strategyIdx, deltaIdx, :));
        papr_at_1e2(strategyIdx, deltaIdx) = percentile_by_sort(samples, 0.99);
        papr_at_1e3(strategyIdx, deltaIdx) = percentile_by_sort(samples, 0.999);
    end
end
for deltaIdx = 1:numDelta
    loss_dB(:, deltaIdx) = papr_at_1e3(:, deltaIdx) - papr_at_1e3(1, deltaIdx);
end

outputDir = fullfile(fileparts(rootDir), 'results');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

results = struct();
results.strategies = strategies;
results.strategy_names = string(cellfun(@(s) s.name, strategies, 'UniformOutput', false));
results.delta_ratio_list = delta_ratio_list;
results.papr_all = papr_all;
results.eval_count_all = eval_count_all;
results.elapsed_time_total = elapsed_time_total;
results.selected_pattern_all = selected_pattern_all;
results.papr_at_1e2 = papr_at_1e2;
results.papr_at_1e3 = papr_at_1e3;
results.avg_eval_count = avg_eval_count;
results.avg_runtime_per_frame = avg_runtime_per_frame;
results.loss_dB = loss_dB;
results.numFrames = numFrames;
results.topK = topK;
results.oversampling_final = oversampling_final;
save(fullfile(outputDir, 'results_papr_search_complexity_study.mat'), 'results');

plot_complexity_study(results, outputDir, timestamp);
print_summary(results);

% TODO: partial_waveform_reuse
% 鍙繘涓€姝ラ璁＄畻姣忎釜 group銆佹瘡涓€欓€夌殑 oversampled partial waveform:
% s_part{v,w}锛岀敤 s_new = s_old - s_part{v,old} + s_part{v,new}
% 閬垮厤姣忔鍊欓€夎瘎浼伴兘瀹屾暣璁＄畻 IDAFT/oversampled IFFT銆?
function [candidateSet, groupIndex] = build_proposed_candidates(N, V, baseC2, delta)
    groupIndex = repelem((1:V).', N / V);
    candidateSet = baseC2 + repmat([-delta, 0, delta], N, 1);
end

function result = proposed_search_with_strategy(symbols, candidateSet, groupIndex, searchOS, finalOS, keepK)
    numGroups = max(groupIndex);
    numCandidates = size(candidateSet, 2);
    N = numel(symbols);

    state.c2 = candidateSet(:, 2); % 浠庝腑蹇?c2 寮€濮?    state.pattern = 2 * ones(1, numGroups);
    state.metric = papr_for_c2(symbols, state.c2, searchOS);
    states = state;
    evalCount = 1;

    for groupId = 1:numGroups
        indices = groupIndex == groupId;
        expanded = repmat(make_empty_state(N, numGroups), 1, numel(states) * numCandidates);
        outIdx = 0;
        for stateIdx = 1:numel(states)
            for candId = 1:numCandidates
                outIdx = outIdx + 1;
                c2Trial = states(stateIdx).c2;
                c2Trial(indices) = candidateSet(indices, candId);
                patternTrial = states(stateIdx).pattern;
                patternTrial(groupId) = candId;
                expanded(outIdx).c2 = c2Trial;
                expanded(outIdx).pattern = patternTrial;
                expanded(outIdx).metric = papr_for_c2(symbols, c2Trial, searchOS);
                evalCount = evalCount + 1;
            end
        end
        [~, order] = sort([expanded.metric], 'ascend');
        keep = order(1:min(keepK, numel(order)));
        states = expanded(keep);
    end

    finalPapr = zeros(1, numel(states));
    for idx = 1:numel(states)
        finalPapr(idx) = papr_for_c2(symbols, states(idx).c2, finalOS);
        evalCount = evalCount + 1;
    end
    [bestPapr, bestIdx] = min(finalPapr);

    result.final_papr_os4 = bestPapr;
    result.selected_c2 = states(bestIdx).c2;
    result.selected_pattern = states(bestIdx).pattern;
    result.eval_count = evalCount;
end

function state = make_empty_state(N, numGroups)
    state = struct('c2', zeros(N, 1), 'pattern', zeros(1, numGroups), 'metric', Inf);
end

function papr = papr_for_c2(symbols, c2Vec, oversampling)
    signal = afdm_oversampled_waveform_for_papr(symbols, c2Vec, oversampling);
    papr = compute_papr(signal);
end

function signal = afdm_oversampled_waveform_for_papr(symbols, c2Vec, oversampling)
    % PAPR 鎼滅储鍙渶瑕佸彂灏勭骞呭害銆俻ost-chirp 涓哄崟浣嶆ā锛屼笉褰卞搷 PAPR銆?    % OS>1 鏃讹紝瀵?pre-chirped DAFT-domain 绗﹀彿鍋氶鍩熼浂濉厖杩戜技杩囬噰鏍枫€?    x = symbols(:);
    N = numel(x);
    n = (0:N-1).';
    if isscalar(c2Vec)
        c2Vec = c2Vec * ones(N, 1);
    else
        c2Vec = c2Vec(:);
    end
    xPre = x .* exp(1i * 2 * pi .* c2Vec .* (n.^2));

    if oversampling == 1
        signal = ifft(xPre) * sqrt(N);
        return;
    end

    Nos = oversampling * N;
    half = N / 2;
    Xos = zeros(Nos, 1);
    Xos(1:half) = xPre(1:half);
    Xos(end-half+1:end) = xPre(half+1:end);
    signal = ifft(Xos) * sqrt(Nos);
end

function value = percentile_by_sort(samples, q)
    sorted = sort(samples(:));
    idx = max(1, min(numel(sorted), ceil(q * numel(sorted))));
    value = sorted(idx);
end

function plot_complexity_study(results, outputDir, timestamp)
    names = cellstr(results.strategy_names);
    [~, deltaIdx01] = min(abs(results.delta_ratio_list - 0.1));
    thresholds = 0:0.05:12;

    figure('Name', 'PAPR CCDF strategies delta=0.1', 'Color', 'w');
    hold on;
    for strategyIdx = 1:numel(names)
        samples = squeeze(results.papr_all(strategyIdx, deltaIdx01, :));
        ccdf = arrayfun(@(x) mean(samples > x), thresholds);
        semilogy(thresholds, ccdf, 'LineWidth', 2);
    end
    grid on;
    xlabel('PAPR threshold (dB)');
    ylabel('Pr(PAPR > threshold)');
    title(sprintf('PAPR CCDF by search strategy, \\delta/c_2=%.3g', results.delta_ratio_list(deltaIdx01)));
    legend(names, 'Location', 'southwest');
    saveas(gcf, fullfile(outputDir, ['fig_search_strategy_ccdf_' timestamp '.png']));

    figure('Name', 'PAPR@1e-3 by strategy', 'Color', 'w');
    bar(results.papr_at_1e3.');
    grid on;
    xticklabels(string(results.delta_ratio_list));
    xlabel('\delta / c_2');
    ylabel('PAPR@CCDF=10^{-3} (dB)');
    title('PAPR tail by search strategy');
    legend(names, 'Location', 'best');
    saveas(gcf, fullfile(outputDir, ['fig_search_strategy_papr_point_' timestamp '.png']));

    figure('Name', 'Average eval count', 'Color', 'w');
    bar(results.avg_eval_count.');
    grid on;
    xticklabels(string(results.delta_ratio_list));
    xlabel('\delta / c_2');
    ylabel('Average candidate evaluations / frame');
    title('Search candidate evaluation count');
    legend(names, 'Location', 'best');
    saveas(gcf, fullfile(outputDir, ['fig_search_strategy_eval_count_' timestamp '.png']));

    figure('Name', 'Runtime per frame', 'Color', 'w');
    bar(results.avg_runtime_per_frame.');
    grid on;
    xticklabels(string(results.delta_ratio_list));
    xlabel('\delta / c_2');
    ylabel('Runtime / frame (s)');
    title('Search runtime per frame');
    legend(names, 'Location', 'best');
    saveas(gcf, fullfile(outputDir, ['fig_search_strategy_runtime_' timestamp '.png']));

    figure('Name', 'Complexity performance tradeoff', 'Color', 'w');
    hold on;
    markers = {'o', 's', '^'};
    for strategyIdx = 1:numel(names)
        scatter(results.avg_eval_count(strategyIdx, :), results.papr_at_1e3(strategyIdx, :), ...
            70, markers{min(strategyIdx, numel(markers))}, 'filled');
    end
    grid on;
    xlabel('Average candidate evaluations / frame');
    ylabel('PAPR@CCDF=10^{-3} (dB)');
    title('PAPR-complexity tradeoff');
    legend(names, 'Location', 'best');
    saveas(gcf, fullfile(outputDir, ['fig_search_complexity_tradeoff_' timestamp '.png']));
end

function print_summary(results)
    fprintf('\nPAPR search complexity summary\n');
    for deltaIdx = 1:numel(results.delta_ratio_list)
        fprintf('delta/c2=%.3g\n', results.delta_ratio_list(deltaIdx));
        for strategyIdx = 1:numel(results.strategy_names)
            fprintf('  %-24s PAPR@1e-3 %.3f dB | loss %.3f dB | eval %.1f | time %.4g s/frame\n', ...
                results.strategy_names(strategyIdx), ...
                results.papr_at_1e3(strategyIdx, deltaIdx), ...
                results.loss_dB(strategyIdx, deltaIdx), ...
                results.avg_eval_count(strategyIdx, deltaIdx), ...
                results.avg_runtime_per_frame(strategyIdx, deltaIdx));
        end
    end
end
