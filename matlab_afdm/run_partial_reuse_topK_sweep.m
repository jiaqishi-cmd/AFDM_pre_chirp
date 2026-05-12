% RUN_PARTIAL_REUSE_TOPK_SWEEP
% 研究 partial waveform reuse 下 topK/beam 宽度变化带来的 PAPR 与复杂度折中。
% 最终 reported PAPR 均使用 final_os=4 重新计算。

rootDir = fileparts(mfilename('fullpath'));
addpath(rootDir);
setup_paths(rootDir);

% ========================
% 参数集中设置
% ========================
rng(1, 'twister');
if ~exist('M', 'var'), M = 64; end
if ~exist('V', 'var'), V = 8; end
if ~exist('W', 'var'), W = 3; end
if ~exist('delta_ratio_list', 'var'), delta_ratio_list = [0.1 0.2]; end
if ~exist('topK_list', 'var'), topK_list = [1 2 4 8 16]; end
if ~exist('search_os', 'var'), search_os = 2; end
if ~exist('final_os', 'var'), final_os = 4; end
if ~exist('numFrames', 'var'), numFrames = 1000; end
if ~exist('ccdfLevels', 'var'), ccdfLevels = [1e-2 1e-3]; end

outputDir = fullfile(fileparts(rootDir), 'results');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

numDelta = numel(delta_ratio_list);
numTopK = numel(topK_list);
papr_at_1e2 = nan(numDelta, numTopK);
papr_at_1e3 = nan(numDelta, numTopK);
avg_eval_count = nan(numDelta, numTopK);
avg_ifft_count = nan(numDelta, numTopK);
avg_runtime = nan(numDelta, numTopK);
avg_complexity_norm = nan(numDelta, numTopK);
papr_all = nan(numDelta, numTopK, numFrames);
eval_count_all = nan(numDelta, numTopK, numFrames);
eval_search_all = nan(numDelta, numTopK, numFrames);
eval_final_all = nan(numDelta, numTopK, numFrames);
ifft_count_all = nan(numDelta, numTopK, numFrames);
runtime_all = nan(numDelta, numTopK, numFrames);

if mod(M, V) ~= 0
    error('M=%d must be divisible by V=%d.', M, V);
end

base_c2 = sqrt(2) / (10 * M);
group_index = repelem((1:V).', M / V);

fprintf('========== Partial reuse topK sweep ==========\n');
fprintf('M=%d V=%d W=%d searchOS=%d finalOS=%d frames=%d\n', M, V, W, search_os, final_os, numFrames);

% Warm-up，避免第一次函数调用影响计时。
warmDelta = delta_ratio_list(1) * base_c2;
warmOffsets = linspace(-warmDelta, warmDelta, W);
warmSymbols = 1 - 2 * randi([0, 1], M, 1);
reuse_search_beam(warmSymbols, base_c2, warmOffsets, group_index, search_os, final_os, topK_list(1), topK_list(1));

for deltaIdx = 1:numDelta
    delta_ratio = delta_ratio_list(deltaIdx);
    delta = delta_ratio * base_c2;
    candidate_offsets = linspace(-delta, delta, W);

    for topIdx = 1:numTopK
        topK = topK_list(topIdx);
        beam_width = topK;

        paprVec = zeros(numFrames, 1);
        evalVec = zeros(numFrames, 1);
        evalSearchVec = zeros(numFrames, 1);
        evalFinalVec = zeros(numFrames, 1);
        ifftVec = zeros(numFrames, 1);
        runtimeVec = zeros(numFrames, 1);

        for frameIdx = 1:numFrames
            rng(20260509 + 100000 * deltaIdx + 1000 * topIdx + frameIdx, 'twister');
            symbols = 1 - 2 * randi([0, 1], M, 1);

            tStart = tic;
            out = reuse_search_beam(symbols, base_c2, candidate_offsets, group_index, search_os, final_os, beam_width, topK);
            runtimeVec(frameIdx) = toc(tStart);

            paprVec(frameIdx) = out.papr;
            evalVec(frameIdx) = out.eval_count;
            evalSearchVec(frameIdx) = out.eval_search;
            evalFinalVec(frameIdx) = out.eval_final;
            ifftVec(frameIdx) = out.num_ifft;
        end

        papr_all(deltaIdx, topIdx, :) = paprVec;
        eval_count_all(deltaIdx, topIdx, :) = evalVec;
        eval_search_all(deltaIdx, topIdx, :) = evalSearchVec;
        eval_final_all(deltaIdx, topIdx, :) = evalFinalVec;
        ifft_count_all(deltaIdx, topIdx, :) = ifftVec;
        runtime_all(deltaIdx, topIdx, :) = runtimeVec;

        papr_at_1e2(deltaIdx, topIdx) = papr_at_ccdf(paprVec, 1e-2);
        papr_at_1e3(deltaIdx, topIdx) = papr_at_ccdf(paprVec, 1e-3);
        avg_eval_count(deltaIdx, topIdx) = mean(evalVec);
        avg_ifft_count(deltaIdx, topIdx) = mean(ifftVec);
        avg_runtime(deltaIdx, topIdx) = mean(runtimeVec);
        avg_complexity_norm(deltaIdx, topIdx) = mean(evalSearchVec) * search_os + mean(evalFinalVec) * final_os;

        fprintf('delta/c2=%.3f topK=%d: PAPR@1e-3=%.3f dB, runtime=%.3f ms/frame, eval=%.1f\n', ...
            delta_ratio, topK, papr_at_1e3(deltaIdx, topIdx), 1e3 * avg_runtime(deltaIdx, topIdx), avg_eval_count(deltaIdx, topIdx));

        save(fullfile(outputDir, 'results_partial_reuse_topK_sweep.mat'), ...
            'M', 'V', 'W', 'delta_ratio_list', 'topK_list', 'search_os', 'final_os', 'numFrames', ...
            'papr_all', 'papr_at_1e2', 'papr_at_1e3', 'avg_eval_count', 'avg_ifft_count', ...
            'avg_runtime', 'avg_complexity_norm', 'eval_count_all', 'eval_search_all', ...
            'eval_final_all', 'ifft_count_all', 'runtime_all');
    end
end

% ========================
% Summary 推荐 topK
% ========================
for deltaIdx = 1:numDelta
    delta_ratio = delta_ratio_list(deltaIdx);
    [bestPapr, bestIdx] = min(papr_at_1e3(deltaIdx, :));
    recTopK = topK_list(bestIdx);
    % 由于有限帧数下 CCDF=1e-3 近似高分位/最大值，相邻 topK 可能有抖动。
    % 这里用“进入 best+0.1 dB 范围的最小 topK”作为饱和推荐，更稳健。
    nearBestIdx = find(papr_at_1e3(deltaIdx, :) <= bestPapr + 0.1, 1, 'first');
    satTopK = topK_list(nearBestIdx);
    fprintf('delta/c2=%.3f: best topK=%d (PAPR@1e-3=%.3f dB), saturation topK~%d\n', ...
        delta_ratio, recTopK, bestPapr, satTopK);
end

% ========================
% 绘图
% ========================
fig1 = figure('Color', 'w');
plot_by_delta(topK_list, papr_at_1e3, delta_ratio_list);
grid on; xlabel('topK'); ylabel('PAPR@CCDF=10^{-3} (dB)');
title('PAPR@10^{-3} versus topK');
saveas(fig1, fullfile(outputDir, 'fig_topK_papr_1e3.png'));

fig2 = figure('Color', 'w');
plot_by_delta(topK_list, papr_at_1e2, delta_ratio_list);
grid on; xlabel('topK'); ylabel('PAPR@CCDF=10^{-2} (dB)');
title('PAPR@10^{-2} versus topK');
saveas(fig2, fullfile(outputDir, 'fig_topK_papr_1e2.png'));

fig3 = figure('Color', 'w');
plot_by_delta(topK_list, 1e3 * avg_runtime, delta_ratio_list);
grid on; xlabel('topK'); ylabel('Runtime per frame (ms)');
title('Runtime versus topK');
saveas(fig3, fullfile(outputDir, 'fig_topK_runtime.png'));

fig4 = figure('Color', 'w');
plot_by_delta(topK_list, avg_eval_count, delta_ratio_list);
grid on; xlabel('topK'); ylabel('Average eval count per frame');
title('Candidate evaluations versus topK');
saveas(fig4, fullfile(outputDir, 'fig_topK_eval_count.png'));

fig5 = figure('Color', 'w'); hold on;
colors = lines(numDelta);
for deltaIdx = 1:numDelta
    scatter(1e3 * avg_runtime(deltaIdx, :), papr_at_1e3(deltaIdx, :), 60, colors(deltaIdx, :), 'filled');
    plot(1e3 * avg_runtime(deltaIdx, :), papr_at_1e3(deltaIdx, :), '-', 'Color', colors(deltaIdx, :), 'LineWidth', 1.4);
    for topIdx = 1:numTopK
        text(1e3 * avg_runtime(deltaIdx, topIdx), papr_at_1e3(deltaIdx, topIdx), ...
            sprintf('  K=%d', topK_list(topIdx)), 'Color', colors(deltaIdx, :));
    end
end
grid on; xlabel('Runtime per frame (ms)'); ylabel('PAPR@CCDF=10^{-3} (dB)');
title('Complexity-performance tradeoff');
legend(arrayfun(@(x) sprintf('\\delta/c_2=%.2f', x), delta_ratio_list, 'UniformOutput', false), 'Location', 'best');
saveas(fig5, fullfile(outputDir, 'fig_topK_complexity_tradeoff.png'));

fprintf('Saved topK sweep results to %s\n', fullfile(outputDir, 'results_partial_reuse_topK_sweep.mat'));

function out = reuse_search_beam(symbols, base_c2, offsets, group_index, search_os, final_os, beam_width, topK)
    V = max(group_index);
    W = numel(offsets);
    M = numel(symbols);
    s_part_search = precompute_partial_waveforms(symbols, base_c2, offsets, group_index, search_os);
    num_ifft = V * W;
    initPattern = ceil(W / 2) * ones(1, V);
    initWaveform = combine_partial_waveform(s_part_search, initPattern);
    states = make_initial_state(M, V, ceil(W / 2));
    states.waveform = initWaveform;
    states.metric = compute_papr(initWaveform);
    eval_search = 0;
    for groupId = 1:V
        expanded = repmat(make_state(M, V), 1, numel(states) * W);
        outIdx = 0;
        for stateIdx = 1:numel(states)
            oldCand = states(stateIdx).pattern(groupId);
            for candId = 1:W
                outIdx = outIdx + 1;
                pattern = states(stateIdx).pattern;
                pattern(groupId) = candId;
                waveform = states(stateIdx).waveform - s_part_search{groupId, oldCand} + s_part_search{groupId, candId};
                expanded(outIdx).pattern = pattern;
                expanded(outIdx).waveform = waveform;
                expanded(outIdx).metric = compute_papr(waveform);
                eval_search = eval_search + 1;
            end
        end
        [~, order] = sort([expanded.metric], 'ascend');
        states = expanded(order(1:min(beam_width, numel(order))));
    end
    s_part_final = precompute_partial_waveforms(symbols, base_c2, offsets, group_index, final_os);
    num_ifft = num_ifft + V * W;
    finalK = min(topK, numel(states));
    finalPapr = zeros(1, finalK);
    eval_final = 0;
    for idx = 1:finalK
        s = combine_partial_waveform(s_part_final, states(idx).pattern);
        finalPapr(idx) = compute_papr(s);
        eval_final = eval_final + 1;
    end
    [bestPapr, bestIdx] = min(finalPapr);
    out.papr = bestPapr;
    out.pattern = states(bestIdx).pattern;
    out.eval_search = eval_search;
    out.eval_final = eval_final;
    out.eval_count = eval_search + eval_final;
    out.num_ifft = num_ifft;
end

function state = make_initial_state(M, V, candId)
    state = make_state(M, V);
    state.pattern = candId * ones(1, V);
    state.metric = Inf;
end

function state = make_state(M, V)
    state = struct('pattern', zeros(1, V), 'metric', Inf, 'waveform', zeros(M, 1));
end

function s_part = precompute_partial_waveforms(symbols, base_c2, offsets, group_index, os)
    V = max(group_index);
    W = numel(offsets);
    M = numel(symbols);
    s_part = cell(V, W);
    for groupId = 1:V
        groupMask = group_index == groupId;
        for candId = 1:W
            c2Vec = base_c2 * ones(M, 1);
            c2Vec(groupMask) = base_c2 + offsets(candId);
            s_part{groupId, candId} = partial_waveform_for_group(symbols, c2Vec, groupMask, os);
        end
    end
end

function s = combine_partial_waveform(s_part, pattern)
    s = zeros(size(s_part{1, 1}));
    for groupId = 1:numel(pattern)
        s = s + s_part{groupId, pattern(groupId)};
    end
end

function s = partial_waveform_for_group(symbols, c2Vec, groupMask, os)
    M = numel(symbols);
    m = (0:M-1).';
    xPre = zeros(M, 1);
    xPre(groupMask) = symbols(groupMask) .* exp(1i * 2 * pi .* c2Vec(groupMask) .* (m(groupMask).^2));
    s = ifft_oversampled(xPre, os);
end

function s = ifft_oversampled(xPre, os)
    M = numel(xPre);
    if os == 1
        s = ifft(xPre) * sqrt(M);
        return;
    end
    Mos = os * M;
    half = M / 2;
    Xos = zeros(Mos, 1);
    Xos(1:half) = xPre(1:half);
    Xos(end-half+1:end) = xPre(half+1:end);
    s = ifft(Xos) * sqrt(Mos);
end

function value = papr_at_ccdf(paprVec, ccdfLevel)
    paprVec = sort(paprVec(:), 'ascend');
    n = numel(paprVec);
    idx = max(1, min(n, ceil((1 - ccdfLevel) * n)));
    value = paprVec(idx);
end

function plot_by_delta(x, yMat, deltaList)
    hold on;
    colors = lines(numel(deltaList));
    for idx = 1:numel(deltaList)
        plot(x, yMat(idx, :), '-o', 'LineWidth', 1.8, 'Color', colors(idx, :));
    end
    legend(arrayfun(@(x) sprintf('\\delta/c_2=%.2f', x), deltaList, 'UniformOutput', false), 'Location', 'best');
end
