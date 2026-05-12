% RUN_PARTIAL_REUSE_THEORY_AND_TIMING
% 在实现大规模 partial waveform reuse 仿真前，先做理论复杂度计数和小规模 runtime 验证。
%
% proposed 是 group-wise c2 perturbation，每个 group 从 {c2-delta,c2,c2+delta}
% 中选择。partial reuse 利用：
%   s(pattern) = sum_v s_part{v, pattern(v)}
% 搜索时：
%   s_new = s_cur - s_part{v, old_w} + s_part{v, new_w}
% 避免每次候选评估都完整生成 IDAFT/IFFT 波形。

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
if ~exist('I', 'var'), I = 1; end
if ~exist('beam_width', 'var'), beam_width = 8; end
if ~exist('topK', 'var'), topK = 8; end
if ~exist('search_os', 'var'), search_os = 2; end
if ~exist('final_os', 'var'), final_os = 4; end
if ~exist('delta_ratio', 'var'), delta_ratio = 0.1; end
if ~exist('numFrames', 'var'), numFrames = 500; end

baseCfg = afdm_config();
base_c2 = baseCfg.pre_chirp.base_c2;
delta = delta_ratio * base_c2;
candidate_offsets = [-delta, 0, delta];
group_index = repelem((1:V).', M / V);

params = struct();
params.M = M;
params.V = V;
params.W = W;
params.I = I;
params.beam_width = beam_width;
params.topK = topK;
params.search_os = search_os;
params.final_os = final_os;
params.delta_ratio = delta_ratio;
params.numFrames = numFrames;
params.base_c2 = base_c2;

fprintf('========== Partial reuse theory and timing ==========\n');
fprintf('M=%d V=%d W=%d beam=%d topK=%d searchOS=%d finalOS=%d frames=%d\n', ...
    M, V, W, beam_width, topK, search_os, final_os, numFrames);

% ========================
% Warm-up，避免首次函数调用影响计时
% ========================
warmBits = randi([0, 1], M, 1);
warmSymbols = 1 - 2 * warmBits;
full_search_beam(warmSymbols, base_c2, candidate_offsets, group_index, search_os, final_os, beam_width, topK);
reuse_search_beam(warmSymbols, base_c2, candidate_offsets, group_index, search_os, final_os, beam_width, topK);

% ========================
% partial waveform 一致性自检
% ========================
rel_err_list = zeros(20, 1);
x = 1 - 2 * randi([0, 1], M, 1);
s_part = precompute_partial_waveforms(x, base_c2, candidate_offsets, group_index, final_os);
for idx = 1:numel(rel_err_list)
    pattern = randi(W, 1, V);
    s_direct = direct_full_waveform(x, base_c2, candidate_offsets, group_index, pattern, final_os);
    s_partial = combine_partial_waveform(s_part, pattern);
    rel_err_list(idx) = norm(s_direct - s_partial) / max(norm(s_direct), eps);
end
fprintf('Partial waveform self-check max rel err = %.3e\n', max(rel_err_list));

% ========================
% runtime 对比
% ========================
papr_full = zeros(numFrames, 1);
papr_reuse = zeros(numFrames, 1);
pattern_full = zeros(numFrames, V);
pattern_reuse = zeros(numFrames, V);
eval_count_full = zeros(numFrames, 1);
eval_count_reuse = zeros(numFrames, 1);
num_ifft_full = zeros(numFrames, 1);
num_ifft_reuse = zeros(numFrames, 1);
runtime_full = zeros(numFrames, 1);
runtime_reuse = zeros(numFrames, 1);

for frameIdx = 1:numFrames
    rng(20260509 + frameIdx, 'twister');
    bits = randi([0, 1], M, 1);
    symbols = 1 - 2 * bits;

    tStart = tic;
    outFull = full_search_beam(symbols, base_c2, candidate_offsets, group_index, search_os, final_os, beam_width, topK);
    runtime_full(frameIdx) = toc(tStart);

    tStart = tic;
    outReuse = reuse_search_beam(symbols, base_c2, candidate_offsets, group_index, search_os, final_os, beam_width, topK);
    runtime_reuse(frameIdx) = toc(tStart);

    papr_full(frameIdx) = outFull.papr;
    papr_reuse(frameIdx) = outReuse.papr;
    pattern_full(frameIdx, :) = outFull.pattern;
    pattern_reuse(frameIdx, :) = outReuse.pattern;
    eval_count_full(frameIdx) = outFull.eval_count;
    eval_count_reuse(frameIdx) = outReuse.eval_count;
    num_ifft_full(frameIdx) = outFull.num_ifft;
    num_ifft_reuse(frameIdx) = outReuse.num_ifft;
end

% ========================
% 理论复杂度计数
% ========================
N_eval_full = I * V * W;
N_eval_search = mean(eval_count_full - topK); % 去掉 final refine 的 topK 次精评
C_full_os4 = N_eval_full * (final_os * M) * log2(final_os * M);
C_full_os4 = N_eval_full * (final_os * M) * log2(final_os * M);
C_2stage_full = N_eval_search * (search_os * M) * log2(search_os * M) + ...
    topK * (final_os * M) * log2(final_os * M);
C_pre_search = V * W * (search_os * M) * log2(search_os * M);
C_eval_search = N_eval_search * (search_os * M);
C_pre_final = V * W * (final_os * M) * log2(final_os * M);
C_refine = topK * V * (final_os * M);
C_2stage_reuse = C_pre_search + C_eval_search + C_pre_final + C_refine;

complexity = struct();
complexity.N_eval_full = N_eval_full;
complexity.N_eval_search = N_eval_search;
complexity.C_full_os4 = C_full_os4;
complexity.C_2stage_full = C_2stage_full;
complexity.C_2stage_reuse = C_2stage_reuse;
complexity.speedup_theory = C_2stage_full / C_2stage_reuse;

fprintf('Average runtime full  = %.3f ms/frame\n', 1e3 * mean(runtime_full));
fprintf('Average runtime reuse = %.3f ms/frame\n', 1e3 * mean(runtime_reuse));
fprintf('Measured speedup = %.3f\n', mean(runtime_full) / mean(runtime_reuse));
fprintf('Theory speedup = %.3f\n', complexity.speedup_theory);
fprintf('Average eval count full/reuse = %.1f / %.1f\n', mean(eval_count_full), mean(eval_count_reuse));
fprintf('Average IFFT count full/reuse = %.1f / %.1f\n', mean(num_ifft_full), mean(num_ifft_reuse));
fprintf('Max |PAPR_full-PAPR_reuse| = %.3e dB\n', max(abs(papr_full - papr_reuse)));
fprintf('Pattern agreement = %.2f%%\n', 100 * mean(all(pattern_full == pattern_reuse, 2)));

outputDir = fullfile(fileparts(rootDir), 'results');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
save(fullfile(outputDir, 'results_partial_reuse_theory_and_timing.mat'), ...
    'params', 'complexity', 'runtime_full', 'runtime_reuse', ...
    'papr_full', 'papr_reuse', 'pattern_full', 'pattern_reuse', ...
    'eval_count_full', 'eval_count_reuse', 'num_ifft_full', 'num_ifft_reuse', ...
    'rel_err_list');

function out = full_search_beam(symbols, base_c2, offsets, group_index, search_os, final_os, beam_width, topK)
    V = max(group_index);
    W = numel(offsets);
    M = numel(symbols);
    states = make_initial_state(M, V, 2);
    eval_count = 0;
    num_ifft = 0;

    for groupId = 1:V
        expanded = repmat(make_state(M, V), 1, numel(states) * W);
        outIdx = 0;
        for stateIdx = 1:numel(states)
            for candId = 1:W
                outIdx = outIdx + 1;
                pattern = states(stateIdx).pattern;
                pattern(groupId) = candId;
                s = direct_full_waveform(symbols, base_c2, offsets, group_index, pattern, search_os);
                expanded(outIdx).pattern = pattern;
                expanded(outIdx).metric = compute_papr(s);
                eval_count = eval_count + 1;
                num_ifft = num_ifft + 1;
            end
        end
        [~, order] = sort([expanded.metric], 'ascend');
        states = expanded(order(1:min(beam_width, numel(order))));
    end

    finalK = min(topK, numel(states));
    finalPapr = zeros(1, finalK);
    for idx = 1:finalK
        s = direct_full_waveform(symbols, base_c2, offsets, group_index, states(idx).pattern, final_os);
        finalPapr(idx) = compute_papr(s);
        eval_count = eval_count + 1;
        num_ifft = num_ifft + 1;
    end
    [bestPapr, bestIdx] = min(finalPapr);
    out.papr = bestPapr;
    out.pattern = states(bestIdx).pattern;
    out.eval_count = eval_count;
    out.num_ifft = num_ifft;
end

function out = reuse_search_beam(symbols, base_c2, offsets, group_index, search_os, final_os, beam_width, topK)
    V = max(group_index);
    W = numel(offsets);
    M = numel(symbols);
    s_part_search = precompute_partial_waveforms(symbols, base_c2, offsets, group_index, search_os);
    num_ifft = V * W;
    initPattern = 2 * ones(1, V);
    initWaveform = combine_partial_waveform(s_part_search, initPattern);
    states = make_initial_state(M, V, 2);
    states.waveform = initWaveform;
    states.metric = compute_papr(initWaveform);
    eval_count = 0;

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
                eval_count = eval_count + 1;
            end
        end
        [~, order] = sort([expanded.metric], 'ascend');
        states = expanded(order(1:min(beam_width, numel(order))));
    end

    s_part_final = precompute_partial_waveforms(symbols, base_c2, offsets, group_index, final_os);
    num_ifft = num_ifft + V * W;
    finalK = min(topK, numel(states));
    finalPapr = zeros(1, finalK);
    for idx = 1:finalK
        s = combine_partial_waveform(s_part_final, states(idx).pattern);
        finalPapr(idx) = compute_papr(s);
        eval_count = eval_count + 1;
    end
    [bestPapr, bestIdx] = min(finalPapr);
    out.papr = bestPapr;
    out.pattern = states(bestIdx).pattern;
    out.eval_count = eval_count;
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

function s = direct_full_waveform(symbols, base_c2, offsets, group_index, pattern, os)
    M = numel(symbols);
    c2Vec = base_c2 * ones(M, 1);
    for groupId = 1:numel(pattern)
        c2Vec(group_index == groupId) = base_c2 + offsets(pattern(groupId));
    end
    s = full_waveform(symbols, c2Vec, os);
end

function s = partial_waveform_for_group(symbols, c2Vec, groupMask, os)
    M = numel(symbols);
    m = (0:M-1).';
    xPre = zeros(M, 1);
    xPre(groupMask) = symbols(groupMask) .* exp(1i * 2 * pi .* c2Vec(groupMask) .* (m(groupMask).^2));
    s = ifft_oversampled(xPre, os);
end

function s = full_waveform(symbols, c2Vec, os)
    M = numel(symbols);
    m = (0:M-1).';
    xPre = symbols(:) .* exp(1i * 2 * pi .* c2Vec(:) .* (m.^2));
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
