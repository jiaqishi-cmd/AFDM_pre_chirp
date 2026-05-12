% RUN_PARTIAL_REUSE_M_SWEEP
% 楠岃瘉涓嶅悓 AFDM block size M 涓嬶紝partial waveform reuse 鐩告瘮 full recompute
% 鐨?runtime 鍔犻€熻秼鍔裤€傝鑴氭湰涓嶄慨鏀瑰凡鏈変富浠跨湡鑴氭湰銆?
rootDir = find_afdm_root(fileparts(mfilename('fullpath')));
addpath(rootDir);
setup_paths(rootDir);

% ========================
% 鍙傛暟闆嗕腑璁剧疆
% ========================
rng(1, 'twister');
if ~exist('M_list', 'var'), M_list = [64 128 256]; end
if ~exist('V', 'var'), V = 8; end
if ~exist('W', 'var'), W = 3; end
if ~exist('delta_ratio', 'var'), delta_ratio = 0.1; end
if ~exist('beam_width', 'var'), beam_width = 8; end
if ~exist('topK', 'var'), topK = 8; end
if ~exist('search_os', 'var'), search_os = 2; end
if ~exist('final_os', 'var'), final_os = 4; end
if ~exist('numFrames', 'var'), numFrames = 200; end
if ~exist('selfCheckTrials', 'var'), selfCheckTrials = 20; end

outputDir = fullfile(fileparts(rootDir), 'results');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

numM = numel(M_list);
avg_runtime_full = nan(numM, 1);
avg_runtime_reuse = nan(numM, 1);
measured_speedup = nan(numM, 1);
avg_ifft_full = nan(numM, 1);
avg_ifft_reuse = nan(numM, 1);
avg_eval_count = nan(numM, 1);
max_papr_diff = nan(numM, 1);
pattern_agreement = nan(numM, 1);
selfcheck_max_rel_err = nan(numM, 1);
theory_speedup = nan(numM, 1);
C_full = nan(numM, 1);
C_reuse = nan(numM, 1);

fprintf('========== Partial reuse M sweep ==========\n');
fprintf('V=%d W=%d topK=%d searchOS=%d finalOS=%d frames=%d\n', ...
    V, W, topK, search_os, final_os, numFrames);

for mIdx = 1:numM
    M = M_list(mIdx);
    if mod(M, V) ~= 0
        warning('Skip M=%d because M must be divisible by V=%d.', M, V);
        continue;
    end

    base_c2 = sqrt(2) / (10 * M);
    delta = delta_ratio * base_c2;
    candidate_offsets = linspace(-delta, delta, W);
    group_index = repelem((1:V).', M / V);

    % Warm-up before timing.
    warmSymbols = 1 - 2 * randi([0, 1], M, 1);
    afdm.search.full_beam_search(warmSymbols, base_c2, candidate_offsets, group_index, search_os, final_os, beam_width, topK);
    afdm.search.reuse_beam_search(warmSymbols, base_c2, candidate_offsets, group_index, search_os, final_os, beam_width, topK);

    % Check that partial-waveform composition matches direct waveform.
    relErr = zeros(selfCheckTrials, 1);
    x = 1 - 2 * randi([0, 1], M, 1);
    sPart = afdm.search.precompute_partial_waveforms(x, base_c2, candidate_offsets, group_index, final_os);
    for trialIdx = 1:selfCheckTrials
        pattern = randi(W, 1, V);
        sDirect = afdm.search.direct_full_waveform(x, base_c2, candidate_offsets, group_index, pattern, final_os);
        sPartial = afdm.search.combine_partial_waveform(sPart, pattern);
        relErr(trialIdx) = norm(sDirect - sPartial) / max(norm(sDirect), eps);
    end
    selfcheck_max_rel_err(mIdx) = max(relErr);

    paprFull = zeros(numFrames, 1);
    paprReuse = zeros(numFrames, 1);
    patternFull = zeros(numFrames, V);
    patternReuse = zeros(numFrames, V);
    evalFull = zeros(numFrames, 1);
    evalReuse = zeros(numFrames, 1);
    ifftFull = zeros(numFrames, 1);
    ifftReuse = zeros(numFrames, 1);
    runtimeFull = zeros(numFrames, 1);
    runtimeReuse = zeros(numFrames, 1);

    for frameIdx = 1:numFrames
        rng(20260509 + 1000 * mIdx + frameIdx, 'twister');
        symbols = 1 - 2 * randi([0, 1], M, 1);

        tStart = tic;
        outFull = afdm.search.full_beam_search(symbols, base_c2, candidate_offsets, group_index, search_os, final_os, beam_width, topK);
        runtimeFull(frameIdx) = toc(tStart);

        tStart = tic;
        outReuse = afdm.search.reuse_beam_search(symbols, base_c2, candidate_offsets, group_index, search_os, final_os, beam_width, topK);
        runtimeReuse(frameIdx) = toc(tStart);

        paprFull(frameIdx) = outFull.papr;
        paprReuse(frameIdx) = outReuse.papr;
        patternFull(frameIdx, :) = outFull.pattern;
        patternReuse(frameIdx, :) = outReuse.pattern;
        evalFull(frameIdx) = outFull.eval_count;
        evalReuse(frameIdx) = outReuse.eval_count;
        ifftFull(frameIdx) = outFull.num_ifft;
        ifftReuse(frameIdx) = outReuse.num_ifft;
    end

    avg_runtime_full(mIdx) = mean(runtimeFull);
    avg_runtime_reuse(mIdx) = mean(runtimeReuse);
    measured_speedup(mIdx) = avg_runtime_full(mIdx) / avg_runtime_reuse(mIdx);
    avg_ifft_full(mIdx) = mean(ifftFull);
    avg_ifft_reuse(mIdx) = mean(ifftReuse);
    avg_eval_count(mIdx) = mean(evalFull);
    max_papr_diff(mIdx) = max(abs(paprFull - paprReuse));
    pattern_agreement(mIdx) = mean(all(patternFull == patternReuse, 2));

    N_eval_search = mean(evalFull - topK);
    C_full(mIdx) = N_eval_search * (search_os * M) * log2(search_os * M) + ...
        topK * (final_os * M) * log2(final_os * M);
    C_pre_search = V * W * (search_os * M) * log2(search_os * M);
    C_eval_search = N_eval_search * (search_os * M);
    C_pre_final = V * W * (final_os * M) * log2(final_os * M);
    C_refine = topK * V * (final_os * M);
    C_reuse(mIdx) = C_pre_search + C_eval_search + C_pre_final + C_refine;
    theory_speedup(mIdx) = C_full(mIdx) / C_reuse(mIdx);

    fprintf('M=%d: full %.3f ms/frame, reuse %.3f ms/frame, measured speedup %.3fx, theory %.3fx, IFFT %.1f/%.1f, pattern agreement %.2f%%\n', ...
        M, 1e3 * avg_runtime_full(mIdx), 1e3 * avg_runtime_reuse(mIdx), ...
        measured_speedup(mIdx), theory_speedup(mIdx), avg_ifft_full(mIdx), ...
        avg_ifft_reuse(mIdx), 100 * pattern_agreement(mIdx));

    save(fullfile(outputDir, 'results_partial_reuse_M_sweep.mat'), ...
        'M_list', 'V', 'W', 'delta_ratio', 'topK', 'beam_width', 'search_os', 'final_os', 'numFrames', ...
        'avg_runtime_full', 'avg_runtime_reuse', 'measured_speedup', 'theory_speedup', ...
        'avg_ifft_full', 'avg_ifft_reuse', 'avg_eval_count', 'max_papr_diff', ...
        'pattern_agreement', 'selfcheck_max_rel_err', 'C_full', 'C_reuse');
end

% ========================
% 缁樺浘
% ========================
fig1 = figure('Color', 'w');
plot(M_list, 1e3 * avg_runtime_full, '-o', 'LineWidth', 1.8); hold on;
plot(M_list, 1e3 * avg_runtime_reuse, '-s', 'LineWidth', 1.8);
grid on; xlabel('M'); ylabel('Runtime per frame (ms)');
title('Runtime versus AFDM block size');
legend({'Full recompute', 'Partial reuse'}, 'Location', 'northwest');
saveas(fig1, fullfile(outputDir, 'fig_partial_reuse_M_runtime.png'));

fig2 = figure('Color', 'w');
plot(M_list, measured_speedup, '-o', 'LineWidth', 1.8); hold on;
plot(M_list, theory_speedup, '-s', 'LineWidth', 1.8);
grid on; xlabel('M'); ylabel('Speedup');
title('Measured and theoretical speedup versus M');
legend({'Measured', 'Theory'}, 'Location', 'best');
saveas(fig2, fullfile(outputDir, 'fig_partial_reuse_M_speedup.png'));

fig3 = figure('Color', 'w');
plot(M_list, avg_ifft_full, '-o', 'LineWidth', 1.8); hold on;
plot(M_list, avg_ifft_reuse, '-s', 'LineWidth', 1.8);
grid on; xlabel('M'); ylabel('Average IFFT count per frame');
title('IFFT count versus AFDM block size');
legend({'Full recompute', 'Partial reuse'}, 'Location', 'best');
saveas(fig3, fullfile(outputDir, 'fig_partial_reuse_M_ifft_count.png'));

fprintf('Saved M sweep results to %s\n', fullfile(outputDir, 'results_partial_reuse_M_sweep.mat'));
