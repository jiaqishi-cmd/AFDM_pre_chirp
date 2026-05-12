function out = full_beam_search(symbols, baseC2, offsets, groupIndex, searchOs, finalOs, beamWidth, topK)
%FULL_BEAM_SEARCH Beam search that recomputes a full waveform per candidate.

    V = max(groupIndex);
    W = numel(offsets);
    M = numel(symbols);
    states = make_initial_state(M, V, ceil(W / 2));
    evalCount = 0;
    numIfft = 0;

    for groupId = 1:V
        expanded = repmat(make_state(M, V), 1, numel(states) * W);
        outIdx = 0;
        for stateIdx = 1:numel(states)
            for candId = 1:W
                outIdx = outIdx + 1;
                pattern = states(stateIdx).pattern;
                pattern(groupId) = candId;
                s = afdm.search.direct_full_waveform(symbols, baseC2, offsets, groupIndex, pattern, searchOs);
                expanded(outIdx).pattern = pattern;
                expanded(outIdx).metric = afdm.tx.compute_papr(s);
                evalCount = evalCount + 1;
                numIfft = numIfft + 1;
            end
        end
        [~, order] = sort([expanded.metric], 'ascend');
        states = expanded(order(1:min(beamWidth, numel(order))));
    end

    finalK = min(topK, numel(states));
    finalPapr = zeros(1, finalK);
    for idx = 1:finalK
        s = afdm.search.direct_full_waveform(symbols, baseC2, offsets, groupIndex, states(idx).pattern, finalOs);
        finalPapr(idx) = afdm.tx.compute_papr(s);
        evalCount = evalCount + 1;
        numIfft = numIfft + 1;
    end

    [bestPapr, bestIdx] = min(finalPapr);
    out.papr = bestPapr;
    out.pattern = states(bestIdx).pattern;
    out.eval_count = evalCount;
    out.eval_search = evalCount - finalK;
    out.eval_final = finalK;
    out.num_ifft = numIfft;
end

function state = make_initial_state(M, V, candId)
    state = make_state(M, V);
    state.pattern = candId * ones(1, V);
    state.metric = Inf;
end

function state = make_state(M, V)
    state = struct('pattern', zeros(1, V), 'metric', Inf, 'waveform', zeros(M, 1));
end
